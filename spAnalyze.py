##############################################################################
#
# author: Aaron T Beaudoin
#
# Spectral line analysis tool to be run in iPython.
#
# Directory structure:
# spAnalyze.py/"Supernova name"/Data/:
#	breakpoints
#	data
#	eqw/
#		errors
#		MonteCarlo/Histograms
#		limits
#	Pictures/
#		hwr
#		norm_vs_res
# 	param
#	include
#	spectra
#
# Classes and functions:
#	spectrum(object)	initializes spectrum class
#	initSpec(file)		Retrieves initial spectrum properties
#	setGlobals(exclude, ions, epochs, order5, bins)		sets parameters used throughout
#	saveInc(spectrum)	Sets, saves, or loads values that determine the wavelength ranges to include in Equivalent Width measurements
#	loadInc(spectrum)
#	setInc(spectrum)	
#	interpolate(spectrum, breakpoints)	Helper function for breakpoint functions which calculates normalized spectrum
#	spNames(file, spectrum)		Records name, supernova, directory, and visit for spectrum
# 	savgol(flux, window_size, order, deriv = 0, rate = 1)	helper function to findErrors()
# 	findErrors(flux)	finds noise error in spectrum	
# 	getSpectrum(file, spectrum)		Gets spectrum from FLM file
# 	findNewBkpts(spectrum, linsp = False)		Finds and records new breakpoints and records normalized spectrum
#	getBkpts(spectrum, linsp = False)		retrieves previously calculated breakpoints from FLM files and records normalized spectrum
#	addBkpts(spectrum, linsp = False)		adds and records new and previously calculate breakpoints
#	rmvBkpts(spectrum, range_minimum, range_maximum)		removes breakpoints between a specified range
# 	saveBkpts(spectrum)		save breakpoints
#	randBreak(breakpoints, shift = 2.5)		helper function to monteCarlo(), randomizes breakpoints for variation in EQWs	
#	sigma3(spectrum, rest_wave, redshift = .01)		inds signal-to-noise 3 sigma limits for a given spectrum
# 	getLimits(spectrum, rest_wave, height, sigma, redshift = .01, plots = False)		Finds signal-to-noise limits for a given spectrum
# 	findEqw(spectrum, plot = False)		find EQWs for lines in spectrum sp	
# 	findErr(randomEqws, actualEqw, spectrum.name, spectrum.dir, ion)		helper function to monteCarlo() which calculates standard deviations	
# 	sepMgDoub(spectrum, fits = True)		Separates MgII doublets at 2800 angstroms from galaxy absorption and Milky Way absorption.
#	sepFeDoub(spectrum, rest_wave1, rest_wave2, indices_of_include_values = [], fit = False)		Separates FeII doublets at 2 rest wavelengths set on function call. 
# 	monteCarlo(numberForMC, spectrum, Mg_triplet = False, plot = False, Fe_triplet = False, rest_wavelengths = [], indices_of_include_values = [])		finds systematic errors for EQWs
# 	calcErr(spectrum)		calculates total error from monte carlo (systematic) and eqw error (statistical)
#	orgMeas(spectra, k, ionNum)		?
#	statusBar(count, total)		prints a status bar to run during loops	
#	saveParam(spectrum)		Saves and loads optimized Gaussian parameters
#	loadParam(spectrum)		
# 	saveData(spectrum)		saves all data to FLM files
# 	loadData(directory)		loads any saved data from FLM files. Requires a directory ("'Supernovae'/Data") ex. ('sn2011fe/Data'). Returns array of spectra with loaded data
#	chiSq(spectra)		Calculates the goodness of fit for each line via Chi-square tests
# 	plotEqw(spectra)		Plots eqws for each spectrum in parameter array. Each species plotted in same subplot.
# 	testPlots(spectra)		Plots Gaussian fits over a normalized continuum, the residual between the two, and parameter variation over time
#	scaleFlux(spectrum1, spectrum2)		Helper function that scales 2 continuum's with user input
#	scaledPlots(spectra)		Allows user to scale two spectra to analyse absorption features
#	proPlot(spectra)		Creates profile plots to double check variability in lines if detected in EQW measurements
#	pickSN()		Allows user to select the SN to work with (initializes initial parameters)
#
# Uses astropy, numpy, scipy, matplotlib, glob, sys, and math packages.
#############################################################################

from astropy.table import Table
import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from scipy import optimize
from scipy import integrate	
from scipy.stats import chi2
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import glob
import sys
import math
from math import factorial
import operator

#initializes spectrum class
class spectrum(object):
	"""Attributes can be added"""

#Retrieves initial spectrum properties
def initSpec(f):
	s = spectrum()
	spNames(f, s)
	getSpectrum(f, s)
	return s
	
#sets parameters used throughout
#param array ex- x values of lines to exclude in normalizing spectrum.
#param array i- List of ions to search for.
#param array ep- Epochs of each visit. Ex. {'visit1': 0.111, ...}
#param array ord- Orders of each visit for b-spline creation
#param array bin- Size of wavelength bins for b-spline creation
def setGlobals(ex, i, ep, ord, bin, sigI, lsp):
	global exclude
	global ions
	global epochs
	global order5
	global bins
	global sigIons
	global linsp
	exclude = ex
	ions = i
	epochs = ep
	order5 = ord
	bins = bin
	sigIons = sigI
	linsp = lsp

#Sets, saves, or loads values that determine the wavelength ranges to include in Equivalent Width measurements
def saveInc(sp):
	zeros = np.zeros(len(sp.inc))
	t = Table([sp.inc, zeros], names = ('col1', 'col2'))
	t.write(sp.dir + 'include/' + sp.name + '.dat', format = 'ascii')
	
def loadInc(sp):
	t = Table.read(sp.dir + 'include/' + sp.name + '.dat', format='ascii')
	sp.inc = t["col1"]
	
def setInc(sp):
	inc = input('Line ranges to include for '  + sp.name + ': ')
	sp.inc = inc
	
#Helper function for breakpoint functions which calculates normalized spectrum
def interpolate(sp, bp):
	bp.sort()
	x = sp.wave
	y = sp.flux
	
	#remove absorption lines from spectrum
	"""
	if (linsp == False):
		exclude = []
		"""
	
	if (linsp):
		lines = np.append([], np.where(x<exclude[0])[0])
		lines = [int(i) for i in lines]
		i = 1
		while (i<len(exclude)-1):
			lines = np.append(lines, np.where((x>exclude[i]) & (x<exclude[i+1]))[0])
			i += 2	
		lines = np.append(lines, np.where(x>exclude[len(exclude)-1])[0])
		
		#create b-spline over spectrum
		m = len(y)
		s = (m-  np.sqrt(2 * m), m + np.sqrt(2 * m))
		sm1 = inter.splrep(x[lines],y[lines], w = sp.var[lines], k = sp.order, s = s[1])
		
		xx = np.linspace(min(x), max(x), round(len(x) / sp.bin, 0))		
		
		linesxx = np.append([], np.where(xx<exclude[0])[0])
		linesxx = [int(i) for i in linesxx]
		i = 1
		while(i<len(exclude)-1):
			linesxx = np.append(linesxx, np.where((xx>exclude[i]) & (xx<exclude[i+1]))[0])
			i += 2	
		linesxx = np.append(linesxx, np.where(xx>exclude[len(exclude)-1])[0])
		xx = xx[linesxx]
		xx = np.append(xx,bp)
		
	else:
		m = len(y)
		s = (m-  np.sqrt(2 * m), m + np.sqrt(2 * m))
		sm1 = inter.splrep(x,y, w = sp.var, k = sp.order, s = s[1])
		xx = np.append([], np.linspace(min(x), 2000, round(len(np.where(x < 2000)[0]) / sp.bin, 0)))
		xx = np.append(xx, bp)
		xx = np.append(xx, np.linspace(3000, max(x), round(len(np.where(x > 3000)[0]) / sp.bin, 0)))
		xx = np.array(xx)
	
	""""
	if (rand):
		shift = 1.5
		if (sp.sn == '1992A'):
			shift = 1
		if ((sp.sn == '2011fe') & (sp.visit == 1)):
			shift = .25
		if (sp.sn == '2013dy'):
			shift = .25
		xx = [i + np.random.normal(0, .68, 1)[0] * shift for i in xx]
		indecies = []
		i = len(xx) - 10
		j = 0
		while (j < 10):
			if (xx[j] < min(x)):
				indecies.append(j)
			if (xx[i] > max(x)):
				indecies.append(i)
			i += 1
			j += 1
		xx = np.delete(xx, indecies)
	"""
	
	xx.sort()
	
	#fits b-spline over wavelength range
	y1 = inter.splev(xx,sm1)
	sm2 = inter.splrep(xx, y1, k = sp.order)
	y2 = inter.splev(x, sm2)
	
	#creates normalized spectrum
	norm = y / y2
	
	return y2, norm
	
#Records name, supernova, directory, and visit for spectrum
def spNames(file, sp):
	name = file[file.index('visit'):]
	name = name[:name.rindex('-')]
	
	sn = file[8:]
	sn = sn[sn.index('sn'):]
	sn = sn[2:]
	sn = sn[:sn.index('-')]
	
	dir = 'sn' + sn + '/Data/'
	
	visit = file[file.index('visit'):]
	visit = visit[5:]
	visit = visit[:visit.rindex('-')]
	
	sp.name = name
	sp.dir = dir
	sp.sn = sn
	try:
		sp.visit = float(visit)
	except ValueError:
		temp = visit[:1] + '.' + visit[-1:]
		sp.visit = float(temp)
	
	sp.epoch = epochs[sp.name]
	if sp.visit in order5:
		sp.order = 5
	else:
		sp.order = 3
	sp.bin = bins[sp.visit]
		
#helper function to findErrors()
def savgol(y, window_size, order, deriv=0, rate=1):
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')	
	
#finds noise error in spectrum	
def findErrors(y):
	win = round(len(y) / 10, 0)
	if win % 2 == 0:
		win += 1
	win = 11
	sf=savgol(y, win, 5)
	res=y-sf
	ares=np.abs(res)
	se = savgol(ares, win, 5)
	#norm_e = se/sflux
	return se

#Gets spectrum from FLM file	
def getSpectrum(file, sp):
	data = Table.read(file,format='ascii')
	sp.wave = data["col1"]
	sp.flux = data["col2"]

#Finds and records new breakpoints and records normalized spectrum
def findNewBkpts(sp):
	global linsp
	""""
	if (sp.visit == 1) & (sp.sn == '2013dy'):
		linsp = True
	"""
	if (sp.sn == '2014J') & (sp.visit < 4):
		linsp = True
	if (sp.sn == '2014J') & (sp.visit == 10):
		linsp = True
	x = sp.wave
	y = sp.flux
	sp.var = findErrors(y)**-1
	k = 1
	xmin1, xmax1, ymin1, ymax1 = 2000, 3000, 0, -1
	xmin2, xmax2, ymin2, ymax2 = 2000, 3000, 0, -1
	while(k == 1):
		points = input('Breakpoints for ' + sp.name + ': ')
		points = set(points)
		bps = [float(i) for i in points]
		sp_flux, norm = interpolate(sp, bps)
		fig1 = plt.figure(1)
		fig1.clf()
		plt.plot(x, y, color = 'k')
		plt.plot(x, sp_flux, color = 'b')
		ax1 = plt.gca()
		if (ymax1 > 0):
			ax1.set_xlim([xmin1, xmax1])
			ax1.set_ylim([ymin1, ymax1])
		plt.show()
		fig2 = plt.figure(2)
		fig2.clf()
		plt.plot(x, norm, linewidth = 2, color = 'k')
		ax2 = plt.gca()
		if (ymax2 > 0):
			ax2.set_xlim([xmin2, xmax2])
			ax2.set_ylim([ymin2, ymax2])
		#plt.ylim(0,1.5)
		plt.show()
		#plt.ion()
		t = True
		while(t):
			cont = raw_input('More breakpoints?: (y/n) ')
			xmin1, xmax1 = ax1.get_xlim()
			ymin1, ymax1 = ax1.get_ylim()
			xmin2, xmax2 = ax2.get_xlim()
			ymin2, ymax2 = ax2.get_ylim()
			if (plt.fignum_exists(1) == False):
				ymax1 = -1
			if (plt.fignum_exists(2) == False):
				ymax2 = -1
			if (cont == 'y'):
				t = False
			elif (cont == 'n'):
				sp.norm_flux = norm
				sp.sflux = sp_flux
				sp.bkpts = bps
				k = 2
				t = False
			else:
				continue
	
	t = Table([sp.bkpts], names = ('col1',))
	t.write('sn' + sp.sn + '/Data/breakpoints/' + sp.name + '.flm', format='ascii')

#retrieves previously calculated breakpoints from FLM files and records normalized spectrum
def getBkpts(sp):
	global linsp
	"""
	if (sp.visit == 1) & (sp.sn == '2013dy'):
		linsp = True
	"""
	if (sp.sn == '2014J') & (sp.visit < 4):
		linsp = True
	if (sp.sn == '2014J') & (sp.visit == 10):
		linsp = True
	t = Table.read('sn' + sp.sn + '/Data/breakpoints/' + sp.name + '.flm', format='ascii')
	sp.bkpts = t["col1"]
	sp.var = (findErrors(sp.flux)**-1)
	sp_flux, norm = interpolate(sp, sp.bkpts)
	sp.norm_flux = norm
	sp.sflux = sp_flux

#adds and records new and previously calculate breakpoints
def addBkpts(sp):
	global linsp
	"""
	if (sp.visit == 1) & (sp.sn == '2013dy'):
		linsp = True
	"""
	if (sp.sn == '2014J') & (sp.visit < 4):
		linsp = True
	if (sp.sn == '2014J') & (sp.visit == 10):
		linsp = True
	x = sp.wave
	y = sp.flux
	sp.var = findErrors(y)**-1
	k = 1
	xmin1, xmax1, ymin1, ymax1 = 2000, 3000, 0, -1
	xmin2, xmax2, ymin2, ymax2 = 2000, 3000, 0, -1
	sflux = 0
	while(k == 1):
		points = input('Breakpoints for ' + sp.name + ': ')
		points = set(points)
		bps = [float(i) for i in points]
		bps = np.append(bps, sp.bkpts)
		bps = set(bps)
		bps = [b for b in bps]
		sp_flux, norm = interpolate(sp, bps)
		fig1 = plt.figure(1)
		fig1.clf()
		plt.plot(x, y, color = 'k')
		plt.plot(x, sp_flux, color = 'b')
		ax1 = plt.gca()
		if (ymax1 > 0):
			ax1.set_xlim([xmin1, xmax1])
			ax1.set_ylim([ymin1, ymax1])
		plt.show()
		fig2 = plt.figure(2)
		fig2.clf()
		plt.plot(x, norm, color = 'k')
		ax2 = plt.gca()
		if (ymax2 > 0):
			ax2.set_xlim([xmin2, xmax2])
			ax2.set_ylim([ymin2, ymax2])
		#plt.ylim(0,1.5)
		plt.show()
		t = True
		while(t):
			cont = raw_input('More breakpoints?: (y/n) ')
			xmin1, xmax1 = ax1.get_xlim()
			ymin1, ymax1 = ax1.get_ylim()
			xmin2, xmax2 = ax2.get_xlim()
			ymin2, ymax2 = ax2.get_ylim()
			if (plt.fignum_exists(1) == False):
				ymax1 = -1
			if (plt.fignum_exists(2) == False):
				ymax2 = -1
			if (cont == 'y'):
				t = False
			elif (cont == 'n'):
				sp.norm_flux = norm
				sp.sflux = sp_flux
				sp.bkpts = bps
				k = 2
				t = False
			else:
				continue
	
	t = Table([sp.bkpts], names = ('col1',))
	t.write('sn' + sp.sn + '/Data/breakpoints/' + sp.name + '.flm', format='ascii')
	
#removes breakpoints between a specified range
def rmvBkpts(sp, min, max):
	if max < min:
		return 'ValueError: Max must be greater than min'
		
	for b in sp.bkpts:
		if b < 0:
			return 'ValueError: All breakpoints must be greater than or equal to 0'
		
	x = sp.wave
	y = sp.flux
	
	sp.bkpts.sort() 
	temp = []
	for b in sp.bkpts:
		if (b < min) | (b > max):
			temp.append(b)
						
	sp.bkpts = temp
	sp.bkpts.sort()
	sp.bkpts = np.array(sp.bkpts)
	
#saves breakpoints
def saveBkpts(sp):
	t = Table([sp.bkpts], names = ('col1',))
	t.write('sn' + sp.sn + '/Data/breakpoints/' + sp.name + '.flm', format = 'ascii')

#helper function to monteCarlo(), randomizes breakpoints for variation in EQWs	
def randBreak(breakpoints, shift = 2.5):
	bkpts = []
	for b in breakpoints:
		bkpts.append(b)
	bshift = np.random.randn(len(bkpts)) * shift
	new_bkpts = bshift + bkpts
	return new_bkpts

#Finds signal-to-noise 3 sigma limits for a given spectrum
def sigma3(sp, rest, red = .01):
	gaussian = lambda p, w, x: p[0]*np.exp(-(x - w*(1 + p[2]))**2 / (2*p[1]**2))
	
	gauss = gaussian([-.3, 3, red], rest, sp.wave)
	ind = np.where(gauss < -.001)
	
	res = 1 - sp.norm_flux[ind]
	temp = []
	for r in res:
		temp.append(r**2)
	rms = np.sqrt(np.sum(temp) / len(res))
	
	sp.sig3 = 9/np.sqrt(2) * rms
	
#Finds signal-to-noise limits for a given spectrum
def getLimits(sp, rest, height, sigma, red = .01, plots = False):
	gaussian = lambda p, w, x: p[0]*np.exp(-(x - w*(1 + p[2]))**2 / (2*p[1]**2))
	
	norm_flux = sp.norm_flux

	#sigma = np.linspace(2.4, 4, 17)
	#sigma = np.linspace(2.6, 4, 15)
	#height = np.linspace(-.02, -.4, 39)
	#height = np.linspace(-.03, -.3, 28)
	
	#sigma = np.linspace(3, 4, 11)
	#height = np.append(np.linspace(-.03, -.05, 5), np.linspace(-.06, -.3, 25))
	
	red = [red]
	
	meas = []
	guess = []
	err = []
	snr = []
	params = []
	for s in sigma:
		for h in height:
			for r in red:
				gauss = gaussian([h, s, r], rest, sp.wave)
				xindecies = np.where(gauss < -.001)
				xindecies = xindecies[0]
				xmin = sp.wave[min(xindecies)] - 1
				xmax = sp.wave[max(xindecies)] + 1
				flux = sp.norm_flux + gauss
				g = integrate.quad(lambda x: h * np.exp(-(x-rest*(1+r))**2/(2*s**2)), min(sp.wave), max(sp.wave))[0] * -1
				sp.norm_flux = flux
				sp.inc = [xmin, xmax]
				setGlobals(sp.inc, sp.inc, ions, epochs, sp.order, sp.bin, sigIons)
				if(plots):
					findEqw(sp, True)
					plt.xlim(2750, 3100)
					plt.ylim(0, 2)
					plt.savefig(sp.dir + 'eqw/limits/fits/' + sp.name + '/' + sp.name + '_' + str(h) + '.' + str(s) + '.' + str(r) + '.png')
					plt.close()
				else:
					findEqw(sp, False)
				ew = sp.eqws[0]
				er = sp.eqw_errs[0]
				sn = ew / er
				if ((g > .05) & (er != 0)):
					guess.append(g)
					meas.append(ew)
					err.append(er)
					snr.append(sn)
					params.append([h, s, r])
				
				sp.norm_flux = norm_flux
		
	font = {'family' : 'serif','color'  : 'black','weight' : 'bold','size' : 10,} 
	plt.figure()
	plt.scatter(guess, meas)
	plt.xlabel('EW In', fontdict = font)
	plt.ylabel('EW Out', fontdict = font)
	fit = np.polyfit(guess, meas, 1)
	x = np.linspace(-1, 7, 9)
	y = fit[0] * x + fit[1]
	plt.plot(x, y, linewidth = 3, color = 'r')
	plt.savefig(sp.dir + 'eqw/limits/' + sp.name + '_meas_vs_gauss' + str(r) + '.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
		
	fit = np.polyfit(guess, snr, 1)
	x = np.linspace(-1, 7, 9)
	y = fit[0] * x + fit[1]
	
	snr1 = (1 - fit[1]) / fit[0]
	snr3 = (3 - fit[1]) / fit[0]
	ew1 = 0
	ew3 = 0
	diff1 = 999999
	diff3 = 999999
	p1 = []
	p3 = []
	for i in range(len(snr)):
		temp = np.abs(1 - snr[i])
		if (temp < diff1):
			diff1 = temp
			ew1 = guess[i]
			#snr1 = snr[i]
			p1 = params[i]
		temp = np.abs(3 - snr[i])
		if (temp < diff3):
			diff3 = temp
			ew3 = guess[i]
			#snr3 = snr[i]
			p3 = params[i]
		
	plt.figure()
	plt.scatter(guess, snr)
	plt.plot(x, y, linewidth = 2.5, color = 'r')
	plt.text(4, 0, 'SNR 1: ' + str(round(snr1, 2)), fontdict = font)
	plt.text(7, 0, 'SNR 3: ' + str(round(snr3, 2)), fontdict = font)
	plt.ylabel('EW Out / EW Err (SNR)')
	plt.xlabel('EW In')
	#plt.xlim(-1,8)
	#plt.ylim(-1,8)
	plt.savefig(sp.dir + 'eqw/limits/' + sp.name + '_SNR_vs_gauss' + str(r) +'.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
			
	plt.figure()
	plt.title('SNR 1', fontdict = font)
	plt.plot(sp.wave, sp.norm_flux + gaussian(p1, rest, sp.wave), color = 'k', linestyle = 'steps', linewidth = 2)
	plt.xlim(2750, 2900)
	plt.ylim(0, 1.5)
	plt.xlabel('Wavelength')
	plt.ylabel('Normalized Flux')
	plt.savefig(sp.dir + 'eqw/limits/' + sp.name + '_SNR1' + str(r) + '.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
	plt.figure()
	plt.title('SNR 3', fontdict = font)
	plt.plot(sp.wave, sp.norm_flux + gaussian(p3, rest, sp.wave), color = 'k', linestyle = 'steps', linewidth = 2)
	plt.xlim(2750, 2900)
	plt.ylim(0, 1.5)
	plt.xlabel('Wavelength')
	plt.ylabel('Normalized Flux')
	plt.savefig(sp.dir + 'eqw/limits/' + sp.name + '_SNR3' + str(r) + '.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
#find EQWs for lines in spectrum sp	
def findEqw(sp, plot = False):
	include = sp.inc
	w = sp.wave
	f = sp.norm_flux
	
	gaussian = lambda p, x: 1+p[0]*np.exp(-(x-p[1])**2/(2*p[2]**2))
	residual = lambda p, x, y: y - gaussian(p, x)
	redshift = lambda obs, emit: (obs / emit) - 1
	
	if(plot):		
		start = np.where(w<include[0])
		end = np.where(w>include[len(include)-1])

		baseline_x = []
		baseline_y = []

		for s in range(len(start[0])):
			baseline_x.append(w[start[0][s]])
			baseline_y.append(1)

	i=0
	eqws = []
	eqw_errs = []
	param = []
	while(i<len(include)):
		inc = np.where((w>include[i]) & (w<include[i+1]))
		wave = w[inc]
		flux = f[inc]
		
		x = np.array(wave)
		y = np.array(flux)
		
		wobs = 0
		fmin = 24352345
		for j in range(len(x)):
			if (y[j] < fmin):
				fmin = y[j]
				wobs = x[j]
				
		width = .68 * (include[i + 1] - include[i])
						
		p = [fmin - 1, wobs, width]
		p1, pcov, infodict, errmsg, success= optimize.leastsq(residual, p[:], args=(x,y), full_output = 1)
		hix = np.linspace(min(x), max(x), len(x) * 10)
		gaus = gaussian(p1, hix)
		for j in p1:
			param.append(j)
		
		if(len(x) > len(p))	 and pcov is not None:
			s_sq = (residual(p1, x, y)**2).sum()
			pcov = pcov * s_sq / (len(x) - len(p))
		else:
			pcov = np.inf
		
		errors = []
		for j in range(len(p1)):
			try:
				errors.append(np.abs(pcov[j][j])**.5)
			except:
				errors.append(0.00)
			
		errors = np.array(errors)
		
		error = np.sqrt(2 * math.pi) * (errors[0]*p1[2] + errors[2]*p1[0])
		
		if(plot):
			for j in range(len(gaus)):
				baseline_x.append(hix[j])
				baseline_y.append(gaus[j])
			
			if(i<len(include)-2):
				ex = np.where((w>include[i+1]) & (w<include[i+2]))
				ex_wave = w[ex]
				
				for k in range(len(ex_wave)):
					baseline_x.append(ex_wave[k])
					baseline_y.append(1)
		
		eqw = integrate.quad(lambda x: p1[0]*np.exp(-(x-p1[1])**2/(2*p1[2]**2)), min(x), max(x))
		eqw = np.abs(eqw[0])
		
		eqws.append(eqw)
		eqw_errs.append(error)
		
		i+=2

	if(plot):
		for e in range(len(end[0])):
			baseline_x.append(w[end[0][e]])
			baseline_y.append(1)
		
		plt.figure()
		plt.plot(w,f,label='cont')
		plt.plot(baseline_x,baseline_y,label='gaus')
		legend = plt.legend(loc='upper right')
		plt.show()
	
	sp.eqws = eqws
	sp.eqw_errs = eqw_errs
	sp.param = param

#helper function to monteCarlo() which calculates standard deviations	
def findErr(eqws, original, name, dir, ion):
	eqws.sort()
	n = len(eqws)
	
	err_lo = eqws[int(round(.16*n))]
	err_hi = eqws[int(round(.84*n))]
	med = np.median(eqws)
	mean = np.mean(eqws)
	
	mc_lo = med - err_lo
	mc_hi = err_hi - med
	
	plt.figure()
	hist,bins=np.histogram(eqws,bins=len(eqws)/2)
	width=.7*(bins[1]-bins[0])
	center=(bins[:-1]+bins[1:])/2
	plt.bar(center,hist,align='center',width=width)
	plt.plot([err_lo],[2], 'go', label='err_lo')
	plt.plot([err_hi],[2],'co',label='err_hi')
	plt.plot([med],[2],'yo', label='median')
	plt.plot([mean], [2], 'mo', label = 'mean')
	plt.plot([original],[2],'ro', label='original')
	legend=plt.legend(loc='upper right')
	plt.savefig(dir + 'eqw/MonteCarlo/Histograms/MCHist' + name + '-' + ion + '.png',orientation='landscape')
	plt.show()
	#plt.close()
	
	return mc_lo, mc_hi
	
#Separates MgII doublets at 2800 angstroms from galaxy absorption and Milky Way absorption.
def sepMgDoub(sp, fit = True):
	if (type(sp.inc) != np.ndarray):
		sp.inc = np.array(sp.inc)
	ind = np.where((sp.inc > 2780) & (sp.inc < 2830))
	inc = sp.inc[ind[0]]
	rest1 = 2796.354
	rest2 = 2803.5311

	gaussian = lambda p, w, x: p[0]*np.exp(-(x - w*(1 + p[2]))**2 / (2*p[1]**2))
	doubGaus = lambda p, rest1, rest2, x: gaussian([p[0],p[1],p[2]], rest1, x) + gaussian([p[0],p[1],p[2]], rest2, x)
	combGaus = lambda p, rest1, rest2, x: doubGaus([p[0],p[1],p[2]], rest1, rest2, x) + doubGaus([p[3],p[4],p[5]], rest1, rest2, x)
	residual = lambda p, rest1, rest2, x, y: y - (1 + combGaus(p, rest1, rest2, x))

	redshift = lambda obs, emit: (obs / emit) - 1

	w = sp.wave
	f = sp.norm_flux

	inc1 = np.where((w>inc[0]) & (w<inc[1]))
	inc3 = np.where((w>inc[2]) & (w<inc[3]))
	inc2 = np.where((w>inc[4]) & (w<inc[5]))
	include = np.where((w>inc[0]) & (w<inc[5]))

	wave1 = w[inc1]
	flux1 = f[inc1]
	wave2 = w[inc2]
	flux2 = f[inc2]
	wave3 = w[inc3]
	flux3 = f[inc3]
	wave = w[include]
	flux = f[include]

	x1 = np.array(wave1)
	y1 = np.array(flux1)
	x2 = np.array(wave2)
	y2 = np.array(flux2)
	x3 = np.array(wave3)
	y3 = np.array(flux3)
	x = np.array(wave)
	y = np.array(flux)

	wobs1 = 0
	fmin1 = 134121324
	for i in range(len(x1)):
		if (y1[i] < fmin1):
			fmin1 = y1[i]
			wobs1 = x1[i]
						
	wobs2 = 0
	fmin2 = 12341234
	for i in range(len(x2)):
		if (y2[i] < fmin2):
			fmin2 = y2[i]
			wobs2 = x2[i]
			
	z1 = redshift(wobs1, rest1)
	z2 = redshift(wobs2, rest2)
	guess = [.68 * (inc[1] - inc[0]), .68 * (inc[len(inc) - 1] - inc[len(inc) - 2])]

	p = [fmin1 - 1, guess[0], z1, fmin2 - 1, guess[1], z2]

	p1, pcov, infodict, errmsg, success = optimize.leastsq(residual, p[:], args = (rest1, rest2, w, f), full_output = 1)

	if (fit):
		sp.wfit = np.linspace(min(w), max(w), len(w) * 10)
		sp.comb = combGaus(p1, rest1, rest2, sp.wfit) + 1
		sp.gausB = gaussian([p1[0], p1[1], p1[2]], rest2, sp.wfit) + 1
		sp.gausC = gaussian([p1[3], p1[4], p1[5]], rest1, sp.wfit) + 1 
		pa = [p1[0], p1[1], p1[2], p1[0], p1[1], p1[2]]
		pb = [p1[3], p1[4], p1[5], p1[3], p1[4], p1[5]]
		sp.doubA = doubGaus(pa, rest1, rest2, sp.wfit) + 1
		sp.doubB = doubGaus(pb, rest1, rest2, sp.wfit) + 1
		sp.sing = gaussian(pa[:3], rest1, sp.wfit) + gaussian(pb[:3], rest2, sp.wfit) + 1

	if (len(w) > len(p)) and pcov is not None:
		s_sq = (residual(p1, rest1, rest2, x3, y3)**2).sum()
		pcov = pcov * s_sq / (len(w) - len(p))
	else:
		pcov = np.inf

	errors = []
	for j in range(len(p1)):
		try:
			errors.append(np.abs(pcov[j][j])**.5)
		except:
			errors.append(0.00)
			
	errors = np.array(errors)

	error = []
	
	error.append(np.sqrt(2 * math.pi) * (errors[0]*p1[1] + errors[1]*p1[0]))
	error.append(np.sqrt(2 * math.pi) * (errors[0]*p1[1] + errors[1]*p1[0]))
	error.append(np.sqrt(2 * math.pi) * (errors[4]*p1[3] + errors[3]*p1[4]))
	error.append(np.sqrt(2 * math.pi) * (errors[4]*p1[3] + errors[3]*p1[4]))
	
	eqws = []
	eqws.append(np.abs(integrate.quad(lambda x: p1[0] * np.exp(-(x-rest1*(1+p1[2]))**2/(2*p1[1]**2)), min(x1), max(x2))[0]))
	eqws.append(np.abs(integrate.quad(lambda x: p1[0] * np.exp(-(x-rest2*(1+p1[2]))**2/(2*p1[1]**2)), min(x1), max(x2))[0]))
	eqws.append(np.abs(integrate.quad(lambda x: p1[3] * np.exp(-(x-rest1*(1+p1[5]))**2/(2*p1[4]**2)), min(x1), max(x2))[0]))
	eqws.append(np.abs(integrate.quad(lambda x: p1[3] * np.exp(-(x-rest2*(1+p1[5]))**2/(2*p1[4]**2)), min(x1), max(x2))[0]))

	oldEqws = sp.eqws
	newEqws = eqws
	newEqws = [float(n) for n in newEqws]
	if ('MgIb' in ions): 
		oldEqws[len(oldEqws) - 5: len(oldEqws) - 2] = newEqws
	else:
		oldEqws[len(oldEqws) - 4: len(oldEqws) - 1] = newEqws
	
	
	oldErr = sp.eqw_errs
	newErr = error	
	newErr = [float(n) for n in newErr]
	if ('MgIb' in ions):
		oldErr[len(oldErr) - 5: len(oldErr) - 2] = newErr
	else:
		oldErr[len(oldErr) - 4: len(oldErr) - 1] = newErr
	
	"""
	temp = []
	tempErr = []
	temp[
	add = True
	i = 0
	while(i < len(sp.inc)):
		if (sp.inc[i] > inc[0]) & (sp.inc[i] < inc[len(inc) - 1]):
			if add:
				temp = np.append(temp, newEqws)
				tempErr = np.append(tempErr, newErr)
				temp.append(float(newEqws[0]))
				temp.append(float(newEqws[1]))
				temp.append(float(newEqws[2]))
				temp.append(float(newEqws[3]))
				tempErr.append(float(newErr[0]))
				tempErr.append(float(newErr[0]))
				tempErr.append(float(newErr[1]))
				tempErr.append(float(newErr[1]))
				add = False
				i += 2
		else:
			temp = np.append(temp, [float(oldEqws[i])])
			tempErr = np.append(temp, [float(oldErr[i])])
			i += 2
	"""
			
	sp.eqws = oldEqws
	sp.eqw_errs = oldErr
	
	sp.param = p1
	
#Separates FeII doublets at 2 rest wavelengths set on function call. 
def sepFeDoub(sp, rest1, rest2, ind = [], fit = False):
	if (type(sp.inc) != np.ndarray):
		sp.inc = np.array(sp.inc)
	if (len(ind) == 0):
		ind = np.where((sp.inc > rest1 - 11) & (sp.inc < rest2 + 20))[0]
		inc = sp.inc[ind]
		ind = []
	elif (len(ind) != 0):
		inc = sp.inc[ind]
	else:
		return 'Include index values invalid'

	gaussian = lambda p, w, x: p[0]*np.exp(-(x - w*(1 + p[2]))**2 / (2*p[1]**2))
	doubGaus = lambda p, rest1, rest2, x: gaussian([p[0],p[1],p[2]], rest1, x) + gaussian([p[0],p[1],p[2]], rest2, x)
	combGaus = lambda p, rest1, rest2, x: doubGaus([p[0],p[1],p[2]], rest1, rest2, x) + doubGaus([p[3],p[4],p[5]], rest1, rest2, x)
	residual = lambda p, rest1, rest2, x, y: y - (1 + combGaus(p, rest1, rest2, x))

	redshift = lambda obs, emit: (obs / emit) - 1

	w = sp.wave
	f = sp.norm_flux

	inc1 = np.where((w>inc[0]) & (w<inc[1]))
	inc3 = np.where((w>inc[2]) & (w<inc[3]))
	inc2 = np.where((w>inc[4]) & (w<inc[5]))
	include = np.where((w>inc[0]) & (w<inc[5]))

	wave1 = w[inc1]
	flux1 = f[inc1]
	wave2 = w[inc2]
	flux2 = f[inc2]
	wave3 = w[inc3]
	flux3 = f[inc3]
	wave = w[include]
	flux = f[include]

	x1 = np.array(wave1)
	y1 = np.array(flux1)
	x2 = np.array(wave2)
	y2 = np.array(flux2)
	x3 = np.array(wave3)
	y3 = np.array(flux3)
	x = np.array(wave)
	y = np.array(flux)

	wobs1 = 0
	fmin1 = 134121324
	for i in range(len(x1)):
		if (y1[i] < fmin1):
			fmin1 = y1[i]
			wobs1 = x1[i]
						
	wobs2 = 0
	fmin2 = 12341234
	for i in range(len(x2)):
		if (y2[i] < fmin2):
			fmin2 = y2[i]
			wobs2 = x2[i]
			
	z1 = redshift(wobs1, rest1)
	z2 = redshift(wobs2, rest2)
	guess = [.68 * (inc[1] - inc[0]), .68 * (inc[len(inc) - 1] - inc[len(inc) - 2])]

	p = [fmin1 - 1, guess[0], z1, fmin2 - 1, guess[1], z2]

	p1, pcov, infodict, errmsg, success = optimize.leastsq(residual, p[:], args = (rest1, rest2, w, f), full_output = 1)
	
	if (fit):
		wfit = np.linspace(min(w), max(w), len(w) * 10)
		comb = combGaus(p1, rest1, rest2, wfit) + 1
		gausB = gaussian([p1[0], p1[1], p1[2]], rest2, wfit) + 1
		gausC = gaussian([p1[3], p1[4], p1[5]], rest1, wfit) + 1 
		pa = [p1[0], p1[1], p1[2], p1[0], p1[1], p1[2]]
		pb = [p1[3], p1[4], p1[5], p1[3], p1[4], p1[5]]
		doubA = doubGaus(pa, rest1, rest2, wfit) + 1
		doubB = doubGaus(pb, rest1, rest2, wfit) + 1
		sing = gaussian(pa[:3], rest1, wfit) + gaussian(pb[:3], rest2, wfit) + 1 
		fits = Table([wfit, comb, gausB, gausC, doubA, doubB, sing], names = ('wfit', 'comb', 'gausB', 'gausC', 'doubA', 'doubB', 'sing'))

	if (len(w) > len(p)) and pcov is not None:
		s_sq = (residual(p1, rest1, rest2, x3, y3)**2).sum()
		pcov = pcov * s_sq / (len(w) - len(p))
	else:
		pcov = np.inf

	errors = []
	for j in range(len(p1)):
		try:
			errors.append(np.abs(pcov[j][j])**.5)
		except:
			errors.append(0.00)
			
	errors = np.array(errors)

	error = []
	
	error.append(np.sqrt(2 * math.pi) * (errors[0]*p1[1] + errors[1]*p1[0]))	
	error.append(np.sqrt(2 * math.pi) * (errors[0]*p1[1] + errors[1]*p1[0]))
	error.append(np.sqrt(2 * math.pi) * (errors[4]*p1[3] + errors[3]*p1[4]))
	error.append(np.sqrt(2 * math.pi) * (errors[4]*p1[3] + errors[3]*p1[4]))
	
	eqws = []
	eqws.append(np.abs(integrate.quad(lambda x: p1[0] * np.exp(-(x-rest1*(1+p1[2]))**2/(2*p1[1]**2)), min(x1), max(x2))[0]))
	eqws.append(np.abs(integrate.quad(lambda x: p1[0] * np.exp(-(x-rest2*(1+p1[2]))**2/(2*p1[1]**2)), min(x1), max(x2))[0]))
	eqws.append(np.abs(integrate.quad(lambda x: p1[3] * np.exp(-(x-rest1*(1+p1[5]))**2/(2*p1[4]**2)), min(x1), max(x2))[0]))
	eqws.append(np.abs(integrate.quad(lambda x: p1[3] * np.exp(-(x-rest2*(1+p1[5]))**2/(2*p1[4]**2)), min(x1), max(x2))[0]))
	
	if (len(ind) == 0):
		inc = np.where((sp.inc > (rest1 - 11)) & (sp.inc < (rest2 + 20)))[0]
	elif (len(ind) != 0):
		inc = ind
	ind = []
	i = 0
	j = 0
	while(i < len(sp.inc)):
		if i in inc:
			ind.append(j)
		i += 2
		j += 1
		
	sp.eqws[ind[0]: ind[len(ind) - 1] + 1] = eqws
	sp.eqw_errs[ind[0]: ind[len(ind) - 1] + 1] = error
	
	if (fit):
		return fits
		
#finds systematic errors for EQWs
def monteCarlo(mc, sp, mgtriplet = False, fetriplet = False, rests = [], inds = [[]], plot = False):
	global linsp
	""""
	if (sp.visit == 1) & (sp.sn == '2013dy'):
		linsp = True
	"""
	if (sp.sn == '2014J') & (sp.visit < 4):
		linsp = True
	if (sp.sn == '2014J') & (sp.visit == 10):
		linsp = True
	shift = .5
	if (sp.sn == '2011by') & (sp.visit == 2):
		shift = .5
	if (sp.sn == '2012cg'):
		shift = 2
		
	tele = 'uv'
	if (sp.sn == '2013dy'):
		tele = 'hst'
		
	print "Finding error for", sp.name
	spectra=[]
	for i in range(mc):
		good = True
		nsp = initSpec('sn' + str(sp.sn) + '/Data/data/sn' + str(sp.sn) + '-visit' + str(int(sp.visit)) + '-' + tele + '.flm')
		nsp.bkpts = [bp + np.random.normal(0, .68, 1)[0] * shift for bp in sp.bkpts]
		nsp.var = findErrors(nsp.flux)**-1
		nsp.sflux, nsp.norm_flux = interpolate(nsp, nsp.bkpts)
		nsp.inc = sp.inc
		"""
		inc = []
		left, right = 0, 0
		j = 0
		while (j < len(nsp.inc)):
			ind = np.where((nsp.wave > nsp.inc[j] - 1.5) & (nsp.wave < nsp.inc[j + 1] + 1.5))[0]
			nf = nsp.norm_flux[ind]
			min, index = 999999, 0
			for f in range(len(nf)):
				if nf[f] < min:
					min, index = nf[f], f
			r, s = -1, 1
			rcont, scont = True, True
			rdiff, sdiff = 999, 999
			while (rcont | scont):
				if (index + r) >= 0:
					rf = nf[index + r]
				else:
					rcont = False
				try:
					sf = nf[index + s]
				except:
					scont = False
				if (np.abs(rf - 1) < rdiff) & rcont:
					left = nsp.wave[ind][index + r]
					rdiff = np.abs(rf - 1)
				if (np.abs(sf - 1) < sdiff) & scont:
					right = nsp.wave[ind][index + s]
					sdiff = np.abs(sf - 1)
				r -= 1
				s += 1
			cont = True
			while(cont):
				temp = np.where((nsp.wave > left) & (nsp.wave < right))[0]
				if len(temp) < 3:
					left -= .1
					right += .1
				else:
					cont = False
			inc = np.append(inc, [round(left, 2), round(right, 2)])
			j += 2
		inc = list(inc)
		for j in range(len(inc)):
			if inc.count(inc[j]) > 1:
				inc[j] -= .01
		inc = np.array(inc)
		nsp.inc = inc
		"""
		if (plot):
			findEqw(nsp, True)
			plt.xlim(2300, 2900)
			plt.ylim(0, 1.5)
		else:
			findEqw(nsp,False)
		if (mgtriplet):
			try:	
				sepMgDoub(nsp, False)
			except:
				good = False
			if (plot):
				plt.plot(sp.wfit, sp.comb, sp.wfit, sp.gausB, sp.wfit, sp.gausC, sp.wfit, sp.doubA, sp.wfit, sp.doubB, sp.wfit, sp.sing)
		if (fetriplet):
			l = len(rests) - 1
			m = len(inds) - 1
			while(l > 0):
				if (plot):
					fits = sepFeDoub(nsp, rests[l - 1], rests[l], inds[m], True)
					wfit = fits['wfit']
					plt.plot(wfit, fits['comb'], wfit, fits['gausB'], wfit, fits['gausC'], wfit, fits['doubA'], wfit, fits['doubB'], wfit, fits['sings'])
				else:
					sepFeDoub(nsp, rests[l - 1], rests[l], inds[m], False)
				l -= 2
				m -= 1
		if (plot):
			plt.savefig(sp.dir + 'eqw/MonteCarlo/lines/' + sp.name + ', rand ' + str(i + 1) + '.png')
			plt.close()
		if (0 in nsp.eqws):
			good = False
		if (good):
			spectra.append(nsp)
		else:
			i -= 1
		statusBar(i + 1, mc)
	sys.stdout.write('\n')
	
	
	mc_lo_errs = []
	mc_hi_errs = []
	for e in range(len(spectra[0].eqws)):
		temp = []
		for s in range(len(spectra)):
			temp.append(spectra[s].eqws[e])
		mc_lo, mc_hi = findErr(temp, sp.eqws[e], sp.name, sp.dir, ions[e])
		mc_lo_errs.append(mc_lo)
		mc_hi_errs.append(mc_hi)
	
	sp.mc_lo = mc_lo_errs
	sp.mc_hi = mc_hi_errs
	print "Errors found for", sp.name

#calculates total error from monte carlo (systematic) and eqw error (statistical)
def calcErr(sp):
	mc_lo = sp.mc_lo
	mc_hi = sp.mc_hi
	eqw_err = sp.eqw_errs
	
	err_lo = []
	err_hi = []
	for a in range(len(eqw_err)):
		err_lo.append((mc_lo[a]**2 + eqw_err[a]**2)**.5)
		err_hi.append((mc_hi[a]**2 + eqw_err[a]**2)**.5)
	
	sp.err_lo = err_lo
	sp.err_hi = err_hi
	
	t = Table([ions, sp.err_lo, sp.err_hi], names = ('col1','col2','col3'))
	t.write('sn' + sp.sn + '/Data/eqw/errors/' + sp.name + '.flm', format='ascii')

#?
def orgMeas(spectra, k, ionNum):
	s = spectra[k]
	i = 0
	j = 0
	rests = input('Rests: ')
	temp = []
	for r in rests:
		sigma3(s, r)
		temp.append(s.sig3)
	s.sig3 = temp
	ew, err, lo, hi = [], [], [], []
	for ion in ions:
		if ion in ionNum:
			ew.append(s.eqws[i])
			err.append(s.eqw_errs[i])
			lo.append(s.mc_lo[i])
			hi.append(s.mc_hi[i])
			i += 1
		else:
			ew.append(s.sig3[j])
			err.append(0)
			lo.append(0)
			hi.append(0)
			j += 1
	s.eqws = ew
	s.eqw_errs = err
	s.mc_lo = lo
	s.mc_hi = hi
	s.sig3 = np.zeros(len(s.eqws))
	spectra[k] = s

#prints a status bar to run during loops	
def statusBar(count, total, counter = ':'):
	sys.stdout.write('\r')
	sys.stdout.write('[')
	if ((float(count) / float(total)) >= .05):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .1):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .15):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .2):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .25):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .3):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .35):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .4):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .45):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .5):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .55):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .6):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .65):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .7):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .75):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .8):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .85):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .9):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if ((float(count) / float(total)) >= .95):
		sys.stdout.write(counter)
	else:
		sys.stdout.write(' ')
	if (count == total):
		sys.stdout.write(counter)
	per = str(int(round(float(count) / float(total), 2) * 100))
	sys.stdout.write('] ' + per + '%')
	
#Saves and loads optimized Gaussian parameters
def saveParam(sp):
	zeros = np.zeros(len(sp.param))
	t = Table([sp.param, zeros], names = ('col1', 'col2'))
	t.write(sp.dir + 'param/' + sp.name + '.dat', format = 'ascii')
	
def loadParam(sp):
	t = Table.read(sp.dir + 'param/' + sp.name + '.dat', format='ascii')
	sp.param = t["col1"]
	
#saves all data to FLM files
def saveData(sp):
	att = dir(sp)
	if('inc' in att):
		saveInc(sp)
	if('param' in att):
		saveParam(sp)
	#write eqws and errors
	attr = ()
	i = 18
	while (i<len(att)):
		if ('eqws' in att[i]) & ('eqws' not in attr):
			attr = attr + (att[i],)
		if ('eqw_errs' in att[i]) & ('eqw_errs' not in attr):
			attr = attr + (att[i],)
		if ('mc_lo' in att[i]) & ('mc_lo' not in attr):
			attr = attr + (att[i],)
		if ('mc_hi' in att[i]) & ('mc_hi' not in attr):
			attr = attr + (att[i],)
		if ('err_lo' in att[i]) & ('err_lo' not in attr):
			attr = attr + (att[i],)
		if ('err_hi' in att[i]) & ('err_hi' not in attr):
			attr = attr + (att[i],)
		if ('sig3' in att[i]) & ('sig3' not in attr):
			attr = attr + (att[i],)
		i += 1
	
	values = []
	for j in range(len(attr)):
		values.append(getattr(sp, attr[j]))
		
	t = Table(values, names = attr)
	t.write(sp.dir + 'eqw/' + sp.name + '.flm', format = 'ascii')
	
	#write spectrum to file
	attr = ()
	i = 18
	while (i<len(att)):
		if ('wave' in att[i]) & ('wave' not in attr):
			attr = attr + (att[i],)
		if ('flux' in att[i]) & ('flux' not in attr):
			attr = attr + (att[i],)
		if ('cont_err' in att[i]) & ('cont_err' not in attr):
			attr = attr + (att[i],)
		if ('var' in att[i]) & ('var' not in attr):
			attr = attr + (att[i],)
		if ('sflux' in att[i]) & ('sflux' not in attr):
			attr = attr + (att[i],)
		if ('norm_flux' in att[i]) & ('norm_flux' not in attr):
			attr = attr + (att[i],)
		i += 1
	
	values = []
	for j in range(len(attr)):
		values.append(getattr(sp, attr[j]))
		
	t = Table(values, names = attr)
	t.write(sp.dir + 'spectra/' + sp.name + '.flm', format = 'ascii')
	
	try:
		attr = ()
		i = 18
		while (i<len(att)):
			if ('wfit' in att[i]) & ('wfit' not in attr):
				attr = attr + (att[i],)
			if ('comb' in att[i]) & ('comb' not in attr):
				attr = attr + (att[i],)
			if ('doubA' in att[i]) & ('doubA' not in attr):
				attr = attr + (att[i],)
			if ('doubB' in att[i]) & ('doubB' not in attr):
				attr = attr + (att[i],)
			if ('sing' in att[i]) & ('sing' not in attr):
				attr = attr + (att[i],)
			i += 1
				
		values = []
		for j in range(len(attr)):
			values.append(getattr(sp, attr[j]))
		
		t = Table(values, names = attr)
		t.write(sp.dir + 'spectra/' + sp.name + '_fit.flm', format = 'ascii')
	except ValueError as e:
		print "No Fits"

#loads any saved data from FLM files. Requires a directory ("'Supernovae'/Data") ex. ('sn2011fe/Data')
#Returns array of spectra with loaded data
def loadData(dir):
	spectra = []
	
	#load spectrum info
	files = glob.glob(dir + '/data/*.flm')
	i = 1
	for f in files:
		sp = spectrum()
		spNames(f, sp)
		getSpectrum(f, sp)
		try:
			getBkpts(sp)
		except IOError as e:
			sys.stdout.write(" No Breakpoints")
		try:
			loadInc(sp)
		except IOError as e:
			sys.stdout.write(" No Include Values")
		try:
			loadParam(sp)
		except IOError as e:
			sys.stdout.write(" No Parameters Saved")
			
		file = sp.dir + 'spectra/' + sp.name + '.flm'
		t = Table.read(file, format = 'ascii')
		attr = t.colnames
		for a in attr:
			setattr(sp, a, t[a])	
		file = sp.dir + 'eqw/' + sp.name + '.flm'
		t = Table.read(file, format = 'ascii')
		attr = t.colnames
		for a in attr:
			setattr(sp, a, t[a])
		try:
			file = sp.dir + 'spectra/' + sp.name + '_fit.flm'
			t = Table.read(file, format = 'ascii')
			attr = t.colnames
			for a in attr:
				setattr(sp, a, t[a])
		except IOError as e:
			sys.stdout.write(" No fits 	")
		
			
		spectra.append(sp)
		statusBar(i, len(files))
		i += 1
	
	return spectra

#Calculates the goodness of fit for each line via Chi-square tests
def chiSq(spectra, order = 0):
	spectra.sort(key = lambda x: x.visit)
	
	eps = epochs.values()
	eps.sort()
	
	eqws = []
	for i in range(len(spectra[0].eqws)):
		ion = []
		for s in spectra:
			if (s.eqws[i] != 0):
				ion.append(s.eqws[i])
			else:
				ion.append(s.sig3[i])
		eqws.append(ion)
		
	errs = []
	for i in range(len(spectra[0].eqws)):
		err = []
		for s in spectra:
			err.append(np.average([s.err_lo[i], s.err_hi[i]]))
		errs.append(err)
		
	weights = []
	for e in errs:
		temp = []
		for i in e:
			if (i != 0):
				temp.append(1/i)
			else:
				temp.append(0)
		weights.append(temp)
		
	if (order != 0):
		fits = []
		for i in range(len(eqws)):
			fit = np.polyfit(eps, eqws[i], order, w = weights[i])
			fits.append(fit)
			
		exp = []
		for i in range(len(eqws)):
			e = []
			for j in range(len(eqws[i])):
				e.append(fits[i][0] * eps[j] + fits[i][1])
			exp.append(e)
			
		ddof = 1
	else:	
		exp = []
		count = 0
		for e in eqws:
			avg = np.average(e, weights = weights[count])
			ex = np.empty(len(e))
			ex[:] = avg
			exp.append(ex)
			count += 1
			
		ddof = 0
		
	chisq = []
	chisqRed = []
	pvalues = []
	for i in range(len(eqws)):
		temp = ddof
		chi = 0
		j = 0
		ews = []
		for o in eqws[i]:
			if (errs[i][j] != 0):
				chi += ((o - exp[i][j]) / errs[i][j])**2
			else:	
				temp += 1
			j += 1
		chisq.append(chi)
		chisqRed.append(chi/ (len(eqws[0]) - 1 - temp))
		pvalues.append(chi2.pdf(chi, (len(eqws[0]) - 1 - temp)))
		
	return chisq, chisqRed, pvalues
	
#Plots eqws for each spectrum in parameter array. Each species plotted in same subplot.
def plotEqw(spectra):
	#sorts spectra by epoch
	sort_sp = []
	temp = 0
	for s in spectra:
		if s.visit < temp:
			sort_sp.insert(len(sort_sp) - 1, s)
		else:
			sort_sp.append(s)
			temp = s.visit
	spectra = sort_sp
	
	epochs = []
	for s in spectra:
		epochs.append(s.epoch)

	#prepares EQWs to be plotted by ion
	eqws = []
	for i in range(len(spectra[0].eqws)):
		ion = []
		for s in spectra:
			ion.append(s.eqws[i])
		eqws.append(ion)
		
	hi_err = []
	lo_err = []
	for i in range(len(spectra[0].eqws)):
		hi = []
		lo = []
		for s in spectra:
			try:
				hi.append(s.err_hi[i])
				lo.append(s.err_lo[i])
			except:
				hi.append(0.00)
				lo.append(0.00)
		hi_err.append(hi)
		lo_err.append(lo)
		
	#finds number of plots to create based on species
	species = []
	for i in ions:
		temp = i[:2]
		#if any(temp in string for string in species):
		if temp in species:
			trash = 0
		else:
			species.append(temp)
	spec_num = []
	for s in species:
		num = 0
		for i in ions:
			if s in i:
				num += 1
		spec_num.append(num)
		
	weights = []
	for i in range(len(hi_err)):
		var = []
		for j in range(len(hi_err[i])):
			avg = (lo_err[i][j] + hi_err[i][j]) / 2
			var.append(1/avg)
		weights.append(var)	
		
	stdDev = []
	medErr = []
		
	font = {'family' : 'serif','color'  : 'black','weight' : 'bold','size' : 10,} 
	colors = ['black', 'blue', 'cyan', 'green', 'red', 'orange', 'yellow', 'purple', 'pink']
		
	fig = plt.figure(num = 5, dpi = 100, figsize = [10, (3*len(spec_num)) + 2], facecolor = 'w')
	gs = gridspec.GridSpec(1,1)
	plt.title(spectra[0].sn, fontdict = font)
	plt.gca().axes.get_xaxis().set_visible(False)
	plt.gca().axes.get_yaxis().set_visible(False)
	count = 0
	for j in range(len(species)):
		k = 0
		ax = fig.add_subplot(len(species), 1, j+1)
		plt.gca().set_color_cycle(colors)
		while (k < spec_num[j]):
			ax.scatter(epochs, eqws[count], color = colors[k])
			ax.errorbar(epochs, eqws[count], yerr = [lo_err[count], hi_err[count]], label = ions[count]) #, linestyle = "None")
			k += 1
			count += 1
		ax.text(.5, .9, species[j], horizontalalignment = 'center', transform = ax.transAxes, fontdict = font)
		plt.xlim(epochs[0] - 1, epochs[len(epochs) - 1] + 1)
		ax.set_xticks(epochs)
		ax.set_xticklabels([])
		plt.gca().axes.get_yaxis().set_visible(True)
		plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
		plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))
		plt.ylabel('EQW', fontdict = font)
		legend = plt.legend(bbox_to_anchor = (1, 1), loc = 2, shadow = True, prop = {'family':'serif'})
	ax.set_xticklabels(np.around(epochs, decimals = 1))
	plt.gca().axes.get_xaxis().set_visible(True)
	plt.xlabel('Days Since B Maximum', fontdict = font)
	plt.subplots_adjust(hspace = 0, right = .8)
	plt.savefig(spectra[0].dir + 'eqw/EQW by Epoch.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
	fig = plt.figure(num = 6, dpi = 100, figsize = [10, (3*len(spec_num)) + 2], facecolor = 'w')
	gs = gridspec.GridSpec(1,1)
	plt.title(spectra[0].sn, fontdict = font)
	plt.gca().axes.get_xaxis().set_visible(False)
	plt.gca().axes.get_yaxis().set_visible(False)
	count = 0
	for j in range(len(species)):
		k = 0
		ax = fig.add_subplot(len(species), 1, j+1)
		plt.gca().set_color_cycle(colors)
		while (k < spec_num[j]):
			ax.scatter(epochs, eqws[count] - np.average(eqws[count], weights = weights[count]), color = colors[k])
			ax.errorbar(epochs, eqws[count] - np.average(eqws[count], weights = weights[count]), yerr = [lo_err[count], hi_err[count]], label = ions[count]) #, linestyle = "None")
			k += 1
			count += 1
		ax.text(.5, .9, species[j], horizontalalignment = 'center', transform = ax.transAxes, fontdict = font)
		plt.xlim(epochs[0] - 1, epochs[len(epochs) - 1] + 1)
		ax.set_xticks(epochs)
		ax.set_xticklabels([])
		plt.gca().axes.get_yaxis().set_visible(True)
		plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
		plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))
		plt.ylabel('$\Delta$EQW', fontdict = font)
		legend = plt.legend(bbox_to_anchor = (1, 1), loc = 2, shadow = True, prop = {'family':'serif'})
	ax.set_xticklabels(np.around(epochs, decimals = 1))
	plt.gca().axes.get_xaxis().set_visible(True)
	plt.xlabel('Days Since B Maximum', fontdict = font)
	plt.subplots_adjust(hspace = 0, right = .8)
	plt.savefig(spectra[0].dir + 'eqw/delta EQW by Epoch.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
	fig = plt.figure(num = 7, dpi = 100, figsize = [10, 5], facecolor = 'w')
	gs = gridspec.GridSpec(1,1)
	stdDev = []
	medErr = []
	for j in range(len(eqws)):
		stdDev.append(np.std(eqws[j]))
		error = []
		for i in range(len(lo_err[j])):
			error.append(np.average([lo_err[j][i], hi_err[j][i]]))
		medErr.append(np.average(error))
	plt.title(spectra[0].sn + 'std dev vs. med Err', fontdict = font)
	for i in range(len(ions)):
		plt.scatter(medErr[i], stdDev[i], label = ions[i], color = colors[i])
	plt.ylabel('Std. Dev.', fontdict = font)
	plt.xlabel('Med. Err.', fontdict = font)
	legend = plt.legend(bbox_to_anchor = (1, 1), loc = 2, shadow = True, prop = {'family': 'serif'})
	plt.savefig(spectra[0].dir + 'eqw/stddev_vs_mederr.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
#Plots Gaussian fits over a normalized continuum, the residual between the two, and parameter variation over time
def testPlots(spectra):
	for s in spectra:
		sm1 = inter.splrep(s.wfit, s.comb)
		comb = inter.splev(s.wave, sm1)
		s.res = s.norm_flux - comb
		
	font = {'family' : 'serif','color'  : 'black','weight' : 'bold','size' : 10,} 
	ticks = np.linspace(.55, 1.15, 13)
	ticklabels = []
	for t in ticks:
		ticklabels.append(str(t))
	
	i = 1
	for s in spectra:
		fig = plt.figure(num = i, dpi = 100, figsize = [10,8], facecolor = 'w')
		#plt.title(s.name, fontdict = font)
		plt.gca().axes.get_xaxis().set_visible(False)
		plt.gca().axes.get_yaxis().set_visible(False)
		ax1 = fig.add_subplot(2, 1, 1)
		ax1.plot(s.wave, s.norm_flux, color = 'k', linestyle = 'steps', linewidth = 3)
		ind = np.where(s.comb < .99)
		ax1.plot(s.wfit[ind], s.comb[ind], color = 'r', label = 'Combined Fit', linewidth = 2, linestyle = '--', drawstyle = 'steps')
		ind = np.where(s.doubA < .99)
		ax1.plot(s.wfit[ind], s.doubA[ind], color = 'g', label = 'Double Gauss A & B', linewidth = 2, linestyle = '--', drawstyle = 'steps')
		ind = np.where(s.doubB < .99)
		ax1.plot(s.wfit[ind], s.doubB[ind], color = '#ff8000', label = 'Double Gauss C & D', linewidth = 2, linestyle = '--', drawstyle = 'steps')
		ind = np.where((s.wfit > 2792.8) & (s.wfit < 2802.85))
		ax1.plot(s.wfit[ind], s.sing[ind], color = 'b', label = 'Single Gauss A & D', linewidth = 2, linestyle = '--', drawstyle = 'steps')
		ind = np.where((s.wfit > 2811.13) & (s.wfit < 2821.09))
		ax1.plot(s.wfit[ind], s.sing[ind], color = 'b', linewidth = 2, linestyle = '--', drawstyle = 'steps')
		plt.xlim(2780, 2880)
		plt.ylim(.5, 1.15)
		ax1.set_xticklabels([])
		ax1.hlines(1, 2780, 2880, linestyles = 'dashed')
		ax1.set_yticks(ticks)
		ax1.set_yticklabels(ticklabels)
		plt.gca().axes.get_yaxis().set_visible(True)
		plt.ylabel('Normalized Flux', fontdict = font)
		plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))
		
		ax2 = fig.add_subplot(2, 1, 2)
		ax2.plot(s.wave, s.res, color = 'k', linestyle = 'steps', linewidth = 3)
		ax2.hlines(0, 2780, 2880, linestyles = 'dashed')
		plt.xlim(2780, 2880)
		plt.ylim(-.2, .2)
		plt.gca().axes.get_yaxis().set_visible(True)
		plt.ylabel('Residual', fontdict = font)

		plt.gca().axes.get_yaxis().set_visible(True)
		plt.xlabel('Wavelength', fontdict = font)
		plt.subplots_adjust(hspace = 0, right = .8)
		plt.savefig(s.dir + 'Pictures/norm_vs_res/norm_vs_res_' + str(i) + '.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
		i += 1
		
	heightMW = []
	sigmaMW = []
	redMW = []
	heightHG = []
	sigmaHG = []
	redHG = []
	epochs = []
	for s in spectra:
		p1 = s.param
		heightMW.append(p1[0])
		heightHG.append(p1[3])
		sigmaMW.append(p1[1])
		sigmaHG.append(p1[4])
		redMW.append(p1[2])
		redHG.append(p1[5])
		epochs.append(s.epoch)
		
	fig = plt.figure(num = i, dpi = 100, figsize = [10,8], facecolor = 'w')
	plt.title(s.sn, fontdict = font)
	plt.gca().axes.get_xaxis().set_visible(False)
	plt.gca().axes.get_yaxis().set_visible(False)
	ax1 = fig.add_subplot(3, 1, 1)
	ax1.scatter(epochs, heightMW - np.average(heightMW), color = 'k', label = 'Milky Way')
	ax1.scatter(epochs, heightHG - np.average(heightHG), color = 'c', label = 'Host Galaxy')
	legend = plt.legend(bbox_to_anchor = (1, 1), loc = 2, shadow = True, prop = {'family':'serif'})
	ax1.set_xticklabels([])
	plt.gca().axes.get_yaxis().set_visible(True)
	plt.ylabel('Heights', fontdict = font)
	plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))

	ax2 = fig.add_subplot(3, 1, 2)
	ax2.scatter(epochs, sigmaMW - np.average(sigmaMW), color = 'b', label = 'Milky Way')
	ax2.scatter(epochs, sigmaHG - np.average(sigmaHG), color = 'g', label = 'Host Galaxy')
	legend = plt.legend(bbox_to_anchor = (1, 1), loc = 2, shadow = True, prop = {'family':'serif'})
	ax2.set_xticklabels([])
	plt.gca().axes.get_yaxis().set_visible(True)
	plt.ylabel('Widths', fontdict = font)
	plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))
	
	ax3 = fig.add_subplot(3, 1, 3)
	ax3.scatter(epochs, redMW - np.average(redMW), color = 'g', label = 'Milky Way')
	ax3.scatter(epochs, redHG - np.average(redHG), color = '#ff8000', label = 'Host Galaxy')
	legend = plt.legend(bbox_to_anchor = (1, 1), loc = 2, shadow = True, prop = {'family':'serif'})
	plt.ylabel('Redshifts', fontdict = font)

	plt.gca().axes.get_yaxis().set_visible(True)
	plt.xlabel('Epochs', fontdict = font)
	plt.subplots_adjust(hspace = 0, right = .8)
	plt.savefig(s.dir + 'Pictures/hwr/hwr_' + s.sn + '.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	i += 1

#Helper function that scales 2 continuum's with user input
def scaleFlux(s1, s2):
	print s1.name + '-' + s2.name
	plt.figure()
	plt.plot(s1.wave, s1.flux)
	plt.plot(s2.wave, s2.flux)
	cont = True
	while(cont):
		scale = raw_input('Input scale factor: ')
		scale = float(scale)
		for i in range(len(s1.flux)):
			s1.flux[i] *= scale
		plt.figure()
		plt.plot(s1.wave, s1.flux)
		plt.plot(s2.wave, s2.flux)
		c = True
		while(c):
			temp = raw_input('Try another scale factor? (y/n): ')
			if (temp == 'y'):
				c = False
				cont = True
			elif(temp == 'n'):
				c = False
				cont = False
			else:
				c = True
				cont = True

#Allows user to scale two spectra to analyse absorption features
def scaledPlots(spectra):
	i = 0
	while(i < len(spectra) - 1):
		s = spectra[i]
		sp = spectra[i + 1]
		scaleFlux(s, sp)
		i += 1
	
#Creates profile plots to double check variability in lines if detected in EQW measurements
def proPlot(spectra, num):
	sort_sp = []
	temp = 0
	for s in spectra:
		if s.visit < temp:
			sort_sp.insert(len(sort_sp) - 1, s)
		else:
			sort_sp.append(s)
			temp = s.visit
	spectra = sort_sp
	
	i = 0
	while i < len(spectra):
		spec = []
		for n in range(num):
			try:
				spec.append(spectra[i + n])
			except:
				continue
		t = spec[0].name + '-' + spec[len(spec) - 1].name
		print t
		j = 0
		while (j < len(spec[0].inc) - 1):
			plt.figure()
			for s in spec:
				plt.plot(s.wave, s.norm_flux, label = s.name, linestyle = 'steps', linewidth = 2)
			plt.xlim(spec[0].inc[j] - 20, spec[0].inc[j+1] + 20)
			
			ymin, ymax = 9999999, -9999999
			for s in spec:
				ind = np.where((s.wave > s.inc[j] - 20) & (s.wave < s.inc[j + 1]))[0]
				tmax = max(s.norm_flux[ind])
				tmin = min(s.norm_flux[ind])
				if tmin < ymin:
					ymin = tmin
				if tmax > ymax:
					ymax = tmax
			ymax += .2
			ymin -= .2
			plt.ylim(ymin, ymax)
			#name = ions[j]
			name = raw_input('Absorption Feature: ')
			plt.title(t + ': ' + name)
			plt.show()
			legend = plt.legend(loc = 'upper right')
			plt.savefig('sn' + spec[0].sn + '/Data/Pictures/' + name	+ '-' + t + '.png')
			plt.close()
			j += 2
		i += num - 1
		
#Allows user to select the SN to work with (initializes initial parameters)
def pickSN(sn = ''):
	cont = True

	while(cont):
		if (sn == ''):
			sn = raw_input('Choose Supernovae (92a, 11by, 11ek, 11fe, 11iv, 12cg, 13dy, 14j): ')

		if (sn == '92a' or sn == '92A'):
			#92a params
			exclude = [2795.6, 2815.67]
			ions = []
			sigIons = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe','MgIIa','MgIIb','MgIa']
			epochs = {'visit1':4.98, 'visit2': 44.83, 'visit3': 289.58}
			order5 = []
			bins = {1:17, 2:17, 3:17}
			linsp = True
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		elif (sn == '06x' or sn == '06X'):
			exclude = []
			ions = ['NaDa', 'NaDb']
			sigIons = []
			epochs = {'visit1': -2, 'visit2': 14, 'visit3': 61, 'visit4': 105}
			order5 = []
			bins = {1:17, 2:17, 3:17, 4:17}
			linsp = True
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		elif(sn == '11by'):
			"""
			#11by params (all lines)
			exclude = [2343.46, 2359.02, 2381.05, 2398.43,2379.71, 2397.94, 2791.84, 2822.68, 2847.22, 2870.18]
			include = [2344.07, 2355.18, 2381.25, 2386.91, 2389.33, 2395.6, 2381.44, 2386.43, 2389.35, 2395.58, 2793.1, 2800.65, 2800.66, 2809.46, 2809.47, 2816.07, 2850.05, 2858.23, 2858.24, 2865.78]
			ions = ['FeIIa', 'FeIIb', 'FeIIc', 'FeIId', 'MgIIa', 'MgIIb', 'MgIIc', 'MgIId', 'MgIa', 'MgIb']
			epochs = {'visit1': 1, 'visit2': 2}
			setGlobals(exclude, include, ions, epochs)
			"""

			#11by params 
			#exclude = [2340.51, 2359.91, 2375.76, 2395, 2579.76, 2605.09, 2788.53, 2822.68, 2849.89, 2869.55]
			exclude = [2792.22, 2818.13, 2849.89, 2869.55]
			ions = ['MgIIa', 'MgIIb', 'MgIIc', 'MgIId', 'MgIa', 'MgIb']
			sigIons = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe']
			epochs = {'visit1': -9.31, 'visit2': -0.52}
			order5 = []
			bins = {1:17, 2:17}
			linsp = True
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		elif(sn == '11ek'):
			#11ek params
			exclude = [2802.45, 2823.52]
			ions = []
			sigIons = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe','MgIIa','MgIIb','MgIa']
			epochs = {'visit1': -3.32, 'visit2': 3.53}
			order5 = []
			bins = {1:17, 2:17}
			linsp = False
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		elif(sn == '11fe'):
			#11fe params
			exclude = [2340.99, 2350.32, 2371.01, 2389.01, 2583.5, 2594.17, 2599.95, 2605.43, 2790.31, 2811.01, 2851.12, 2860.45]
			#exclude = [2341.34, 2350.47, 2372.48, 2390.04, 2583.58, 2607.58, 2789.51, 2813.78, 2848.5, 2862.16]
			ions = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe','MgIIa','MgIIb','MgIa']
			sigIons = []
			epochs = {'visit1': -13.110599,'visit2': -10.041886, 'visit3': -6.9081770, 'visit4': -2.9512300, 'visit5': 0.039233502, 'visit6': 3.2484383, 'visit7': 9.1496879, 'visit8': 20.687449, 'visit9': 26.685855, 'visit10': 40.447814}
			order5 = []
			bins = {1:17, 2:17, 3:17, 4:17, 5:17, 6:17, 7:17, 8:17, 9:17, 10:17} 
			linsp = True
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		elif(sn == '11iv'):
			#11iv params	
			exclude = [2791.63, 2812.83]
			ions = ['MgIIa', 'MgIIb']
			sigIons = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe','MgIIa','MgIIb','MgIa']
			epochs = {'visit1': .6, 'visit2': 5.41, 'visit3': 10.15, 'visit4': 13.99, 'visit5': 17.82, 'visit6': 21.72, 'visit7': 29.52}
			order5 = []
			bins = {1:17, 2:17, 3:17, 4:17, 5:17, 6:17, 7:17} 
			linsp = False
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		elif(sn == '12cg'):
			#12cg params
			exclude = [2790.48, 2814.70, 2850.22, 2865.17]
			ions = ['MgIIa', 'MgIIb', 'MgIIc', 'MgIId', 'MgIa',]
			sigIons = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe']
			epochs = {'visit1': 2.5, 'visit2':16.37}
			order5 = []
			bins = {1:17, 2:17}
			linsp = True
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		elif(sn == '13dy'):
			#13dy params (w/out SiII)
			exclude = [2339, 2360, 2375, 2400, 2585, 2615, 2790.77, 2824.9, 2846.13, 2873.7]
			ions = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe', 'MgIIa', 'MgIIb', 'MgIIc', 'MgIId', 'MgIa', 'MgIb']
			ionsLim = ['MgIIa', 'MgIIb', 'MgIIc', 'MgIId']
			sigIons = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe']
			epochs = {'visit1': -6.17, 'visit2': -2.09, 'visit3':  -0.37, 'visit4':  1.61, 'visit5':  4.85, 'visit6':  8.76, 'visit7':  12.36, 'visit8':  14.38, 'visit9':  18.22, 'visit10':  21.17}
			order5 = []
			bins = {1:10, 2:10, 3:17, 4:10, 5:17, 6:10, 7:10, 8:12, 9:17, 10:10} 
			linsp = False
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		elif(sn == '14j' or sn == '14J'):
			#14j params
			exclude = [2791.63, 2812.83, 2846.1, 2864.73]
			ions = ['MgIIa', 'MgIIb', 'MgIa']
			sigIons = ['FeIIa','FeIIb','FeIIc','FeIId','FeIIe']
			epochs = {'visit1': -6.4956025 , 'visit2': -4.4969556, 'visit3':  -2.4983087, 'visit4':  -0.49966173 , 'visit5':  2.4983087, 'visit6':  6.4956025, 'visit7':  8.4942495, 'visit8':  11.492220 , 'visit9':  14.490190, 'visit10':  24.052866}	
			order5 = []
			bins = {1:17, 2:17, 3:17, 4:17, 5:17, 6:17, 7:17, 8:17, 9:17, 10:17} 
			linsp = False
			setGlobals(exclude, ions, epochs, order5, bins, sigIons, linsp)
			cont = False
		else:
			continue
			
exclude = []
ions = []
epochs = {}

rests = [2344.2128, 2374.4601, 2382.7641, 2586.6494, 2600.1722, 2796.354, 2803.5311, 2852.9628]
