##############################################################################
#
# author: Aaron Beaudoin
#
# Spectral line analysis tool to be run in ipython.
#
# Directory structure:
# spAnalyze.py/"Supernova name"/Data/:
#	breakpoints
#	data
#	eqw/
#		errors
#		MonteCarlo/Histograms
#	include
#	spectra
#
# Example code in Run.py in parent directory
#
# Uses astropy, numpy, scipy, matplotlib, glob, sys, and math packages.
#############################################################################

from astropy.table import Table
import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from scipy import optimize
from scipy import integrate	
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import glob
import sys
import math

#initializes spectrum class
class spectrum(object):
	"""Attributes can be added"""

#sets parameters used throughout
#param array ex- x values of lines to exclude in normalizing spectrum.
#param array g- guesses for height, wavelength, width for each absorption line.
#param array inc- x values of lines to include in Gaussian fitting.
#param array i- List of ions to search for.
#param array ep- Epochs of each visit. Ex. {'visit1': 0.111, ...}
def setGlobals(ex, g, inc, i, ep):
	global exclude
	global guesses
	global include
	global ions
	global epochs
	exclude = ex
	guesses = g
	include = inc
	ions = i
	epochs = ep

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
def interpolate(sp, bp, rand = False):
	x = sp.wave
	y = sp.flux
	#remove absorption lines from spectrum
	lines = np.append([], np.where(x<exclude[0]))
	lines = [int(i) for i in lines]
	i = 1
	while (i<len(exclude)-1):
		lines = np.append(lines, np.where((x>exclude[i]) & (x<exclude[i+1])))
		i += 2	
	lines = np.append(lines, np.where(x>exclude[len(exclude)-1]))
	
	#create b-spline over spectrum
	sm1 = inter.splrep(x[lines],y[lines])
	xx = linspace(min(x), max(x), len(x) / 10)
	
	if (rand):
		xx += np.random.randn(len(xx))
		indecies = []
		i = len(xx) - 10
		while (i < len(xx)):
			if (xx[i] < min(x)):
				indecies.append(i)
			if (xx[i] > max(x)):
				indecies.append(i)
			i += 1
		xx = np.delete(xx, indecies)
	
	linesxx = np.append([], np.where(xx<exclude[0]))
	linesxx = [int(i) for i in linesxx]
	i = 1
	while(i<len(exclude)-1):
		linesxx = np.append(linesxx, np.where((xx>exclude[i]) & (xx<exclude[i+1])))
		i += 2	
	linesxx = np.append(linesxx, np.where(xx>exclude[len(exclude)-1]))
	xx = xx[linesxx]
	
	#add breakpoints to b-spline
	xx = np.append(xx,bp)
	xx.sort()
	
	#fits b-spline over wavelength range
	y1 = inter.splev(xx,sm1)
	sm2 = inter.splrep(xx, y1)
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
	
	dir = file[:14]
	
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
		
#helper function to findErrors()		
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1:
			raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.size < window_len:
			raise ValueError, "Input vector needs to be bigger than window size."
	if window_len<3:
			return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
			raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
	if window == 'flat': #moving average
			w=np.ones(window_len,'d')
	else:  
			w=eval('np.'+window+'(window_len)')
	y=np.convolve(w/w.sum(),s,mode='same')
	return y[window_len:-window_len+1]

#finds noise error in spectrum	
def findErrors(y,sflux):	
	sf=smooth(y)
	res=y-sf
	ares=np.abs(res)
	se = smooth(ares)
	norm_e = se/sflux
	return norm_e

#Gets spectrum from FLM file	
def getSpectrum(file, sp):
	data = Table.read(file,format='ascii')
	sp.wave = data["col1"]
	sp.flux = data["col2"]

#Finds and records new breakpoints and records normalized spectrum
def findNewBkpts(sp):
	x = sp.wave
	y = sp.flux
	k = 1
	sflux = 0
	while(k == 1):
		points = input('Breakpoints for ' + sp.name + ': ')
		bps = [float(i) for i in points]
		sp_flux, norm = interpolate(sp, bps)
		plt.figure(1)
		plt.plot(x, y)
		plt.plot(x, sp_flux)
		plt.show()
		plt.figure(2)
		plt.plot(x, norm)
		plt.ylim(0,1.5)
		plt.show()
		cont = raw_input('More breakpoints?: (y/n) ')
		if (cont == 'y'):
			continue
		else:
			sp.norm_flux = norm
			sp.sflux = sp_flux
			sp.bkpts = bps
			k = 2
	
	temp = np.zeros(len(sp.bkpts))
	t = Table([sp.bkpts, temp], names = ('col1','col2'))
	t.write('sn' + sp.sn + '/Data/breakpoints/' + sp.name + '.flm', format='ascii')
	sp.cont_err = findErrors(y, sp.sflux)

#retrieves previously calculated breakpoints from FLM files and records normalized spectrum
def getBkpts(sp):
	t = Table.read('sn' + sp.sn + '/Data/breakpoints/' + sp.name + '.flm', format='ascii')
	sp.bkpts = t["col1"]
	sp_flux, norm = interpolate(sp, sp.bkpts)
	sp.norm_flux = norm
	sp.sflux = sp_flux
	sp.cont_err = findErrors(sp.flux, sp.sflux)

#helper function to monteCarlo(), randomizes breakpoints for variation in EQWs	
def randBreak(breakpoints):
	bkpts = []
	for b in breakpoints:
		bkpts.append(b)
	bshift = np.random.randn(len(bkpts)) * 2.5
	new_bkpts = bshift + bkpts
	return new_bkpts

#find EQWs for lines in spectrum sp	
def findEqw(sp, plot):
	include = sp.inc
	w = sp.wave
	f = sp.norm_flux
	e = sp.cont_err
	
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
	g=0
	eqws = []
	eqw_errs = []
	param = []
	while(g<len(guesses)):
		inc = np.where((w>include[i]) & (w<include[i+1]))
		wave = w[inc]
		flux = f[inc]
		errs = e[inc]
		
		x = np.array(wave)
		y = np.array(flux)
		err = np.array(errs)
		
		rest1 = 2796.354
		rest2 = 2803.5311
		
		wobs = 0
		fmin = 24352345
		for j in range(len(x)):
			if (y[j] < fmin):
				fmin = y[j]
				wobs = x[j]
						
		p = [fmin - 1, wobs, guesses[g+2]]
		p1, pcov, infodict, errmsg, success= optimize.leastsq(residual, p[:], args=(x,y), full_output = 1)
		gaus = gaussian(p1, x)
		for j in p1:
			param.append(j)
		
		if(len(w) > len(p))	 and pcov is not None:
			s_sq = (residual(p1, x, y)**2).sum()/(len(w)-len(p))
			pcov = pcov * s_sq
		else:
			pcov = np.inf
		
		errors = []
		for j in range(len(p1)):
			try:
				errors.append(np.abs(pcov[j][j])**.5)
			except:
				errors.append(0.00)
			
		errors = np.array(errors)
		
		error = np.sqrt(np.sqrt(math.pi) * ((errors[0]*p1[2])**2 + (errors[2]*p1[0])**2))
		
		if(plot):
			for j in range(len(gaus)):
				baseline_x.append(x[j])
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
		g+=3

	if(plot):
		for e in range(len(end[0])):
			baseline_x.append(w[end[0][e]])
			baseline_y.append(1)
		
		plt.figure()
		plt.plot(baseline_x,baseline_y,label='gaus')
		plt.plot(w,f,label='cont')
		legend = plt.legend(loc='upper right')
		plt.show()
	
	sp.eqws = eqws
	sp.eqw_errs = eqw_errs
	sp.param = param

#helper function to monteCarlo() which calculates standard deviations	
def findErr(temp, original, name, sn, ion):
	eqws = temp
	eqws.sort()
	n = len(eqws)
	
	err_lo = eqws[int(round(.16*n))]
	err_hi = eqws[int(round(.84*n))]
	med = np.median(eqws)
	
	mc_lo = med - err_lo
	mc_hi = err_hi - med
	
	plt.figure(0)
	hist,bins=np.histogram(eqws,bins=len(eqws)/2)
	width=.7*(bins[1]-bins[0])
	center=(bins[:-1]+bins[1:])/2
	plt.bar(center,hist,align='center',width=width)
	plt.plot([err_lo],[2], 'go', label='err_lo')
	plt.plot([err_hi],[2],'co',label='err_hi')
	plt.plot([med],[2],'ko', label='median')
	plt.plot([original],[2],'ro', label='original')
	legend=plt.legend(loc='upper right')
	plt.savefig('sn' + sn + '/Data/eqw/MonteCarlo/Histograms/MCHist' + name + '-' + ion + '.png',orientation='landscape')
	plt.close()
	
	return mc_lo, mc_hi

#def find_red(sp):	
	gaussian = lambda p, w, x: 1+p[0]*np.exp(-(x-w)**2/(2*p[1]**2))
	redshift = lambda obs, emit: (obs / emit) - 1
	wave_obs = lambda z, emit: (z + 1) * emit
	
	w = sp.wave
	f = sp.norm_flux
	
	inc = np.where((w > include[0]) & (w<include[1]))
	wave = w[inc]
	flux = f[inc]
	
	x = np.array(wave)
	y = np.array(flux)
	
	wmin = 0
	fmin = 12341234
	for i in range(len(x)):
		if(y[i] < fmin):
			fmin = y[i]
			wmin = x[i]
	
#Separates MgII doublets at 2800 angstroms from galaxy absorption and Milky Way absorption.
def sepMgDoub(sp):
	#inc = [2790.72, 2799.98, 2799.99, 2811.93, 2811.94, 2822.54]
	#inc = [2792.49, 2802.9, 2802.91, 2813.53, 2813.54, 2823.94]
	inc = sp.inc
	rest1 = 2796.354
	rest2 = 2803.5311

	guess = [2.02, 3.27]
	#guess = [sp.param[2], sp.param[5]]

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

	wave1 = w[inc1]
	flux1 = f[inc1]
	wave2 = w[inc2]
	flux2 = f[inc2]
	wave3 = w[inc3]
	flux3 = f[inc3]

	x1 = np.array(wave1)
	y1 = np.array(flux1)
	x2 = np.array(wave2)
	y2 = np.array(flux2)
	x3 = np.array(wave3)
	y3 = np.array(flux3)

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
	
	

	p = [fmin1 - 1, guess[0], z1, fmin2 - 1, guess[1], z2]

	p1, pcov, infodict, errmsg, success = optimize.leastsq(residual, p[:], args = (rest1, rest2, x3, y3), full_output = 1)

	sp.comb = combGaus(p1, rest1, rest2, w) + 1
	sp.gausB = gaussian([p1[0], p1[1], p1[2]], rest2, w) + 1
	sp.gausC = gaussian([p1[3], p1[4], p1[5]], rest1, w) + 1 

	if (len(w) > len(p)) and pcov is not None:
		s_sq = (residual(p1, rest1, rest1, x3, y3)**2).sum()/(len(w) - len(p))
		pcov = pcov * s_sq
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
	
	error.append(np.sqrt(np.sqrt(math.pi) * ((errors[0]*p1[1])**2 + (errors[1]*p1[0])**2)))	
	error.append(np.sqrt(np.sqrt(math.pi) * ((errors[4]*p1[3])**2 + (errors[3]*p1[4])**2)))
	
	eqws = []	
	eqws.append(np.abs(integrate.quad(lambda x: p1[0] * np.exp(-(x-rest2*(1+p1[2]))**2/(2*p1[1]**2)), min(x1), max(x2))[0]))
	eqws.append(np.abs(integrate.quad(lambda x: p1[3] * np.exp(-(x-rest1*(1+p1[5]))**2/(2*p1[4]**2)), min(x1), max(x2))[0]))
	
	count = len(sp.eqws)
	oldEqws = sp.eqws
	newEqws = eqws
	
	oldErr = sp.eqw_errs
	newErr = error
	
	sp.eqws = [float(oldEqws[0]), float(newEqws[0]), float(newEqws[1]), float(oldEqws[2])]
	sp.eqw_errs = [float(oldErr[0]), float(newErr[0]), float(newErr[1]), float(oldErr[2])]
	
	i = 3
	while (i < count):
		sp.eqws.append(float(oldEqws[i]))
		sp.eqw_errs.append(float(oldErr[i]))
	
#finds systematic errors for EQWs
def monteCarlo(mc, sp, triplet = False):
	print "Finding error for", sp.name
	spectra=[]
	for i in range(mc):
		nsp = spectrum()
		nsp.wave = sp.wave
		nsp.flux = sp.flux
		nsp.bkpts = randBreak(sp.bkpts)
		nsp.inc = sp.inc
		sp_flux, norm = interpolate(nsp, nsp.bkpts, True)
		nsp.norm_flux = norm
		nsp.sflux = sp_flux
		nsp.cont_err = findErrors(nsp.flux, nsp.sflux)
		findEqw(nsp,False)
		if (triplet):
			sepMgDoub(nsp)
		spectra.append(nsp)
	
	mc_lo_errs = []
	mc_hi_errs = []
	for i in range(len(spectra[0].eqws)):
		temp = []
		for j in range(len(spectra)):
			temp.append(spectra[j].eqws[i])
		mc_lo, mc_hi = findErr(temp, sp.eqws[i], sp.name, sp.sn, ions[i])
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

#saves all data to FLM files
def saveData(sp):
	att = dir(sp)
	if('inc' in att):
		saveInc(sp)
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

#loads any saved data from FLM files. Requires a directory
#Returns array of spectra with loaded data
def loadData(dir):
	spectra = []
	
	#load spectrum info
	files = glob.glob(dir + '/data/*.flm')
	for f in files:
		sp = spectrum()
		spNames(f, sp)
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
			
		getBkpts(sp)
		loadInc(sp)
		spectra.append(sp)
	
	return spectra

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
		
	font = {'family' : 'serif','color'  : 'black','weight' : 'bold','size' : 10,} 
		
	fig = plt.figure(num = 5, dpi = 100, figsize = [10, (3*len(spec_num)) + 2], facecolor = 'w')
	gs = gridspec.GridSpec(1,1)
	plt.title(spectra[0].sn, fontdict = font)
	plt.gca().axes.get_xaxis().set_visible(False)
	plt.gca().axes.get_yaxis().set_visible(False)
	count = 0
	for j in range(len(species)):
		k = 0
		ax = fig.add_subplot(len(species), 1, j+1)
		plt.gca().set_color_cycle(['black', 'blue', 'cyan', 'green'])
		while (k < spec_num[j]):
			ax.errorbar(epochs, eqws[count], yerr = [lo_err[count], hi_err[count]], label = ions[count])
			k += 1
			count += 1
		ax.text(.5, .9, species[j], horizontalalignment = 'center', transform = ax.transAxes, fontdict = font)
		plt.xlim(epochs[0], epochs[len(epochs) - 1])
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
		plt.gca().set_color_cycle(['black', 'blue', 'cyan', 'green'])
		while (k < spec_num[j]):
			ax.errorbar(epochs, eqws[count] - np.average(eqws[count]), yerr = [lo_err[count], hi_err[count]], label = ions[count])
			k += 1
			count += 1
		ax.text(.5, .9, species[j], horizontalalignment = 'center', transform = ax.transAxes, fontdict = font)
		plt.xlim(epochs[0], epochs[len(epochs) - 1])
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

"""
exclude = [2789.67, 2826.71, 2843.33, 2878.64, 6025.46, 6344.41]
include = [2785.68, 2800.03, 2800.03, 2812.02, 2812.03, 2822.5, 6040.96, 6316.41]
guesses = [-.3, 2797.1, 3.4, -.45, 2807.77, 7.4, -.34, 2815.21, 4.2, -.6199, 6170.77, 119.02]
ions = ['MgIIa', 'MgIIb', 'MgIIc', 'SiIIa']
epochs = [.5023, 4.61405, 6.33734, 8.32473, 11.5869, 15.5064, 19.1053, 21.1529, 25.0071, 27.9399]

"""
#13dy params (w/out SiII)
exclude = [2789.67, 2826.71, 2843.33, 2878.64]
include = [2790.72, 2799.98, 2799.99, 2811.93, 2811.94, 2822.54]
guesses = [-.3, 2797.1, 3.4, -.45, 2807.77, 7.4, -.34, 2815.21, 4.2]
ions = ['MgIIa', 'MgIIb', 'MgIIc', 'MgIId']
epochs = {'visit1': .5023, 'visit2': 4.61405, 'visit3':  6.33734, 'visit4':  8.32473, 'visit5':  11.5869, 'visit6':  15.5064, 'visit7':  19.1053, 'visit8':  21.1529, 'visit9':  25.0071, 'visit10':  27.9399}
"""
files = glob.glob('sn2013dy/Data/data/*.flm')
spectra = []
for f in files:
	sp = spectrum()
	spNames(f, sp)
	getSpectrum(f, sp)
	getBkpts(sp)
	findEqw(sp, False)
	spectra.append(sp)

"""
#11fe params
"""
exclude = [2341.34, 2350.47, 2372.48, 2390.04, 2583.58, 2607.58, 2789.51, 2813.78, 2848.5, 2862.16]
include = [2340.35, 2353.69, 2371.61, 2379.75, 2380, 2390.5, 2581.61, 2593.2, 2597.06, 2607.95, 2786.39, 2801.95, 2802, 2815.82, 2848.86, 2860.93]
guesses = [-.3675, 2345, 3.7, -.205, 2375.98, 3.62, -.466, 2385.51, 3.99, -.3563, 2588.31, 4.5, -.464, 2602, 4.4, -.6416, 2798.05, 4.37, -.596, 2806.09, 4.1, -.1342, 2854.61, 4.29]
ions=['FeIIa','FeIIb','FeIIc','MnIIa','MnIIb','MgIIa','MgIIb','MgIa']
epochs = {'visit1': -13.110599,'visit2': -10.041886, 'visit3': -6.9081770, 'visit4': -2.9512300, 'visit5': 0.039233502, 'visit6': 3.2484383, 'visit7-obs1': 9.1496879, 'visit7-obs2+3': 9.1496879, 'visit8': 20.687449, 'visit9': 26.685855, 'visit10': 40.447814}
"""
#dir = 'sn2011fe/Data'
#spectra = loadData(dir)
"""

#!!!example code!!!
#Reads in existing calculations and plots them
#dir = 'sn2011fe/Data'
#spectra = loadData(dir)
#plotEqw(spectra)



#Finds the equivalent widths and total error for the lines in visit1
file = 'sn2011fe/Data/data/sn2011fe-visit1-uv.flm'
sp = spectrum()
spNames(file, sp)
getSpectrum(file, sp)
getBkpts(sp)
findEqw(sp, False)
monteCarlo(5, sp)
calcErr(sp)

#Finds the equivalent widths and total error for the lines in every visit for 11fe
files = glob.glob('sn2011fe/Data/data/*.flm')
spectra = []
for f in files:
	sp = spectrum()
	spNames(f, sp)
	getSpectrum(f, sp)
	getBkpts(sp)
	findEqw(sp, False)
	monteCarlo(500, sp)
	calcErr(sp)
	saveData(sp)
	spectra.append(sp)

plotEqw(spectra)


# 13dy params
exclude13 = [2785,2826,2843,2874,2964,3036,3113,3354,3580,4044,4188,4516,4800,5144,5953,6339]
exclude = [2789.67, 2826.71, 2843.33, 2878.64]
include = [2790.72, 2799.98, 2799.99, 2811.93, 2811.94, 2822.54]
2811.94-> 2807, 2800 -> 2807.5
guesses = [-.3, 2797.1, 3.4, -.45, 2807.77, 7.4, -.34, 2815.21, 4.2]
ions = ['MgIIa', 'MgIIb', 'MgIIc', 'MgIId', 'SiIIa']
epochs = [.5023, 4.61405, 6.33734, 8.32473, 11.5869, 15.5064, 19.1053, 21.1529, 25.0071, 27.9399]
"""