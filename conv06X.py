import spAnalyze as a
import numpy as np
import glob
import math
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve as con
from astropy.table import Table
from scipy.interpolate import interp1d
from scipy.stats import chi2

center = {1.0: 5928.08, 2.0: 5928.21, 3.0: 5928.63, 4.0: 5928.92}
waver = {1.0: [2803.25, 2803.73], 2.0: [2802.65, 2803.76], 3.0: [2802.03, 2803.85], 4.0: [2801.55, 2803.74]}

files = glob.glob('sn2006X/Data/data/*.dat')
spectra = []
a.pickSN('06x')

for f in files:
	spectra.append(a.initSpec(f))
	
 
#convolves spectrum (sp) to the specified resolution (res) and pixel scale (rebin)
def convolve(sp, res = 1.3, rebin = 1.5, params = True):
	pix  = np.round(sp.wave[1] - sp.wave[0], 4)    #pixel scale of sp
	num = (3000 - 2600) / pix     #number of pixels between 2600-3000 angstroms
	num = int(np.round(num, 0))

	#creates synthesised spectrum with coordinates (swave (wavelength), sflux (flux))
	swave = np.arange(2600, 2600 + num * pix, pix)
	sflux = []
	for s in swave:
		sflux.append(1)
	
	#Scales the wavelength values to simulate the Na doublet as the Mg doublet at 2796 $ 2802 angstroms
	sp.norm_flux = sp.flux
	plt.figure()
	plt.plot(sp.wave, sp.norm_flux)
	cont = True
	while(cont):
		if (params):
			try:
				temp = raw_input('Select wavelength of center of Na line: ')
				temp = float(temp)
				cont = False
			except NameError:
				print 'Invalid'
		else:
			temp = center[sp.visit]
			cont = False
	plt.close()
	sp.wave = sp.wave / (temp / 2803.5311)

	#finds the wavelength range for the new Na lines
	plt.figure()
	plt.plot(sp.wave, sp.norm_flux)
	plt.ylim(0, 1.5)
	cont = True
	while(cont):
		if params:
			try:
				cent = input('Select wavelength range of Na line: ')
				cont = False
			except NameError:
				print 'Invalid'
		else:
			cent = waver[sp.visit]
			cont = False
	plt.close()
	inc = np.where((sp.wave > cent[0]) & (sp.wave < cent[1]))[0]
	na = sp.flux[inc]      # flux values of sodium lines
	
	num = len(na)
	start = inc[num / 2] - inc[0]      #number of pixels between middle and start of absorption feature
	
	#finds value and index of pixel closest to the absorption peak
	mgiia = 2796.354
	mgiib = 2803.5311
	mida = 0
	ia = 0
	midb = 0
	ib = 0
	diffa = 9999
	diffb = 9999
	for i in range(len(swave)):
		w = swave[i]
		if np.abs(mgiia - w) < diffa:
			mida = w
			diffa = np.abs(mgiia - w)
			ia = i
		if np.abs(mgiib - w) < diffb:
			midb = w
			diffb = np.abs(mgiib - w)
			ib = i
			

	#adds absorption features to synthesized spectrum
	sa = ia - start
	sb = ib - start
	for n in na:
		sflux[sa] = n
		sflux[sb] = n
		sa += 1
		sb += 1
		
	#finds the number of pixels in the new resolution from the original spectral resolution
	pix_res = 0
	pix_rebin = 0
	for w in range(len(swave)):
		diff = swave[w] - swave[0]
		if diff < res:
			pix_res += 1
		if diff < rebin:
			pix_rebin += 1
			
		
	#convolves the spectrum to the given resolution (res)
	ker = Gaussian1DKernel(pix_res)
	conv = con(sflux, ker, boundary = 'extend')
		
	#rebins the convolved data to the given pixel scale (rebin)
	rf = []
	rw = []
	i = 0
	while (i < len(swave)):
		tempw = []
		tempf = []
		it = 0
		while(it < pix_rebin):
			tempw.append(swave[i])
			tempf.append(conv[i])
			it += 1
		rw.append(np.average(tempw))
		rf.append(np.average(tempf))
		i += pix_rebin
		
	#creates a new spectrum object with the synthesized spectrum
	s = a.spectrum()
	s.name = sp.name + '_Conv' 
	s.visit = sp.visit
	s.sn = sp.sn
	s.dir = sp.dir
	s.epoch = sp.epoch
	s.wave = np.array(rw)
	s.norm_flux = np.array(rf)
	s.flux = s.norm_flux
	return s
	

#finds linear parameters for epoch vs equivalent width. spectra is an array of 06X spectra
#convolved to a new resolution and pixel scale. loadData determines if old EW values should
#be used or if new ones should be measured.
def linFit(spec, loadData = False):
	#gets equivalent width values of convolved data
	if ((glob.glob(spec[0].dir + 'eqw/*.flm') == []) & loadData):
		loadData = False
	if loadData:
		files = glob.glob(spec[0].dir + 'eqw/*.flm')
		i = 0
		for f in files:
			data = Table.read(f, format = 'ascii')
			spec[i].eqws = np.array(data['eqws'])
			spec[i].eqw_errs = np.array(data['eqw_errs'])
			a.loadInc(spec[i])
			i += 1
	else:
		for s in spec:
			plt.figure()
			plt.plot(s.wave, s.norm_flux)
			a.setInc(s)
			plt.close()
		for s in spec:
			a.findEqw(s)
			a.saveData(s)

	#finds linear fit parameters between the first-second and second-third epochs
	m1 = (spec[1].eqws[0] - spec[0].eqws[0]) / (spec[1].epoch - spec[0].epoch)
	b1 = spec[0].eqws[0] - m1 * spec[0].epoch
	m2 = (spec[2].eqws[0] - spec[1].eqws[0]) / (spec[2].epoch - spec[1].epoch)
	b2 = spec[1].eqws[0] - m2 * spec[1].epoch
	return [[m1, b1], [m2, b2]]

def testMC(spec, spectra, mc = 20):
	epochs = [spec[0].epoch - 18]	#epochs used for linear fit, starts 18 days before first epoch of 06X
	eqws = [list(spec[0].eqws)]	#EWs used for linear fit, assumes no change before first spectrum
	for s in spec:
		epochs.append(s.epoch)
		eqws.append(list(s.eqws))
	epochs.append(s.epoch + 18)
	eqws.append(list(s.eqws))
	eqws = np.array(eqws)

	#creates synthesized spectra based on epoch
	#assumes a linear change in EW strength between epochs
	#assumes no change before first epoch and after last epoch
	wave = []
	flux = []
	incs = []
	
	for sp in spectra:
		cont = True
		if (sp.epoch < spec[0].epoch):
			wave.append(spec[0].wave)
			flux.append(spec[0].norm_flux)
			incs.append(spec[0].inc)
			cont = False
		i = 0
		while ((i < len(spec) - 1) & (cont)):
			if((sp.epoch > spec[i].epoch) & (sp.epoch < spec[i + 1].epoch)):
				diff = spec[i + 1].epoch - spec[i].epoch
				diff1 = sp.epoch - spec[i].epoch
				diff2 = spec[i + 1].epoch - sp.epoch
				p1 = 1 - diff1/ diff
				p2 = 1 - diff2 / diff
				w = spec[i].wave
				f = spec[i].norm_flux * p1 + spec[i + 1].norm_flux * p2
				wave.append(w)
				flux.append(f)
				inc = spec[i + 1].inc
				incs.append(inc)
				cont = False
			i += 1
		if (cont):
			wave.append(spec[len(spec) - 1].wave)
			flux.append(spec[len(spec) - 1].norm_flux)
			incs.append(spec[len(spec) - 1].inc)
	
			
	#finds the RMS from 'spectra'. Assumes no absorption between 2650 and 2750 angstroms.
	#RMS = sqrt(sum(residual^2) / len(residual)) where residual = 1 - normalized flux
	rms = []
	for sp in spectra:
		inc = np.where((sp.wave > 2870) & (sp.wave < 2950))[0]
		res = 1 - sp.norm_flux[inc]
		temp = []
		for r in res:
			temp.append(r**2)
		rms.append(np.sqrt(np.sum(temp) / (len(res))))

	#adds random noise based on RMS values. Random values chosen from normal distribution.
	print 'Adding random noise to synthesized spectra'
	nspec = []
	for m in range(mc):
		i = 0
		ns = []
		for r in rms:
			n = a.spectrum()
			n.wave = wave[i]
			n.inc = incs[i]
			n.norm_flux = []
			n.sn = spec[0].sn
			n.epoch = spectra[i].epoch
			for f in range(len(flux[i])):
				n.norm_flux.append(flux[i][f] + np.random.normal(0, .68, 1)[0] * r)
			n.norm_flux = np.array(n.norm_flux)
			a.findEqw(n)
			ns.append(n)
			i += 1
		nspec.append(ns)
		a.statusBar(m + 1, mc)
	
	for ns in nspec:
		ns.sort(key = lambda x: x.epoch)
	
	#plots scatter of EW vs. epochs
	eps = []
	ews = []
	errs = []
	
	for ns in nspec:
		ep = []
		ew = []
		err = []
		for n in ns:
			ep.append(n.epoch)
			ew.append(n.eqws)
			err.append(n.eqw_errs)
		eps.append(ep)
		ews.append(ew)
		errs.append(err)
		
	chisq = []
	for e in range(len(epochs) - 1):
		chis = []
		
		temp = []
		ddof = 0
		for n in ns:
			if ((n.epoch > epochs[e]) & (n.epoch < epochs[e + 1])):
				if n.eqw_errs[0] != 0:
					temp.append(n)
				else:
					ddof += 1
		if (e == 0 | e == len(epochs) - 1):
			ddof += 0
		else:
			ddof += 1
		
		exp = []
		if eqws[e][0] == eqws[e + 1][0]:
			for t in temp:
				exp.append(eqws[e][0])
		else:
			for t in temp:
				diff = epochs[e + 1] - epochs[e]
				diff1 = t.epoch - epochs[e]
				diff2 = epochs[e + 1] - t.epoch
				p1 = 1 - diff1 / diff
				p2 = 1 - diff2 / diff
				exp.append(p1 * eqws[e][0] + p2 * eqws[e + 1][0])
		sum = 0
		for i in range(len(temp)):
			obs = temp[i].eqws[0]
			sig = temp[i].eqw_errs[0]
			sum += ((obs - exp[i]) / sig)**2
		chis.append([sum, len(temp) - 1 - ddof])
		chisq.append(chis)
		
		temp = []
		ddof = 0
		for n in ns:
			if ((n.epoch > epochs[e]) & (n.epoch < epochs[e + 1])):
				if n.eqw_errs[0] != 0:
					temp.append(n)
				else:
					ddof += 1
		if (e == 0 | e == len(epochs) - 1):
			ddof += 0
		else:
			ddof += 1
		
		exp = []
		if eqws[e][1] == eqws[e + 1][1]:
			for t in temp:
				exp.append(eqws[e][1])
		else:
			for t in temp:
				diff = epochs[e + 1] - epochs[e]
				diff1 = t.epoch - epochs[e]
				diff2 = epochs[e + 1] - t.epoch
				p1 = 1 - diff1 / diff
				p2 = 1 - diff2 / diff
				exp.append(p1 * eqws[e][1] + p2 * eqws[e + 1][1])
		sum = 0
		chi = []
		for i in range(len(temp)):
			obs = temp[i].eqws[1]
			sig = temp[i].eqw_errs[1]
			sum += ((obs - exp[i]) / sig)**2
			chi.append([sum, len(temp) - 1 - ddof])
		chis.append(chi)
		chisq.append(chis)
	
	"""
	FIX
	"""
	
	alpha = 0.01
	count1 = 0
	count2 = 0
	for chi in chisq:
		for c in chi:
			good = True
			for set in c:
				if chi2.pdf(set[0], set[1] < alpha):
					good = False
			if good:
				count1 += 1
			
	
#spec: perfect convolved spectra from 06X (created earlier)
#creates synthesized spectra of 06X at given epochs (from spectra)
#finds rms of spectra between 2650 and 2750 angstroms
#adds in random noise to synthesized spectrum mc times, determines how often we see a
#3-sigma change
def monteCarlo(spec, spectra, mc = 20):
	epochs = [spec[0].epoch - 18]	#epochs used for linear fit, starts 18 days before first epoch of 06X
	eqws = [list(spec[0].eqws)]	#EWs used for linear fit, assumes no change before first spectrum
	for s in spec:
		epochs.append(s.epoch)
		eqws.append(list(s.eqws))
	epochs.append(s.epoch + 18)
	eqws.append(list(s.eqws))
	eqws = np.array(eqws)

	#creates synthesized spectra based on epoch
	#assumes a linear change in EW strength between epochs
	#assumes no change before first epoch and after last epoch
	wave = []
	flux = []
	incs = []
	
	for sp in spectra:
		cont = True
		if (sp.epoch < spec[0].epoch):
			wave.append(spec[0].wave)
			flux.append(spec[0].norm_flux)
			incs.append(spec[0].inc)
			cont = False
		i = 0
		while ((i < len(spec) - 1) & (cont)):
			if((sp.epoch > spec[i].epoch) & (sp.epoch < spec[i + 1].epoch)):
				diff = spec[i + 1].epoch - spec[i].epoch
				diff1 = sp.epoch - spec[i].epoch
				diff2 = spec[i + 1].epoch - sp.epoch
				p1 = 1 - diff1/ diff
				p2 = 1 - diff2 / diff
				w = spec[i].wave
				f = spec[i].norm_flux * p1 + spec[i + 1].norm_flux * p2
				wave.append(w)
				flux.append(f)
				inc = spec[i + 1].inc
				incs.append(inc)
				cont = False
			i += 1
		if (cont):
			wave.append(spec[len(spec) - 1].wave)
			flux.append(spec[len(spec) - 1].norm_flux)
			incs.append(spec[len(spec) - 1].inc)
	
			
	#finds the RMS from 'spectra'. Assumes no absorption between 2650 and 2750 angstroms.
	#RMS = sqrt(sum(residual^2) / len(residual)) where residual = 1 - normalized flux
	rms = []
	for sp in spectra:
		inc = np.where((sp.wave > 2870) & (sp.wave < 2950))[0]
		res = 1 - sp.norm_flux[inc]
		temp = []
		for r in res:
			temp.append(r**2)
		rms.append(np.sqrt(np.sum(temp) / (len(res))))

	#adds random noise based on RMS values. Random values chosen from normal distribution.
	print 'Adding random noise to synthesized spectra'
	nspec = []
	for m in range(mc):
		i = 0
		ns = []
		for r in rms:
			n = a.spectrum()
			n.wave = wave[i]
			n.inc = incs[i]
			n.norm_flux = []
			n.sn = spec[0].sn
			n.epoch = spectra[i].epoch
			for f in range(len(flux[i])):
				n.norm_flux.append(flux[i][f] + np.random.normal(0, .68, 1)[0] * r)
			n.norm_flux = np.array(n.norm_flux)
			a.findEqw(n)
			ns.append(n)
			i += 1
		nspec.append(ns)
		a.statusBar(m + 1, mc)
	
	for ns in nspec:
		ns.sort(key = lambda x: x.epoch)
	
	#plots scatter of EW vs. epochs
	eps = []
	ews = []
	errs = []
	
	for ns in nspec:
		ep = []
		ew = []
		err = []
		for n in ns:
			ep.append(n.epoch)
			ew.append(n.eqws)
			err.append(n.eqw_errs)
		eps.append(ep)
		ews.append(ew)
		errs.append(err)
		
	print '\nCalculating goodness of fit for linear and consistent models'
	chisq1 = []	
	chisq0 = []
	count = 1
	for ns in nspec:
		chi = []
		
		observed = []
		errors = []
		weights = []
		ddof = 0
		for n in ns:
			observed.append(n.eqws[0])
			errors.append(n.eqw_errs[0])
			try:
				w = 1 / n.eqw_errs[0]
			except:
				w = 0
				ddof += 1
			weights.append(w)
			
		fit = np.polyfit(eps[0], observed, 1, w = weights)
		expected = []
		for e in eps[0]:
			expected.append(fit[0] * e + fit[1])
		sum = 0
		for i in range(len(observed)):
			sum += ((observed[i] - expected[i]) / errors[i])**2
		chi.append((sum, len(observed) - 2 - ddof))
		
		observed = []
		errors = []
		weights = []
		ddof = 0
		for n in ns:
			observed.append(n.eqws[1])
			errors.append(n.eqw_errs[1])
			try:
				w = 1 / n.eqw_errs[1]
			except:
				w = 0
				ddof += 1
			weights.append(w)
			
		fit = np.polyfit(eps[0], observed, 1, w = weights)
		expected = []
		for e in eps[0]:
			expected.append(fit[0] * e + fit[1])
		sum = 0
		for i in range(len(observed)):
			sum += ((observed[i] - expected[i]) / errors[i])**2
		chi.append((sum, len(observed) - 2 - ddof))
		
		chisq1.append(chi)
		
		chi = []
		
		observed = []
		errors = []
		weights = []
		ddof = 0
		for n in ns:
			observed.append(n.eqws[0])
			errors.append(n.eqw_errs[0])
			try:
				w = 1 / n.eqw_errs[0]
			except:
				w = 0
				ddof += 1
			weights.append(w)
			
		avg = np.average(observed, weights = weights)
		expected = []
		for o in observed:
			expected.append(avg)
		sum = 0
		for i in range(len(observed)):
			sum += ((observed[i] - expected[i]) / errors[i])**2
		chi.append((sum, len(observed) - 2 - ddof))
		
		observed = []
		errors = []
		weights = []
		ddof = 0
		for n in ns:
			observed.append(n.eqws[1])
			errors.append(n.eqw_errs[1])
			try:
				w = 1 / n.eqw_errs[1]
			except:
				w = 0
				ddof += 1
			weights.append(w)
			
		avg = np.average(observed, weights = weights)
		expected = []
		for e in eps[0]:
			expected.append(avg)
		sum = 0
		for i in range(len(observed)):
			sum += ((observed[i] - expected[i]) / errors[i])**2
		chi.append((sum, len(observed) - 2 - ddof))
		
		chisq0.append(chi)	
		
		a.statusBar(count, len(nspec))
		count += 1		
	
	fig = plt.figure(dpi = 100, figsize = [12, 16], facecolor = 'w')
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)
	for i in range(len(nspec)):
		max1 = np.array(ews[i])[:,0].max()
		min1 = np.array(ews[i])[:,0].min()
		fit = np.polyfit(eps[i], np.array(ews[i])[:,0], 1, w = np.array(errs[i])[:,0])
		y = [fit[0] * e + fit[1] for e in eps[i]]
		ax1.plot(eps[i], y, color = 'c')
		ax1.plot(epochs, eqws[:,0], color = 'k', linestyle = 'solid', linewidth = 2.5)
		ax1.errorbar(eps[i], np.array(ews[i])[:,0], yerr = np.array(errs[i])[:,0], color = 'b', linestyle = 'none')
		ax1.scatter(eps[i], np.array(ews[i])[:,0], color = 'b')
		
		xmin, xmax = ax1.get_xlim()
		ax1.set_xlim(xmin - 2.5, xmax + 2.5)
		ax1.set_ylim(min1 - .2, max1 + .2)
		
		
		max2 = np.array(ews[i])[:,1].max()
		min1 = np.array(ews[i])[:,1].min()
		fit = np.polyfit(eps[i], np.array(ews[i])[:,1], 1, w = np.array(errs[i])[:,1])
		y = [fit[0] * e + fit[1] for e in eps[i]]
		ax2.plot(eps[i], y, color = 'c')
		ax2.plot(epochs, eqws[:,1], color = 'k', linestyle = 'solid', linewidth = 2.5)
		ax2.errorbar(eps[i], np.array(ews[i])[:,1], yerr = np.array(errs[i])[:,1], color = 'b', linestyle = 'none')
		ax2.scatter(eps[i], np.array(ews[i])[:,1], color = 'b')
		
		ax2.set_xlim(xmin - 2.5, xmax + 2.5)
		ax2.set_ylim(min1 - .2, max2 + .2)
		
	
	alpha = 0.01
	count1 = 0
	count2 = 0
	for chi in chisq1:
		p = chi2.pdf(chi[0][0], chi[0][1])
		if (p > alpha):
			count1 += 1.0 
		
		p = chi2.pdf(chi[1][0], chi[1][1])
		if (p > alpha):
			count2 += 1.0
			
	per1 = str(round(count1 / len(chisq1), 4) * 100) + ' %/ '
	per2 = str(round(count2 / len(chisq1), 4) * 100) + ' %/ '
	
	count1 = 0
	count2 = 0
	for chi in chisq0:
		p = chi2.pdf(chi[0][0], chi[0][1])
		if (p > alpha):
			count1 += 1.0 
		
		p = chi2.pdf(chi[1][0], chi[1][1])
		if (p > alpha):
			count2 += 1.0
			
	per1 += str(round(count1 / len(chisq0), 4) * 100) + ' %'
	per2 += str(round(count2 / len(chisq0), 4) * 100) + ' %'
		
	title = 'Probability of Variation/ Consistency: ' + per1
	ax1.text(-2, (math.ceil(max1 * 10) / 10) + .05, title, fontsize = 18, color = 'k', family = 'serif')
	title = 'Probability of Variation, Consistency: ' + per2
	ax2.text(-2, (math.ceil(max2 * 10) / 10) + .05, title, fontsize = 18, color = 'k', family = 'serif')
	ax1.set_xlim(-20, 125)
	ax2.set_xlim(-20, 125)
	
	fig.savefig(spec[0].dir + 'Pictures/' + spectra[0].sn +' _vs_' + spec[0].sn + '.png')
	fig.show()
	
	return chisq1, chisq0
	
	"""
	sig3, sig5 = [{}, {}] , [{}, {}]
	i = 0
	while(i < len(nspec)):
		j = i
		n = nspec[i]
		while(j < len(nspec)):
			s = nspec[j]
			if n.epoch != s.epoch:
				diff = round(np.abs(n.epoch - s.epoch), 3)
				change = np.abs(n.eqws[0] - s.eqws[0])
				sigma3 = 3 * np.sqrt(n.eqw_errs[0]**2 + s.eqw_errs[0]**2)
				sigma5 = 5 * np.sqrt(n.eqw_errs[0]**2 + s.eqw_errs[0]**2)
				if change > sigma3:
					add = 'y'
				else:
					add = 'n'
				try:
					sig3[0][diff].append(add)
				except:
					sig3[0][diff] = [add]
				if change > sigma5:
					add = 'y'
				else:
					add = 'n'
				try:
					sig5[0][diff].append(add)
				except:
					sig5[0][diff] = [add]
					
				change = np.abs(n.eqws[1] - s.eqws[1])
				sigma3 = 3 * np.sqrt(n.eqw_errs[1]**2 + s.eqw_errs[1]**2)
				sigma5 = 5 * np.sqrt(n.eqw_errs[1]**2 + s.eqw_errs[1]**2)
				if change > sigma3:
					add = 'y'
				else:
					add = 'n'
				try:
					sig3[1][diff].append(add)
				except:
					sig3[1][diff] = [add]
				if change > sigma5:
					add = 'y'
				else:
					add = 'n'
				try:
					sig5[1][diff].append(add)
				except:
					sig5[1][diff] = [add]
			j += 1
		i += 1

	for e in range(len(eqws[0])):
		plt.figure(e)
		plt.plot(epochs, eqws[:,e], color = 'k', linestyle = 'solid', linewidth = 2.5)
		plt.errorbar(eps, ews[:,e], yerr = errs[:,e], color = 'b', linestyle = 'none')
		plt.scatter(eps, ews[:,e], color = 'b')
		plt.show()

		
	temp = set(eps)
	temp = list(temp)
	temp.sort()
	ymin1 = min(eqws[0]) - .15
	ymin2 = min(eqws[1]) - .15
	i = 0
	while(i < len(temp)):
		j = i
		n = temp[i]
		while(j < len(temp)):
			s = temp[j]
			diff = round(np.abs(n - s), 3)
			if diff > 0:
				per3 = float(sig3[0][diff].count('y')) / float(len(sig3[0][diff]))
				per3 = str(int(round(per3, 2) * 100)) + '%'
				per5 = float(sig5[0][diff].count('y')) / float(len(sig5[0][diff]))
				per5 = str(int(round(per5, 2) * 100)) + '%'
				fig = plt.figure()
				fig.clf()
				plt.ylim(ymin1 - .08, 1.1)
				plt.plot(epochs, eqws[:,0], color = 'k', linestyle = 'solid', linewidth = 2.5)
				plt.errorbar(eps, ews[:,0], yerr = errs[:,0], color = 'b', linestyle = 'none')
				plt.scatter(eps, ews[:,0], color = 'b')
				plt.hlines(ymin1, n, s, color = 'k', linestyle = 'solid', linewidth = 1.5)
				plt.vlines(n, ymin1, ymin1 + .03, color = 'k', linestyle = 'solid', linewidth = 1.5)
				plt.vlines(s, ymin1, ymin1 + .03, color = 'k', linestyle = 'solid', linewidth = 1.5)
				plt.text(n - 6, ymin1 - .03, per3, fontsize = 14, color = 'b', family = 'serif')
				plt.text(n + 6, ymin1 - .03, per5, fontsize = 14, color = 'g', family = 'serif')
				plt.savefig(spec[0].dir + 'Pictures/' + spectra[0].sn + '_MgIIa_(' + str(round(n, 2)) + ',' + str(round(s, 2)) + ').png')
				plt.close()
				
				per3 = float(sig3[1][diff].count('y')) / float(len(sig3[1][diff]))
				per3 = str(int(round(per3, 2) * 100)) + '%'
				per5 = float(sig5[1][diff].count('y')) / float(len(sig5[1][diff]))
				per5 = str(int(round(per5, 2) * 100)) + '%'
				fig = plt.figure()
				fig.clf()
				plt.ylim(ymin2 - .08, 1.1)
				plt.plot(epochs, eqws[:,1], color = 'k', linestyle = 'solid', linewidth = 2.5)
				plt.errorbar(eps, ews[:,1], yerr = errs[:,1], color = 'b', linestyle = 'none')
				plt.scatter(eps, ews[:,1], color = 'b')
				plt.hlines(ymin1, n, s, color = 'k', linestyle = 'solid', linewidth = 1.5)
				plt.vlines(n, ymin2, ymin2 + .03, color = 'k', linestyle = 'solid', linewidth = 1.5)
				plt.vlines(s, ymin2, ymin2 + .03, color = 'k', linestyle = 'solid', linewidth = 1.5)
				plt.text(n - 6, ymin2 - .03, per3, fontsize = 14, color = 'b', family = 'serif')
				plt.text(n + 6, ymin2 - .03, per5, fontsize = 14, color = 'g', family = 'serif')
				plt.savefig(spec[0].dir + 'Pictures/' + spectra[0].sn + '_MgIIb_(' + str(round(n, 2)) + ',' + str(round(s, 2)) + ').png')
				plt.close()
			j += 1
		i += 1
	
	
	#return sig3, sig5, eps
	"""