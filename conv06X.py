import spAnalyze as a
import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve as con
from scipy.interpolate import interp1d

files = glob.glob('sn2006X/Data/data/*.dat')
spectra = []
a.pickSN('06x')

for f in files:
	spectra.append(a.initSpec(f))
	
def convolve(sp, res = 1.3, rebin = 1.5):
	pix  = np.round(sp.wave[1] - sp.wave[0], 4)
	num = (3000 - 2600) / pix
	num = int(np.round(num, 0))
	swave = np.arange(2600, 2600 + num * pix, pix)
	sflux = []
	for s in swave:
		sflux.append(1)
	
	sp.norm_flux = sp.flux
	plt.figure()
	plt.plot(sp.wave, sp.norm_flux)
	cont = True
	while(cont):
		try:
			temp = raw_input('Select wavelength of center of Na line: ')
			temp = float(temp)
			cont = False
		except NameError:
			print 'Invalid'
	plt.close()
	sp.wave = sp.wave / (temp / 2803.5311)
	plt.figure()
	plt.plot(sp.wave, sp.norm_flux)
	plt.ylim(0, 1.5)
	cont = True
	while(cont):
		try:
			cent = input('Select wavelength range of Na line: ')
			cont = False
		except NameError:
			print 'Invalid'
	plt.close()
	inc = np.where((sp.wave > cent[0]) & (sp.wave < cent[1]))[0]
	na = sp.flux[inc]
	
	num = len(na)
	start = inc[num / 2] - inc[0]
	
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
			
	sa = ia - start
	sb = ib - start
	for n in na:
		sflux[sa] = n
		sflux[sb] = n
		sa += 1
		sb += 1
		
	pix_res = 0
	pix_rebin = 0
	for w in range(len(swave)):
		diff = swave[w] - swave[0]
		if diff < res:
			pix_res += 1
		if diff < rebin:
			pix_rebin += 1
			
		
	ker = Gaussian1DKernel(pix_res)
	conv = con(sflux, ker, boundary = 'extend')
		
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
		
	s = a.spectrum()
	s.name = sp.name + '_Conv' 
	s.visit = sp.visit
	s.sn = sp.sn
	s.dir = sp.dir
	s.epoch = sp.epoch
	s.wave = np.array(rw)
	s.norm_flux = np.array(rf)
	return s	
	
def monteCarlo(sp1, sp2, s_rms1, s_rms2, mc = 20): 
	a.loadInc(sp1)
	a.loadInc(sp2)
	s_rms = s_rms1
	inc = np.where((s_rms.wave > 2650) & (s_rms.wave < 2750))[0]
	res = 1 - s_rms.norm_flux[inc]
	temp = []
	for r in res:
		temp.append(r**2)
	rms1 = np.sqrt(np.sum(temp) / len(res))
	
	s_rms = s_rms2
	inc = np.where((s_rms.wave > 2650) & (s_rms.wave < 2750))[0]
	res = 1 - s_rms.norm_flux[inc]
	temp = []
	for r in res:
		temp.append(r**2)
	rms2 = np.sqrt(np.sum(temp) / len(res))
		
	ews1 = []
	errs1 = []
	ews2 = []
	errs2 = []
	for m in range(mc):
		s1 = a.spectrum()
		s2 = a.spectrum()
		s1.wave = sp1.wave
		s1.inc = sp1.inc
		s1.norm_flux = []
		s1.sn = sp1.sn
		s2.wave = sp2.wave
		s2.inc = sp2.inc
		s2.norm_flux = []
		s2.sn = sp2.sn
		for f in range(len(sp1.norm_flux)):
			s1.norm_flux.append(sp1.norm_flux[f] + np.random.normal(0, .68, 1)[0] * rms1)
		for g in range(len(sp2.norm_flux)):
			s2.norm_flux.append(sp2.norm_flux[g] + np.random.normal(0, .68, 1)[0] * rms2)
			
		s1.norm_flux = np.array(s1.norm_flux)
		s2.norm_flux = np.array(s2.norm_flux)
		
		a.findEqw(s1)
		a.findEqw(s2)
		ews1.append(s1.eqws)
		errs1.append(s1.eqw_errs)
		ews2.append(s2.eqws)
		errs2.append(s2.eqw_errs)
		
	counta = 0
	countb = 0
	for i in range(len(ews1)):
		change = ews2[i][0] - ews1[i][0]
		sigma3 = np.sqrt(errs2[i][0]**2 + errs1[i][0]**2)
		if change > sigma3:
			counta += 1
			
		change = ews2[i][1] - ews1[i][1]
		sigma3 = np.sqrt(errs2[i][1]**2 + errs1[i][1]**2)
		if change > sigma3:
			countb += 1	
		
	return ([counta, countb])
	
def findSpec(spec):
	diff1 = 9999
	diff2 = 9999
	
	for s in spec:
		if np.abs(s.epoch - spectra[0].epoch) < diff1:
			diff1 = np.abs(s.epoch - spectra[0].epoch)
			sp1 = s
		if np.abs(s.epoch - spectra[1].epoch) < diff2:
			diff2 = np.abs(s.epoch - spectra[1].epoch)
			sp2 = s
			
	return [sp1, sp2]
	
def linFit(nspec):
	ews1 = []
	errs1 = []
	ews2 = []
	errs2 = []
	epochs = []
	
	i = 0
	cont = True
	while(cont):
		n = nspec[i]
		ews1.append(n.eqws[0])
		ews2.append(n.eqws[1])
		errs1.append(n.eqw_errs[0])
		errs2.append(n.eqw_errs[1])
		epochs.append(n.epoch)
		i += 1
		if (i == 3):
			cont = False
	
	
	fit1 = interp1d(x = epochs, y = ews1)
	fit2 = np.polyfit(x = epochs, y = ews2)
	