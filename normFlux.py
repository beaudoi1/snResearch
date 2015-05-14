import spAnalyze as a
import glob
import matplotlib.pyplot as plt
import sys
import numpy as np
import plot


files = glob.glob('norm_flux/*.flm')
spectra = []

for f in files:
	sp = a.spectrum()
	name = f[10:]
	name = name[:-8]
	sp.name = name
	sp.sn = name[2:]
	
	sp.dir = 'norm_flux/' + sp.sn[2:] + '/'
	
	sp.bkpts = []
	a.getSpectrum(f, sp)
	
	spectra.append(sp)

for s in spectra:
	s.norm_flux = s.flux
	print s.name
	a.pickSN()
	plt.figure()
	plt.plot(s.wave, s.norm_flux)
	plt.ylim(0, 1.5)
	plt.xlim(2300, 3000)
	try:
		a.loadInc(s)
		print s.inc
		cont = True
		while(cont):
			temp = raw_input('Use these include values? (y/n): ')
			if (temp == 'y'):
				cont = False
			elif (temp == 'n'):
				a.setInc(s)
				cont = False
			else:
				cont = True
	except:
		print 'No include values'
		a.setInc(s)
	a.saveInc(s)
	plt.close()
	s.cont_err = np.zeros(len(s.wave))
	if (len(s.inc) > 0):
		lines = True
	else:
		lines = False
	if(lines):	
		a.findEqw(s, True)
		plt.xlim(2300, 3000)
		plt.ylim(0, 1.5)
		temp = raw_input('Press any key to continue: ')
		plt.close()
		cont = True
		while(cont):
			trip = raw_input('Is there a Mg II triplet? (y/n): ')
			if trip == 'y':
				a.sepMgDoub(s)
				a.monteCarlo(100, s, True)
				cont = False
			elif trip == 'n':
				a.monteCarlo(100, s, False)
				cont = False
			else:
				cont = True
	limIons = input('Limit ions: ')
	if(len(limIons) > 0):
		rests = input('Limit rest wavelengths: ')
		temp = []
		for r in rests:
			a.sigma3(s, r)
			temp.append(s.sig3)
		s.sig3 = temp
	if ((len(s.inc) > 0) & (len(limIons) > 0)):
		zeros = np.zeros(len(s.eqws))
		s.sig3 = np.append(s.sig3, zeros)
		zeros = np.zeros(len(limIons))
		s.eqws = np.append(zeros, s.eqws)
		s.eqw_errs = np.append(zeros, s.eqw_errs)
		s.mc_lo = np.append(zeros, s.mc_lo)
		s.mc_hi = np.append(zeros, s.mc_hi)
		a.ions = np.append(limIons, a.ions)
		a.calcErr(s)
	
	a.saveData(s)
			
a.plotEqw(spectra)
			
"""
['FeIIa','FeIIb','FeIIc','FeIId','FeIIe','MgIIa','MgIIb','MgIa']
[2344.2128, 2374.4601, 2382.7641, 2586.6494, 2600.1722, 2796.354, 2803.5311, 2852.9628]
['FeIIa','FeIIb','FeIIc','FeIId','FeIIe']
[2344.2128, 2374.4601, 2382.7641, 2586.6494, 2600.1722]
['FeIIa','FeIIb','FeIIc','FeIId','FeIIe','MgIIa','MgIIb','MgIa']
[2344.2128, 2374.4601, 2382.7641, 2586.6494, 2600.1722, 2796.354, 2803.5311, 2852.9628]
[]
['FeIIa','FeIIb','FeIIc','FeIId','FeIIe','MgIIa','MgIIb','MgIa']
[2344.2128, 2374.4601, 2382.7641, 2586.6494, 2600.1722, 2796.354, 2803.5311, 2852.9628]
['FeIIa','FeIIb','FeIIc','FeIId','FeIIe']
[2344.2128, 2374.4601, 2382.7641, 2586.6494, 2600.1722]
['FeIIa','FeIIb','FeIIc','FeIId','FeIIe']
[2344.2128, 2374.4601, 2382.7641, 2586.6494, 2600.1722]
['FeIIa','FeIIb','FeIIc','FeIId','FeIIe']
[2344.2128, 2374.4601, 2382.7641, 2586.6494, 2600.1722]
"""

	

	
