import spAnalyze as a
a.pickSN()
spectra = a.loadData('sn2013dy/Data')
for s in spectra:
	a.getBkpts(s)
guess6 = a.guesses
guess10 = [-.3, 2797.1, 3.4, -.45, 2807.77, 7.4, -.34, 2815.21, 4.2, -.2, 2854.56, 2.94, -.189, 2864.77, 2.61]
a.exclude = [2789.67, 2826.71, 2843.33, 2878.64]

ion10 = a.ions
ion6 = ion10[:4]
for s in spectra:
	if(len(s.inc) < 7):
		a.guesses = guess6
		a.ions = ion6
	else:
		a.guesses = guess10
		a.ions = ion10
	a.findEqw(s, True)
	xlim(2750, 2900)
	ylim(0, 1.5)
	print s.eqws
	a.sepMgDoub(s)
	print s.eqws

rests = [a.rests[len(a.rests) - 1], a.rests[len(a.rests) - 1]]
for s in spectra:
		if((s.visit == 6) | (s.visit == 9)):
			temp = []
			for r in rests:
				a.sigma3(s, r)
				temp.append(s.sig3)
			s.sig3 = temp
for s in spectra:
	if ((s.visit == 6) | (s.visit == 9)):
		s.eqws = np.append(s.eqws, s.sig3)
		print s.eqws

for s in spectra:
	if ((s.visit == 6) | (s.visit == 9)):
		s.eqw_errs.append(0.00)
		s.eqw_errs.append(0.00)
		
rests = a.rests[:5]

for s in spectra:
	temp = []
	for r in rests:
		a.sigma3(s, r)
		temp.append(s.sig3)
	s.sig3 = temp
	
for s in spectra:
	if ((s.visit == 6) | (s.visit == 9)):
		a.guesses = guess6
		a.ions = ion6
	else:
		a.guesses = guess10
		a.ions = ion10
	a.monteCarlo(500, s, True)
	
spectra[6].mc_lo.append(0.00)
spectra[6].mc_lo.append(0.00)
spectra[6].mc_hi.append(0.00)
spectra[6].mc_hi.append(0.00)
spectra[9].mc_lo.append(0.00)
spectra[9].mc_lo.append(0.00)
spectra[9].mc_hi.append(0.00)
spectra[9].mc_hi.append(0.00)

a.ions = ['MgIIa', 'MgIIb', 'MgIIc', 'MgIId', 'MgIa', 'MgIb']
