from pyspeckit import mpfit
from astropy.table import Table
import matplotlib.pylab as plt
import numpy as np
import math
from scipy import optimize
from scipy import integrate
import glob

FeII_obs = [2346.17,2376.62,2384.9]
MnII_obs = [2588.77,2601.14]
MgII_obs = [2798.7,2805.98]
MgI_obs = [2854.97]

include=[2334.97,2355.38,2369.27,2379.75,2380.51,2390.63,2581.73,2592,2597,2606.61,2792.1,2801.51,2802.49,2811.56,2848.32,2860.29]
guesses = [-.41,FeII_obs[0],1.65,-.246,FeII_obs[1],1.31,-.518,FeII_obs[2],1.58,-.38,MnII_obs[0],1.6,-.484,MnII_obs[1],2.07,-.696,MgII_obs[0],1.699,-.638,MgII_obs[1],1.6,-.161,MgI_obs[0],1.62]

gaussian = lambda p, x: 1+p[0]*np.exp(-(x-p[1])**2/(2*p[2]**2))
residual = lambda p, x, y: y - gaussian(p, x)

def find_eqw(file):
	data=Table.read(file,format='ascii')
	w=data["col1"]
	f=data["col2"]
	e=data["col3"]
		
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
	while(g<len(guesses)):
		inc = np.where((w>include[i]) & (w<include[i+1]))
		wave = w[inc]
		flux = f[inc]
		errs = e[inc]
		
		x = np.array(wave)
		y = np.array(flux)
		err = np.array(errs)

		p = [guesses[g],guesses[g+1],guesses[g+2]]
		p1, pcov, infodict, errmsg, success= optimize.leastsq(residual, p[:], args=(x,y), full_output = 1)
		gaus = gaussian(p1, x)	
		
		s_sq = (residual(p1, x, y)**2).sum()/(len(y)-len(p))
		pcov = pcov * s_sq		
		
		errors = []
		for j in range(len(p1)):
			errors.append(np.abs(pcov[j][j])**.5)
			
		errors = np.array(errors)
		
		error = np.sqrt(np.sqrt(math.pi) * ((errors[0]*p1[2])**2 + (errors[2]*p1[0])**2))
		
		for j in range(len(gaus)):
			baseline_x.append(x[j])
			baseline_y.append(gaus[j])
		
		if(i<14):
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

	for e in range(len(end[0])):
		baseline_x.append(w[end[0][e]])
		baseline_y.append(1)
	
	plt.figure()
	plt.plot(baseline_x,baseline_y,label='gaus')
	plt.plot(w,f,label='cont')
	legend = plt.legend(loc='upper right')
	#plt.show()
	
	return eqws, eqw_errs
	
files = glob.glob('Data/residuals/*.flm')
res_names = []
eqws = []
errors = []
for k in range(len(files)):
	name = files[k][15:]
	name = name[:-4]
	res_names.append(name)
	eqw, error = find_eqw(files[k])
	eqws.append(eqw)
	errors.append(error)
	
ions = ['FeIIa','FeIIb','FeIIc','MnIIa','MnIIb','MgIIa','MgIIb','MgIa']
for k in range(len(res_names)):
	t = Table([ions,eqws[k],errors[k]], names = ('col1', 'col2', 'col3'))
	t.write('Data/eqw/' + res_names[k] + '.dat', format='ascii')

