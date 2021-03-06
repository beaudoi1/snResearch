import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from astropy.table import Table
import matplotlib.pyplot as plt
import glob
import sys

files = glob.glob("Data/data/*.flm")

def interpolate(file,bp):
	#read in file and seperate wavelength(x) and flux(y)
	data=Table.read(file,format='ascii')
	x=data["col1"]
	y=data["col2"]

	#remove the absorption lines from the spectrum
	nolines=np.where((x<2341) | ((x>2352) & (x<2372)) | ((x>2390) & (x<2584)) | ((x>2607) & (x<2791)) | ((x>2813) & (x<2850)) | (x>2861))
	#nolines=np.where((x<2341) | ((x>2351) & (x<2373)) | ((x>2380) & (x<2381)) | ((x>2390) & (x<2582)) | ((x>2593) & (x<2597)) | ((x>2607) & (x<2791)) | ((x>2812) & (x<2850)) | (x>2860))

	
	#creates b-spline over spectrum
	sm1=inter.splrep(x[nolines],y[nolines])
	xx=linspace(min(x),max(x),100)
	nolinesxx=np.where((xx<2341) | ((xx>2353) & (xx<2372)) | ((xx>2390) & (xx<2584)) | ((xx>2607) & (xx<2791)) | ((xx>2813) & (xx<2850)) | (xx>2861))
	xx=xx[nolinesxx]

	#breakpoints in spectrum
	aa=bp
	bb=np.concatenate((xx,aa)) #add breakpoints to array
	bb.sort()

	#creates and fits b-spline over wavelength range
	y1=inter.splev(bb,sm1)
	sm2=inter.splrep(bb,y1)
	y2=inter.splev(x,sm2)

	#creates residual for spectrum
	norm=y/y2
	
	xadd=[max(x)+1,max(x)+2]
	normadd=[1,1]
	new_x=np.concatenate((x,xadd))
	new_norm=np.concatenate((norm,normadd))
	return x,y,y2, new_x,new_norm
	
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
	
def find_errors(file,sflux):
	data=Table.read(file,format='ascii')
	w=data["col1"]
	f=data["col2"]
	
	sf=smooth(f)
	res=f-sf
	ares=np.abs(res)
	se = smooth(ares)
	norm_e = se/sflux
	erradd=[0,0]
	new_err=np.concatenate((norm_e,erradd))
	return new_err
	
sp_names=[]
for i in range(len(files)):
	name=''
	file = files[i]
	name = file[19:]
	name = name[:-4]
	sp_names.append(name)
waves=[]
res_waves=[]
fluxes=[]
bsplines=[]
residuals=[]
errs=[]

num=len(files)

for i in range(num):
	k=1
	sflux=0
	while (k==1):
		points=input('Breakpoints for ' + sp_names[i] + ': ')
		bps=[]
		for j in range(len(points)):
			bps.append(float(points[j]))
		wave,old_flux,sp_flux,res_wave,res=interpolate(files[i],bps)
		"""plot1,=plt.plot(wave,old_flux)
		plot2,=plt.plot(wave,sp_flux)
		plt.show()
		plot1,=plt.plot(res_wave,res)
		plt.ylim(0,1.5)
		plt.show()"""
		cont = raw_input('More breakpoints?: (y/n) ')
		if(cont=='y'):
			continue
		else:
			waves.append(wave)
			res_waves.append(res_wave)
			fluxes.append(old_flux)
			bsplines.append(sp_flux)
			residuals.append(res)
			sflux=sp_flux
			k=2	
	err=find_errors(files[i],sflux)
	errs.append(err)
	
for i in range(num):
	save_to = "Data/Pictures/"+ sp_names[i]+ ".png"
	plt.figure(i)
	plt.title(sp_names[i])
	plt.subplot(211)
	plt.plot(waves[i],fluxes[i])
	plt.plot(waves[i],bsplines[i])
	plt.subplot(212)
	plt.errorbar(res_waves[i],residuals[i],yerr=errs[i])
	plt.ylim(0,1.5)
	plt.savefig(save_to)
	
for i in range(num):
	t=Table([res_waves[i],residuals[i],errs[i]],names=('col1','col2','col3'))
	write_to = "Data/residuals/" + sp_names[i] + ".flm"
	t.write(write_to,format='ascii')



