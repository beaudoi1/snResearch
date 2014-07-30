import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from astropy.table import Table
import matplotlib.pyplot as plt
import glob
import sys
import pyspeckit
import math

files = glob.glob("Data/data/*.flm") #Data/data/test_data/*.flm

FeII_obs = [2346.17,2376.62,2384.9]
MnII_obs = [2588.77,2601.14]
MgII_obs = [2798.7,2805.98]
MgI_obs = [2854.97]

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
	return new_x,new_norm
	
def rand_break(breakpoints):
	bkpts=[]
	for b in breakpoints:
		bkpts.append(b)
	febkpts,mnbkpts,mg2bkpts,mg1bkpts,other=[],[],[],[],[]
	for b in bkpts:
		if ((b>2331)&(b<2397)):
			febkpts.append(b)
		elif ((b>2578)&(b<2614)):
			mnbkpts.append(b)
		elif ((b>2783)&(b<2820)):
			mg2bkpts.append(b)
		elif ((b>2842)&(b<2865)):
			mg1bkpts.append(b)
		else:
			other.append(b)
	new_bkpts=[]
	feshift=np.random.randn(len(febkpts)) * 3
	mnshift=np.random.randn(len(mnbkpts)) * 4.25
	mg2shift=np.random.randn(len(mg2bkpts)) * 3
	mg1shift=np.random.randn(len(mg1bkpts)) * 6
	oshift=np.random.randn(len(other)) * 1.5
	new_fe=feshift+febkpts
	new_mn=mnshift+mnbkpts
	new_mg2=mg2shift+mg2bkpts
	new_mg1=mg1shift+mg1bkpts
	new_other=oshift+other
	for fe in new_fe:
		new_bkpts.append(fe)
	for mn in new_mn:
		new_bkpts.append(mn)
	for mg2 in new_mg2:
		new_bkpts.append(mg2)
	for mg1 in new_mg1:
		new_bkpts.append(mg1)
	for o in new_other:
		new_bkpts.append(o)
	new_bkpts.sort()
	return new_bkpts
	
def find_eqw(x,y,guesses,i):
	#Creates a spectrum from the residual
	sp=pyspeckit.Spectrum(xarr=x,data=y)

	#parts of the continuum to exclude (absorption lines)
	exclude=[2334.97,2355.38,2369.27,2379.75,2380.51,2390.63,2581.73,2592,2597,2606.61,2792.1,2801.51,2802.49,2811.56,2848.32,2860.29]

	#guess of height, rest wavelength, width of each absorption line
	guesses=guesses
	
	#fits the baseline to the residual
	xmin=max(x)-2
	xmax=max(x)
	sp.baseline(xmin=xmin,xmax=xmax,order=1,exclude=exclude,subtract=False,reset_selection=True)
	#sp.baseline(xmin=2331,xmax=2900,order=1,exclude=exclude,subtract=False,reset_selection=True)

	
	#Fits a gaussian to each line (height, rest wavelength, width)
	sp.plotter(xmin=2200,xmax=3100,ymin=0,ymax=1.5,figure=1)
	sp.specfit(guesses=guesses,annotate=False)
	sp.plotter.savefig('Data/eqw/MonteCarlo/lines/' +sp_names[a]+',rand ' + str(i) +'.png')
	

	#Gives the EW of each line
	eq=sp.specfit.EQW(components=True)
	
	return eq

def find_err(eqws,ion,name):
	eqws.sort()
	sort_eqw = eqws
	n=len(sort_eqw)

	"""
	figure out how to plot histograms
	plot median, lo and hi error, and best from before
	"""
	err_lo = sort_eqw[int(round(.16*n))]
	err_hi = sort_eqw[int(round(.84*n))]
	std = np.std(sort_eqw)
	med = np.median(eqws)
	
	plt.figure()
	hist,bins=np.histogram(sort_eqw,bins=len(sort_eqw)/2)
	width=.7*(bins[1]-bins[0])
	center=(bins[:-1]+bins[1:])/2
	plt.bar(center,hist,align='center',width=width)
	plt.plot([err_lo],[5], 'go', label='err_lo')
	plt.plot([err_hi],[5],'co',label='err_hi')
	plt.plot([med],[5],'ko', label='median')
	plt.plot([ion],[5],'ro', label='original')
	legend=plt.legend(loc='upper right')
	#plt.xlim(2.27256,2.272565)
	#plt.ylim(0,15)
	plt.savefig('Data/eqw/MonteCarlo/Histograms/MCHist' + sp_names[a] + '-' + name + '.png',orientation='landscape')
	#plt.show()
	
	return med-err_lo,err_hi-med, std
	"""
	find chi square, lmfit?
	"""
	
	
sp_names=[]
for i in range(len(files)):
	name=''
	file = files[i]
	name = file[19:]
	#name = name[10:]
	name = name[:-4]
	sp_names.append(name)

num=len(files)
bkpts_files=glob.glob('Data/breakpoints/*.flm') #Data/breakpoints/test_bkpts/*.flm
breakpoints=[]
for i in range(len(bkpts_files)):
	bkpts=Table.read(bkpts_files[i],format='ascii')
	breakpoints.append(bkpts["col1"])
	
print "breakpoints read in"

for a in range(num):
	file=files[a]
	print "Finding error for", sp_names[a] 
	mc=500
	res_waves=[]
	residuals=[]	
	for i in range(mc):
		bps=rand_break(breakpoints[a])
		res_wave,res=interpolate(file,bps)
		res_waves.append(res_wave)
		residuals.append(res)
		print i+1, "set of breakpoints created for", sp_names[a]

	eqws=[]	
	for i in range(len(res_waves)):
		guesses=[-.41,FeII_obs[0],1.65,-.246,FeII_obs[1],1.31,-.518,FeII_obs[2],1.58,-.38,MnII_obs[0],1.6,-.484,MnII_obs[1],2.07,-.696,MgII_obs[0],1.699,-.638,MgII_obs[1],1.6,-.161,MgI_obs[0],1.62]
		eqws.append(find_eqw(res_waves[i],residuals[i],guesses,i))
		print i+1, "eqws found"
	
	fe2a=[]
	fe2b=[]
	fe2c=[]
	mn2a=[]
	mn2b=[]
	mg2a=[]
	mg2b=[]
	mg1a=[]
	for i in range(len(eqws)):
		set=eqws[i]
		fe2a.append(set[0])
		fe2b.append(set[1])
		fe2c.append(set[2])
		mn2a.append(set[3])
		mn2b.append(set[4])
		mg2a.append(set[5])
		mg2b.append(set[6])
		mg1a.append(set[7])
	
	eqw_file='Data/eqw/'+sp_names[a]+'.dat'
	lo_errors=[]
	hi_errors=[]
	stds=[]
	fe2a_lo,fe2a_hi,fe2a_std=find_err(fe2a,Table.read(eqw_file,format='ascii')["col2"][0],'FeIIa')
	lo_errors.append(fe2a_lo)
	hi_errors.append(fe2a_hi)
	stds.append(fe2a_std)
	fe2b_lo,fe2b_hi,fe2b_std=find_err(fe2b,Table.read(eqw_file,format='ascii')["col2"][1],'FeIIb')
	lo_errors.append(fe2b_lo)
	hi_errors.append(fe2b_hi)
	stds.append(fe2b_std)
	fe2c_lo,fe2c_hi,fe2c_std=find_err(fe2c,Table.read(eqw_file,format='ascii')["col2"][2],'FeIIc')
	lo_errors.append(fe2c_lo)
	hi_errors.append(fe2c_hi)
	stds.append(fe2c_std)
	mn2a_lo,mn2a_hi,mn2a_std=find_err(mn2a,Table.read(eqw_file,format='ascii')["col2"][3],'mn2a')
	lo_errors.append(mn2a_lo)
	hi_errors.append(mn2a_hi)
	stds.append(mn2a_std)
	mn2b_lo,mn2b_hi,mn2b_std=find_err(mn2b,Table.read(eqw_file,format='ascii')["col2"][4],'mn2b')
	lo_errors.append(mn2b_lo)
	hi_errors.append(mn2b_hi)
	stds.append(mn2b_std)
	mg2a_lo,mg2a_hi,mg2a_std=find_err(mg2a,Table.read(eqw_file,format='ascii')["col2"][5],'mg2a')
	lo_errors.append(mg2a_lo)
	hi_errors.append(mg2a_hi)
	stds.append(mg2a_std)
	mg2b_lo,mg2b_hi,mg2b_std=find_err(mg2b,Table.read(eqw_file,format='ascii')["col2"][6],'mg2b')
	lo_errors.append(mg2b_lo)
	hi_errors.append(mg2b_hi)
	stds.append(mg2b_std)
	mg1_lo,mg1_hi,mg1_std=find_err(mg1a,Table.read(eqw_file,format='ascii')["col2"][7],'mg1a')
	lo_errors.append(mg1_lo)
	hi_errors.append(mg1_hi)
	stds.append(mg1_std)
	
	
	#find error for each ion
	print "errors found for ", sp_names[a]
	
	ions=['FeIIa','FeIIb','FeIIc','MnIIa','MnIIb','MgIIa','MgIIb','MgIa']	
	for i in range(len(lo_errors)):
		t=Table([ions,lo_errors,hi_errors,stds],names=('col1','col2','col3','col4'))
		t.write('Data/eqw/MonteCarlo/error_files/' + sp_names[a] + '.dat',format='ascii')
	print "file written for ", sp_names[a]
	
#print "All eqws found, arranging by epoch"
	
print "done"