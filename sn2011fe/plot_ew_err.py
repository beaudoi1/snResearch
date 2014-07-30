import matplotlib.pyplot as plt
import glob
from astropy.table import Table
import numpy as np

def chi_square(epochs,ion,std):
	func=np.polyfit(epochs,ion,0)
	chi=0
	for i in range(len(exp)):
		chi+=(ion[i]-func)**2/std[i]**2
	return chi

ew_files=glob.glob("Data/eqw/*.dat")
err_files=glob.glob("Data/eqw/MonteCarlo/error_files/*.dat")

names=[]
lines=[]
eqws=[]

for i in range(len(ew_files)):
	file=ew_files[i]
	name=file[:-4]
	name=name[9:]
	names.append(name)
	eqw=Table.read(ew_files[i],format='ascii')
	lines=eqw["col1"]
	eqws.append(eqw["col2"])
	
eqw_epoch=[]
eqw_epoch.append(eqws[0])
eqw_epoch.append(eqws[10])
eqw_epoch.append(eqws[1])
eqw_epoch.append(eqws[2])
eqw_epoch.append(eqws[3])
eqw_epoch.append(eqws[4])
eqw_epoch.append(eqws[5])
eqw_epoch.append(eqws[6])
eqw_epoch.append(eqws[7])
eqw_epoch.append(eqws[8])
eqw_epoch.append(eqws[9])

epochs=[-13.110599,-10.041886,-6.9081770,-2.9512300,0.039233502,3.2484383,9.1496879,9.1496879,20.687449,26.685855,40.447814]
ions=[]
for j in range(len(lines)):
	ion=[]
	for i in range(len(eqw_epoch)):
		eqw=eqw_epoch[i]
		ion.append(eqw[j])
	ions.append(ion)
	
deltaion=[]
for i in range(len(ions)):
	dion=[]
	for j in range(len(ions[i])):
		dion.append(ions[i][j]-ions[i][4])
	deltaion.append(dion)
	
errs_lo=[]
errs_hi=[]
errs_stds=[]
for i in range(len(err_files)):
	errs=Table.read(err_files[i],format='ascii')
	errs_lo.append(errs["col2"])
	errs_hi.append(errs["col3"])
	errs_stds.append(errs["col4"])
	
eqw_lo_errs=[]
eqw_hi_errs=[]
eqw_stds=[]
eqw_lo_errs.append(errs_lo[0])
eqw_hi_errs.append(errs_hi[0])
eqw_stds.append(errs_stds[0])
eqw_lo_errs.append(errs_lo[10])
eqw_hi_errs.append(errs_hi[10])
eqw_stds.append(errs_stds[10])
eqw_lo_errs.append(errs_lo[1])
eqw_hi_errs.append(errs_hi[1])
eqw_stds.append(errs_stds[1])
eqw_lo_errs.append(errs_lo[2])
eqw_hi_errs.append(errs_hi[2])
eqw_stds.append(errs_stds[2])
eqw_lo_errs.append(errs_lo[3])
eqw_hi_errs.append(errs_hi[3])
eqw_stds.append(errs_stds[3])
eqw_lo_errs.append(errs_lo[4])
eqw_hi_errs.append(errs_hi[4])
eqw_stds.append(errs_stds[4])
eqw_lo_errs.append(errs_lo[5])
eqw_hi_errs.append(errs_hi[5])
eqw_stds.append(errs_stds[5])
eqw_lo_errs.append(errs_lo[6])
eqw_hi_errs.append(errs_hi[6])
eqw_stds.append(errs_stds[6])
eqw_lo_errs.append(errs_lo[7])
eqw_hi_errs.append(errs_hi[7])
eqw_stds.append(errs_stds[7])
eqw_lo_errs.append(errs_lo[8])
eqw_hi_errs.append(errs_hi[8])
eqw_stds.append(errs_stds[8])
eqw_lo_errs.append(errs_lo[9])
eqw_hi_errs.append(errs_hi[9])
eqw_stds.append(errs_stds[9])
	
ion_lo_errs=[]
ion_hi_errs=[]
for j in range(len(lines)):
	ion_lo_err=[]
	ion_hi_err=[]
	for i in range(len(eqw_lo_errs)):
		lo_errs=eqw_lo_errs[i]
		hi_errs=eqw_hi_errs[i]
		ion_lo_err.append(lo_errs[j])
		ion_hi_err.append(hi_errs[j])
	ion_lo_errs.append(ion_lo_err)
	ion_hi_errs.append(ion_hi_err)
	
a=[epochs[0],epochs[10]]
b=[0,0]

chi_squares=[]
for i in range(len(ions)):
	chi=chi_square(epochs,ions[i],eqw_stds[i])
	print 'Chi Square for visit',i+1,'=',chi
	chi_squares.append(chi)
	
plt.figure(1)
plt.subplot(311)
plt.title('FeII')
plt.errorbar(epochs,ions[0],yerr=[ion_lo_errs[0],ion_hi_errs[0]],label=lines[0],color='k')
plt.errorbar(epochs,ions[1],yerr=[ion_lo_errs[1],ion_hi_errs[1]],label=lines[1],color='b')
plt.errorbar(epochs,ions[2],yerr=[ion_lo_errs[2],ion_hi_errs[2]],label=lines[2],color='c')
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylabel('EQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplot(312)
plt.title('MnII')
plt.errorbar(epochs,ions[3],yerr=[ion_lo_errs[3],ion_hi_errs[3]],label=lines[3],color='k')
plt.errorbar(epochs,ions[4],yerr=[ion_lo_errs[4],ion_hi_errs[4]],label=lines[4],color='b')
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylim(0,4)
plt.ylabel('EQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplot(313)
plt.title('MgII & MgI')
plt.errorbar(epochs,ions[5],yerr=[ion_lo_errs[5],ion_hi_errs[5]],label=lines[5],color='k')
plt.errorbar(epochs,ions[6],yerr=[ion_lo_errs[6],ion_hi_errs[6]],label=lines[6],color='b')
plt.errorbar(epochs,ions[7],yerr=[ion_lo_errs[7],ion_hi_errs[7]],label=lines[7],color='c')
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylim(0,4)
plt.xlabel('Days since B maximum')
plt.ylabel('EQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplots_adjust(hspace=.35,right=.8)
plt.savefig('Data/eqw/EQW by Epoch with Error,ion.png',orientation='landscape')
plt.show()

plt.figure(2)
for i in range(len(ions)):
	plt.errorbar(epochs,ions[i],yerr=[ion_lo_errs[i],ion_hi_errs[i]],label=lines[i])
legend=plt.legend(loc='upper right')
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.savefig('Data/eqw/EQW by Epoch with Error,all.png',orientation='landscape')
plt.show()
	
plt.figure(3)
plt.subplot(311)
plt.title('FeII')
plt.errorbar(epochs,deltaion[0],yerr=[ion_lo_errs[0],ion_hi_errs[0]],label=lines[0],color='k')
plt.errorbar(epochs,deltaion[1],yerr=[ion_lo_errs[1],ion_hi_errs[1]],label=lines[1],color='b')
plt.errorbar(epochs,deltaion[2],yerr=[ion_lo_errs[2],ion_hi_errs[2]],label=lines[2],color='c')
plt.plot(a,b,'g-')
#plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylabel('deltaEQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplot(312)
plt.title('MnII')
plt.errorbar(epochs,deltaion[3],yerr=[ion_lo_errs[3],ion_hi_errs[3]],label=lines[3],color='k')
plt.errorbar(epochs,deltaion[4],yerr=[ion_lo_errs[4],ion_hi_errs[4]],label=lines[4],color='b')
plt.plot(a,b,'g-')
#plt.xlim(epochs[0],epochs[len(epochs)-1])
#plt.ylim(0,4)
plt.ylabel('deltaEQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplot(313)
plt.title('MgII & MgI')
plt.errorbar(epochs,deltaion[5],yerr=[ion_lo_errs[5],ion_hi_errs[5]],label=lines[5],color='k')
plt.errorbar(epochs,deltaion[6],yerr=[ion_lo_errs[6],ion_hi_errs[6]],label=lines[6],color='b')
plt.errorbar(epochs,deltaion[7],yerr=[ion_lo_errs[7],ion_hi_errs[7]],label=lines[7],color='c')
plt.plot(a,b,'g-')
#plt.xlim(epochs[0],epochs[len(epochs)-1])
#plt.ylim(0,4)
plt.xlabel('Days since B maximum')
plt.ylabel('deltaEQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplots_adjust(hspace=.35,right=.8)
plt.savefig('Data/eqw/deltaEQW by Epoch with Error,ion.png',orientation='landscape')
plt.show()

plt.figure(4)
for i in range(len(ions)):
	plt.errorbar(epochs,deltaion[i],yerr=[ion_lo_errs[i],ion_hi_errs[i]],label=lines[i])
plt.plot(a,b)
legend=plt.legend(loc='upper right')
#plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.savefig('Data/eqw/deltaEQW by Epoch with Error,all.png',orientation='landscape')
plt.show()

	
	
	
	
	