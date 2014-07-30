import matplotlib.pyplot as plt
import glob
from astropy.table import Table
import numpy as np

files=glob.glob("Data/eqw/*.dat")

names=[]
lines=[]
eqws=[]
errors=[]

for i in range(len(files)):
	file=files[i]
	name=file[:-4]
	name=name[9:]
	names.append(name)
	eqw=Table.read(files[i],format='ascii')
	lines=eqw["col1"]
	eqws.append(eqw["col2"])
	errors.append(eqw["col3"])
	
print errors[0]	
	
eqw_epoch=[]
err_epoch=[]
eqw_epoch.append(eqws[0])
err_epoch.append(errors[0])
eqw_epoch.append(eqws[10])
err_epoch.append(errors[10])
eqw_epoch.append(eqws[1])
err_epoch.append(errors[1])
eqw_epoch.append(eqws[2])
err_epoch.append(errors[2])
eqw_epoch.append(eqws[3])
err_epoch.append(errors[3])
eqw_epoch.append(eqws[4])
err_epoch.append(errors[4])
eqw_epoch.append(eqws[5])
err_epoch.append(errors[5])
eqw_epoch.append(eqws[6])
err_epoch.append(errors[6])
eqw_epoch.append(eqws[7])
err_epoch.append(errors[7])
eqw_epoch.append(eqws[8])
err_epoch.append(errors[8])
eqw_epoch.append(eqws[9])
err_epoch.append(errors[9])

print err_epoch

	
epochs=[-13.110599,-10.041886,-6.9081770,-2.9512300,0.039233502,3.2484383,9.1496879,9.1496879,20.687449,26.685855,40.447814]
ions=[]
for j in range(len(lines)):
	ion=[]
	for i in range(len(eqw_epoch)):
		eqw=eqw_epoch[i]
		ion.append(eqw[j])
	print ion
	ions.append(ion)
	
errs = []
for j in range(len(lines)):
	err = []
	for i in range(len(err_epoch)):
		e = err_epoch[i]
		err.append(e[j])
	print err
	errs.append(err)

deltaion=[]
for i in range(len(ions)):
	dion=[]
	for j in range(len(ions[i])):
		dion.append(ions[i][j]-ions[i][4])
	deltaion.append(dion)
	
print ions[0]	

plt.figure(1)
plt.subplot(311)
plt.title('FeII')
plt.errorbar(epochs,ions[0][0],yerr=ions[0][1],label=lines[0])
plt.errorbar(epochs,ions[1][0],yerr=ions[1][1],label=lines[1])
plt.errorbar(epochs,ions[2][0],yerr=ions[2][1],label=lines[2])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylabel('EQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplot(312)
plt.title('MnII')
plt.errorbar(epochs,ions[3][0],yerr=ions[3][1],label=lines[3])
plt.errorbar(epochs,ions[4][0],yerr=ions[4][1],label=lines[4])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylim(0,4)
plt.ylabel('EQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplot(313)
plt.title('MgII & MgI')
plt.errorbar(epochs,ions[5][0],yerr=ions[5][1],label=lines[5])
plt.errorbar(epochs,ions[6][0],yerr=ions[6][1],label=lines[6])
plt.errorbar(epochs,ions[7][0],yerr=ions[7][1],label=lines[7])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylim(0,4)
plt.xlabel('Days since B maximum')
plt.ylabel('EQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplots_adjust(hspace=.35,right=.8)
plt.savefig('Data/eqw/EQW by Epoch,ion.png',orientation='landscape')
plt.show()

plt.figure(2)
for i in range(len(ions)):
	plt.plot(epochs,ions[i][0],ions[i][1],label=lines[i])
legend=plt.legend(loc='upper right')
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.savefig('Data/eqw/EQW by Epoch,all.png',orientation='landscape')
plt.show()


plt.figure(3)
plt.subplot(311)
plt.title('FeII')
plot1,=plt.plot(epochs,deltaion[0],yerr=ions[0][1],label=lines[0])
plot2,=plt.plot(epochs,deltaion[1],yerr=ions[1][1],label=lines[1])
plot3,=plt.plot(epochs,deltaion[2],yerr=ions[2][1],label=lines[2])
#plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylabel('deltaEQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplot(312)
plt.title('MnII')
plot1,=plt.plot(epochs,deltaion[3],yerr=ions[3][1],label=lines[3])
plot2,=plt.plot(epochs,deltaion[4],yerr=ions[4][1],label=lines[4])
#plt.xlim(epochs[0],epochs[len(epochs)-1])
#plt.ylim(0,4)
plt.ylabel('deltaEQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplot(313)
plt.title('MgII & MgI')
plot1,=plt.plot(epochs,deltaion[5],yerr=ions[5][1],label=lines[5])
plot2,=plt.plot(epochs,deltaion[6],yerr=ions[6][1],label=lines[6])
plot3,=plt.plot(epochs,deltaion[7],yerr=ions[7][1],label=lines[7])
#plt.xlim(epochs[0],epochs[len(epochs)-1])
#plt.ylim(0,4)
plt.xlabel('Days since B maximum')
plt.ylabel('deltaEQW')
legend=plt.legend(bbox_to_anchor=(1.05,1), loc=2,shadow=True)
plt.subplots_adjust(hspace=.35,right=.8)
plt.savefig('Data/eqw/deltaEQW by Epoch,ion.png',orientation='landscape')
plt.show()

plt.figure(4)
for i in range(len(ions)):
	plt.plot(epochs,deltaion[i],yerr=ions[i][1],label=lines[i])
legend=plt.legend(loc='upper right')
#plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.savefig('Data/eqw/deltaEQW by Epoch,all.png',orientation='landscape')
plt.show()


