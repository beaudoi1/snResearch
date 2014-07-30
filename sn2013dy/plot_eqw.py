import matplotlib.pyplot as plt
import glob
from astropy.table import Table
import numpy as np

files = glob.glob("Data/eqw/*.dat")
err_files = glob.glob("Data/eqw/MonteCarlo/error_files/*.dat")

names = []
lines = []
eqws = []
error = []

for i in range(len(files)):
	file = files[i]
	name = file[:-4]
	name = name[9:]
	names.append(name)
	eqw = Table.read(files[i], format='ascii')
	lines = eqw["col1"]
	eqws.append(eqw["col2"])
	error.append(eqw["col3"])

errs_lo = []
errs_hi = []
for i in range(len(err_files)):
	err = Table.read(err_files[i],format='ascii')
	errs_lo.append(err["col2"])
	errs_hi.append(err["col3"])

temp_lo_errs = []
temp_hi_errs = []
for i in range(len(errs_lo)):
	temp_lo_errs.append((error[i]**2 + errs_lo[i]**2)**.5)
	temp_hi_errs.append((error[i]**2 + errs_hi[i]**2)**.5)
	
eqw_order = []
lo_err_order = []
hi_err_order = []
n_order = []

eqw_order.append(eqws[0])
lo_err_order.append(temp_lo_errs[0])
hi_err_order.append(temp_hi_errs[0])
eqw_order.append(eqws[2])
lo_err_order.append(temp_lo_errs[2])
hi_err_order.append(temp_hi_errs[2])
eqw_order.append(eqws[3])
lo_err_order.append(temp_lo_errs[3])
hi_err_order.append(temp_hi_errs[3])
eqw_order.append(eqws[4])
lo_err_order.append(temp_lo_errs[4])
hi_err_order.append(temp_hi_errs[4])
eqw_order.append(eqws[5])
lo_err_order.append(temp_lo_errs[5])
hi_err_order.append(temp_hi_errs[5])
eqw_order.append(eqws[6])
lo_err_order.append(temp_lo_errs[6])
hi_err_order.append(temp_hi_errs[6])
eqw_order.append(eqws[7])
lo_err_order.append(temp_lo_errs[7])
hi_err_order.append(temp_hi_errs[7])
eqw_order.append(eqws[8])
lo_err_order.append(temp_lo_errs[8])
hi_err_order.append(temp_hi_errs[8])
eqw_order.append(eqws[9])
lo_err_order.append(temp_lo_errs[9])
hi_err_order.append(temp_hi_errs[9])
eqw_order.append(eqws[10])
lo_err_order.append(temp_lo_errs[10])
hi_err_order.append(temp_hi_errs[10])
eqw_order.append(eqws[1])
lo_err_order.append(temp_lo_errs[1])
hi_err_order.append(temp_hi_errs[1])
n_order.append(names[0])
n_order.append(names[2])
n_order.append(names[3])
n_order.append(names[4])
n_order.append(names[5])
n_order.append(names[6])
n_order.append(names[7])
n_order.append(names[8])
n_order.append(names[9])
n_order.append(names[10])
n_order.append(names[1])

epochs=[-13.110599,-10.041886,-6.9081770,-2.9512300,0.039233502,3.2484383,9.1496879,9.1496879,20.687449,26.685855,40.447814]

ions = []
errs = []
lo_errs = []
hi_errs = []
for j in range(len(lines)):
	ion = []
	err = []
	lo_err = []
	hi_err = []
	for k in range(len(eqw_order)):
		eq = eqw_order[k]
		lo_e = lo_err_order[k]
		hi_e = hi_err_order[k]
		ion.append(eq[j])
		lo_err.append(lo_e[j])
		hi_err.append(hi_e[j])
	ions.append(ion)
	lo_errs.append(lo_err)
	hi_errs.append(hi_err)
	
deltaion = []
for j in range(len(ions)):
	dion = []
	for k in range(len(ions[j])):
		dion.append(ions[j][k] - ions[j][4])
	deltaion.append(dion)
	
plt.figure(1)
plt.subplot(311)
plt.title('FeII')
plt.errorbar(epochs,ions[0],yerr = [lo_errs[0],hi_errs[0]], label = lines[0])
plt.errorbar(epochs,ions[1],yerr = [lo_errs[1],hi_errs[1]], label = lines[1])
plt.errorbar(epochs,ions[2],yerr = [lo_errs[2],hi_errs[2]], label = lines[2])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylim(0,3)
plt.ylabel('EQW')
legend = plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, shadow = True)
plt.subplot(312)
plt.title('MnII')
plt.errorbar(epochs,ions[3],yerr = [lo_errs[3],hi_errs[3]], label = lines[3])
plt.errorbar(epochs,ions[4],yerr = [lo_errs[4],hi_errs[4]], label = lines[4])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylabel('EQW')
legend = plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, shadow = True)
plt.subplot(313)
plt.title('MgII & MgI')
plt.errorbar(epochs,ions[5],yerr = [lo_errs[5],hi_errs[5]], label = lines[5])
plt.errorbar(epochs,ions[6],yerr = [lo_errs[6],hi_errs[6]], label = lines[6])
plt.errorbar(epochs,ions[7],yerr = [lo_errs[7],hi_errs[7]], label = lines[7])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.xlabel('Days since B maximum')
plt.ylabel('EQW')
legend = plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, shadow = True)
plt.subplots_adjust(hspace = .35, right = .8)
plt.savefig('Data/eqw/EQW by Epoch with Error,ion.png', orientation = 'landscape')
plt.show()

plt.figure(2)
for i in range(len(ions)):
	plt.errorbar(epochs, ions[i], yerr = [lo_errs[i],hi_errs[i]], label = lines[i])
legend = plt.legend(loc = 'upper right')
plt.ylim(0,3)
plt.xlim(epochs[0], epochs[len(epochs) - 1])
plt.savefig('Data/eqw/EQW by Epoch with Error,all.png', orientation = 'landscape')
plt.show()
	
plt.figure(3)
plt.subplot(311)
plt.title('FeII')
plt.errorbar(epochs,deltaion[0],yerr = [lo_errs[0],hi_errs[0]], label = lines[0])
plt.errorbar(epochs,deltaion[1],yerr = [lo_errs[1],hi_errs[1]], label = lines[1])
plt.errorbar(epochs,deltaion[2],yerr = [lo_errs[2],hi_errs[2]], label = lines[2])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylim(-1,1)
plt.ylabel('delta EQW')
legend = plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, shadow = True)
plt.subplot(312)
plt.title('MnII')
plt.errorbar(epochs,deltaion[3],yerr = [lo_errs[3],hi_errs[3]], label = lines[3])
plt.errorbar(epochs,deltaion[4],yerr = [lo_errs[4],hi_errs[4]], label = lines[4])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.ylabel('delta EQW')
legend = plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, shadow = True)
plt.subplot(313)
plt.title('MgII & MgI')
plt.errorbar(epochs,deltaion[5],yerr = [lo_errs[5],hi_errs[5]], label = lines[5])
plt.errorbar(epochs,deltaion[6],yerr = [lo_errs[6],hi_errs[6]], label = lines[6])
plt.errorbar(epochs,deltaion[7],yerr = [lo_errs[7],hi_errs[7]], label = lines[7])
plt.xlim(epochs[0],epochs[len(epochs)-1])
plt.xlabel('Days since B maximum')
plt.ylabel('delta EQW')
legend = plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, shadow = True)
plt.subplots_adjust(hspace = .35, right = .8)
plt.savefig('Data/eqw/deltaEQW by Epoch with Error,ion.png', orientation = 'landscape')
plt.show()

plt.figure(4)
for i in range(len(ions)):
	plt.errorbar(epochs, deltaion[i], yerr = [lo_errs[i],hi_errs[i]], label = lines[i])
legend = plt.legend(loc = 'upper right')
plt.ylim(-1,1)
plt.xlim(epochs[0], epochs[len(epochs) - 1])
plt.savefig('Data/eqw/deltaEQW by Epoch with Error,all.png', orientation = 'landscape')
plt.show()
	
	
	
	
	
	