from astropy.table import Table, Column
import sys
import glob

files = glob.glob("Data/data/*.flm")
num= len(files)
sp_names=[]
for i in range(num):
	name=''
	file = files[i]
	name = file[19:]
	name = name[:-4]
	sp_names.append(name)


breakpoints=[]
for i in range(num):
	points=input('Breakpoints for ' + sp_names[i] + ': ')
	bps=[]
	for j in range(len(points)):
		bps.append(float(points[j]))
	breakpoints.append(bps)

zeros=[]
for i in range(len(breakpoints)):
	zero=[]
	for j in range(len(breakpoints[i])):
		zero.append(0)
	zeros.append(zero)
	
for i in range(num):		
	t=Table([breakpoints[i],zeros[i]],names=('col1','col2'))
	t.write('Data/breakpoints/'+sp_names[i]+'_breakpoints.flm',format='ascii')
	
	