import pyspeckit
from astropy.table import Table
import glob
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#deredshift data

# Rest wavelengths of the lines we are fitting - use as initial guesses
FeII_emit = [2344.2128,2374.4601,2382.7641]
MnII_emit = [2594.4967,2606.4588]
MgII_emit = [2796.3540,2803.5311]
MgI_emit = [2852.9628]
FeII_obs = [2346.17,2376.62,2384.9]
MnII_obs = [2588.77,2601.14]
MgII_obs = [2798.7,2805.98]
MgI_obs = [2854.97]

files = glob.glob("Data/residuals/*.flm")

#Reads in residual for each spectrum
res_names=[]
waves=[]
fluxes=[]
for i in range(len(files)):
	res=Table.read(files[i],format='ascii')
	name=files[i][15:]
	name=name[:-4]
	res_names.append(name)
	waves.append(res["col1"])
	fluxes.append(res["col2"])
	
def find_eqw(x,y,guesses,spectra_name):
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
	sp.plotter.savefig('Data/Pictures/Lines/' +spectra_name+'.png')
	

	#Gives the EW of each line
	eq=sp.specfit.EQW(components=True,annotate=True)
	
	return eq

eqws=[]
for i in range(len(res_names)):
	guesses=[-.41,FeII_obs[0],1.65,-.246,FeII_obs[1],1.31,-.518,FeII_obs[2],1.58,-.38,MnII_obs[0],1.6,-.484,MnII_obs[1],2.07,-.696,MgII_obs[0],1.699,-.638,MgII_obs[1],1.6,-.161,MgI_obs[0],1.62]
	print "Eqw of " + res_names[i]
	eqws.append(find_eqw(waves[i],fluxes[i],guesses,res_names[i]))
	print eqws[i]

ions=['FeIIa','FeIIb','FeIIc','MnIIa','MnIIb','MgIIa','MgIIb','MgIa']	
for i in range(len(res_names)):
	t=Table([ions,eqws[i]],names=('col1','col2'))
	print res_names[i]
	print t
	t.write('Data/eqw/' + res_names[i] + '.dat',format='ascii')

"""

def find_redshift():
	lines=input('Wavelength of lines: ')
	waves=[]
	FeRed=[]
	MnRed=[]
	Mg2Red=[]
	Mg1Red=[]
	for i in range(len(lines)):
		waves.append(float(lines[i]))
	for k in range(len(waves)):
		for j in range(len(FeII)):
			diff=math.fabs(waves[k]-FeII[j])
			if (diff<3.5):
				red=diff/FeII[j]
				FeRed.append(red)
		for a in range(len(MnII)):
			diff=math.fabs(waves[k]-MnII[a])
			if (diff<3.5):
				red=diff/MnII[a]
				MnRed.append(red)
		for b in range(len(MgII)):
			diff=math.fabs(waves[k]-MgII[b])
			if (diff<3.5):
				red=diff/MgII[b]
				Mg2Red.append(red)
		for c in range(len(MgI)):
			diff=math.fabs(waves[k]-MgI[c])
			if (diff<3.5):
				red=diff/MgI[c]
				Mg1Red.append(red)
	
	Fe=np.average(FeRed)
	Mn=np.average(MnRed)
	Mg2=np.average(Mg2Red)
	Mg1=np.average(Mg1Red)
	
	return Fe,Mn,Mg2,Mg1
	
def find_guesses(Fe,Mn,Mg2,Mg1):
	widths=[]
	heights=[]
	waves=[]
	fe_lines=raw_input("Are there FeII lines?: (y/n) ")
	if fe_lines=='y':
		for i in range(len(FeII)):
			wave=FeII[i]*(Fe+1)
			waves.append(wave)
		g_height=input("Input three heights: ")
		for j in range(len(g_height)):
			heights.append(float(g_height[j]))
		g_width=input("Input three widths: ")
		for k in range(len(g_width)):
			widths.append(float(g_width[k]))
	mn_lines=raw_input("Are there MnII lines?: (y/n) ")
	if mn_lines=='y':
		for i in range(len(MnII)):
			wave=MnII[i]*(Mn+1)
			waves.append(wave)
		g_height=input("Input two heights: ")
		for j in range(len(g_height)):
			heights.append(float(g_height[j]))
		g_width=input("Input two widths: ")
		for k in range(len(g_width)):
			widths.append(float(g_width[k]))
	mg2_lines=raw_input("Are the MgII lines?: (y/n) ")
	if mg2_lines=='y':
		for i in range(len(MgII)):
			wave=MgII[i]*(Mg2+1)
			waves.append(wave)
		g_height=input("Input two heights: ")
		for j in range(len(g_height)):
			heights.append(float(g_height[j]))
		g_width=input("Input two widths: ")
		for k in range(len(g_width)):
			widths.append(float(g_width[k]))
	mg1_lines=raw_input("Are the MgI lines?: (y/n) ")
	if mg1_lines=='y':
		for i in range(len(MgI)):
			wave=MgI[i]*(Mg1+1)
			waves.append(wave)
		g_height=input("Input two heights: ")
		for j in range(len(g_height)):
			heights.append(float(g_height[j]))
		g_width=input("Input two widths: ")
		for k in range(len(g_width)):
			widths.append(float(g_width[k]))
	guesses=[]
	for a in range(len(widths)):
		guesses.append(heights[a])
		guesses.append(waves[a])
		guesses.append(widths[a])
	return guesses
	
for i in range(len(res_names)):	
	print "Find wavelengths"
	plt.figure(i)
	plt.title(res_names[i])
	plt.plot(waves[i],fluxes[i])
	plt.show()
	fe,mn,mg2,mg1=find_redshift()
	Fe.append(fe)
	Mn.append(mn)
	Mg2.append(mg2)
	Mg1.append(mg1)
	guesses=find_guesses(fe,mn,mg2,mg1)
	eqws.append(find_eqw(waves[i],fluxes[i],guesses))
"""
