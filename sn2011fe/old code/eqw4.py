import pyspeckit
from astropy.table import Table

#deredshift data

# Rest wavelengths of the lines we are fitting - use as initial guesses
FeIIa = 2346.17
FeIIb = 2376.62
FeIIc = 2384.9
MnIIa = 2588.77
MnIIb = 2601.14
MgIIa = 2798.7
MgIIb = 2805.98
MgIa = 2854.97

#Reads in residual for the spectrum
res=Table.read('data/Visit 4/residual.flm', format='ascii')

#Creates a spectrum from the residual
sp=pyspeckit.Spectrum(xarr=res["col1"],data=res["col2"])

#parts of the continuum to exclude (absorption lines)
exclude=[2334.97,2355.38,2369.27,2379.75,2380.51,2390.63,2581.73,2592,2597,2606.61,2792.1,2801.51,2802.49,2811.56,2848.32,2860.29]

#height, rest wavelength, width of each absorption line
guesses=[-.392,FeIIa,1.46,-.243,FeIIb,1.38,-.504,FeIIc,1.51,-.372,MnIIa,1.7,-.460,MnIIb,1.71,-.685,MgIIa,1.72,-.630,MgIIb,1.6,-.155, MgIa, 1.62]

#fits the baseline to the residual
sp.baseline(xmin=2331,xmax=2900,exclude=exclude,subtract=False,reset_selection=True)

#Fits a gaussian to each line (height, rest wavelength, width)
sp.specfit(guesses=guesses)

#Gives the EW of each line
ions=['FeIIa','FeIIb','FeIIc','MnIIa','MnIIb','MgIIa','MgIIb','MgIa']	
eq=sp.specfit.EQW(components=True)
t=Table([ions,eq],names=('col1','col2'))
t.write('data/Visit 4/EQW.dat',format='ascii')
