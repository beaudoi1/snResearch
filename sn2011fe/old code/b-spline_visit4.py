import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from astropy.table import Table
import matplotlib.pyplot as plt


#read in file and seperate wavelength(x) and flux(y)
data=Table.read('data/Visit 4/sn2011fe-visit4-uv.flm',format='ascii')
x=data["col1"]
y=data["col2"]

#remove the absorption lines from the spectrum
nolines=np.where((x<2340) | ((x>2351) & (x<2371)) | ((x>2392) & (x<2584)) | ((x>2608) & (x<2791)) | ((x>2817) & (x<2849)) | (x>2863))

#creates b-spline over spectrum
sm1=inter.splrep(x[nolines],y[nolines])
xx=linspace(min(x),max(x),65)

#breakpoints in spectrum
aa=np.array([2294,2305,2331,2393,2515,2553,2580,2613,2673,2763,2783,2790,2819,2833,2838,2848,2865,2872.2898,3011,3043,3056])
bb=np.concatenate((xx,aa)) #add breakpoints to array
bb.sort()

#creates and fits b-spline over wavelength range
y1=inter.splev(bb,sm1)
sm2=inter.splrep(bb,y1)
y2=inter.splev(x,sm2)

#creates residual for spectrum
norm=y/y2

#plots b-spline over spectrum
plot1,=plt.plot(x,y)
plot2,=plt.plot(x,y2)
plt.savefig('data/Visit 4/b-spline.png')
plt.show()

#plots residual
plot1,=plt.plot(x,norm)
plt.ylim(0,2)
plt.savefig('data/Visit 4/residual.png')
plt.show()

#records residual plot
t=Table([x,norm],names=('col1','col2'))
t.write('data/Visit 4/residual.flm',format='ascii')


