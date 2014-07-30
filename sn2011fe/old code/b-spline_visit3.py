import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from astropy.table import Table
import matplotlib.pyplot as plt


#read in file and seperate wavelength(x) and flux(y)
data=Table.read('data/Visit 3/sn2011fe-visit3-uv.flm',format='ascii')
x=data["col1"]
y=data["col2"]

#remove the absorption lines from the spectrum
nolines=np.where((x<2341) | ((x>2351) & (x<2373)) | ((x>2380) & (x<2381)) | ((x>2390) & (x<2585)) | ((x>2593) & (x<2597)) | ((x>2607) & (x<2791)) | ((x>2812) & (x<2848)) | (x>2860))


#creates b-spline over spectrum
sm1=inter.splrep(x[nolines],y[nolines])
xx=linspace(min(x),max(x),65)

#breakpoints in spectrum
aa=np.array([2360,2369,2390,2788,2848,2863])
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
plt.savefig('data/Visit 3/b-spline.png')
plt.show()

#plots residual
plot1,=plt.plot(x,norm)
plt.ylim(0,2)
plt.savefig('data/Visit 3/residual.png')
plt.show()

#records residual plot
t=Table([x,norm],names=('col1','col2'))
t.write('data/Visit 3/residual.flm',format='ascii')