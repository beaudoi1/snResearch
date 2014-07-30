import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from astropy.table import Table
import matplotlib.pyplot as plt


#read in file and seperate wavelength(x) and flux(y)
data=Table.read('data/Visit 6.2/sn2011fe-visit6-fuv.flm',format='ascii')
x=data["col1"]
y=data["col2"]

#remove the absorption lines from the spectrum
nolines=np.where((x<2335) | ((x>2352) & (x<2367.5)) | ((x>2389) & (x<2381)) | ((x>2391) & (x<2580)) | ((x>2591) & (x<2594)) | ((x>2605) & (x<2788)) | ((x>2809) & (x<2847)) | (x>2858))


#creates b-spline over spectrum
sm1=inter.splrep(x[nolines],y[nolines])
xx=linspace(min(x),max(x),65)

#breakpoints in spectrum
aa=np.array([2334,2363,2422,2490,2512,2579,2609,2843,2882])
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
plt.savefig('data/Visit 6.2/b-spline.png')
plt.show()

#plots residual
plot1,=plt.plot(x,norm)
plt.ylim(0,2)
plt.savefig('data/Visit 6.2/residual.png')
plt.show()

#records residual plot
t=Table([x,norm],names=('col1','col2'))
t.write('data/Visit 6.2/residual.flm',format='ascii')