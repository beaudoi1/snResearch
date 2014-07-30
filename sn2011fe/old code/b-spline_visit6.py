import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from astropy.table import Table
import matplotlib.pyplot as plt


#read in file and seperate wavelength(x) and flux(y)
data=Table.read('Data/data/other_data/sn2011fe-visit6-uv.flm',format='ascii')
x=data["col1"]
y=data["col2"]

#remove the absorption lines from the spectrum
nolines=np.where((x<2334) | ((x>2352) & (x<2368)) | ((x>2391) & (x<2581)) | ((x>2605) & (x<2788)) | ((x>2810) & (x<2846)) | (x>2861))


#creates b-spline over spectrum
sm1=inter.splrep(x[nolines],y[nolines])
xx=linspace(min(x),max(x),65)
nolinesxx=np.where((xx<2334) | ((xx>2352) & (xx<2368)) | ((xx>2391) & (xx<2581)) | ((xx>2605) & (xx<2788)) | ((xx>2810) & (xx<2846)) | (xx>2861))
xx=xx[nolinesxx]

#breakpoints in spectrum
aa=np.array([2351,2391,2580,2812,2845])
bb=np.concatenate((xx,aa)) #add breakpoints to array
bb.sort()

#creates and fits b-spline over wavelength range
y1=inter.splev(bb,sm1)
sm2=inter.splrep(bb,y1)
y2=inter.splev(x,sm2)

#creates residual for spectrum
norm=y/y2

#plots b-spline over spectrum
plt.figure(1)
plt.title("visit6-uv")
#plt.subplot(211)
plt.plot(x,y)
plt.plot(x,y2)
#plt.subplot(212)
#plot1,=plt.plot(x,norm)
#plt.ylim(0,1.5)
#plt.savefig('Data/Pictures/visit6_uv.png')
plt.show()

#records residual plot
#t=Table([x,norm],names=('col1','col2'))
#t.write('Data/residuals/other_residuals/visit6_uv.flm',format='ascii')