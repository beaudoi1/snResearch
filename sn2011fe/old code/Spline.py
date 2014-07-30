import numpy
import scipy
from numpy import linspace
import scipy.interpolate as inter
from astropy.table import Table

data=Table.read('sn2011fe-visit5-uv.flm',format='ascii')
x=data["col1"]
y=data["col2"]
nolines=np.where((x<2500) | ((x>2600) & (x<2700)) | (x>2800))
sm1=inter.InterpolatedUnivariateSpline(x[nolines],y[nolines])
xx=linspace(min(x),max(x),100)
y1=sm1(xx)
sm2=inter.InterpolatedUnivariateSpline(xx,y1)
y2=sm2(x)
norm = y/y2

line1=np.where(x>..,x<)
EW=1-sum(norm[line1]*deltax)