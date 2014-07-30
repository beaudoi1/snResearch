import numpy
from astropy.table import Table
import matplotlib.pyplot as plt
import glob

def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=numpy.ones(window_len,'d')
        else:  
                w=eval('numpy.'+window+'(window_len)')
        y=numpy.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]
	
files=glob.glob("Data/data/*.flm")
	
sp_names=[]
for i in range(len(files)):
	name=''
	file = files[i]
	name = file[19:]
	name = name[:-4]
	sp_names.append(name)	
	
waves=[]
errs=[]
for i in range(len(files)):
	data=Table.read(files[i],format='ascii')
	w=data["col1"]
	f=data["col2"]

	sf=smooth(f)
	res=f-sf
	ares=numpy.abs(res)
	se = smooth(ares)
	#res_err=res/se
	wadd=[max(w)+1,max(w)+2]
	erradd=[0,0]
	new_w=numpy.concatenate((w,wadd))
	new_err=numpy.concatenate((se,erradd))
	waves.append(new_w)
	errs.append(new_err)	

for i in range(len(files)):
	t=Table([waves[i],errs[i]],names=("col1","col2"))
	t.write('Data/data/errors/'+sp_names[i]+'.dat',format='ascii')

"""
plt.figure(1)
plt.plot(w,f)
plt.plot(w,sf)
plt.show()

plt.figure(2)
plt.plot(w,res)
plt.plot(w,se)
plt.show()
"""


