import numpy as np
import glob
from astropy.table import Table

bkpts_files=glob.glob('Data/breakpoints/*.flm')
breakpoints=[]
for i in range(len(bkpts_files)):
	bkpts=Table.read(bkpts_files[i],format='ascii')
	breakpoints.append(bkpts["col1"])
	
print breakpoints[0]

bshift =  np.random.randn(len(breakpoints[0])) * 2

print bshift