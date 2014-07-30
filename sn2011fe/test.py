import numpy as np
from numpy import linspace
import scipy.interpolate as inter
from astropy.table import Table
import matplotlib.pyplot as plt
import glob
import sys
import pyspeckit
import math

def rand_break(breakpoints):
	bkpts=[]
	for b in breakpoints:
		bkpts.append(b)
	print bkpts
	febkpts,mnbkpts,mg2bkpts,mg1bkpts,other=[],[],[],[],[]
	for b in bkpts:
		if ((b>2331)&(b<2397)):
			febkpts.append(b)
		if ((b>2578)&(b<2614)):
			mnbkpts.append(b)
		if ((b>2783)&(b<2820)):
			mg2bkpts.append(b)
		if ((b>2842)&(b<2865)):
			mg1bkpts.append(b)
		else:
			other.append(b)
	print febkpts,mnbkpts,mg2bkpts,mg1bkpts,other
	new_bkpts=[]
	feshift=np.random.randn(len(febkpts)) * 4
	mnshift=np.random.randn(len(mnbkpts)) * 5
	mg2shift=np.random.randn(len(mg2bkpts)) * 3.25
	mg1shift=np.random.randn(len(mg1bkpts)) * 3.5
	new_fe=feshift+febkpts
	new_mn=mnshift+mnbkpts
	new_mg2=mg2shift+mg2bkpts
	new_mg1=mg1shift+mg1bkpts
	print new_fe, new_mn, new_mg2,new_mg1
	for fe in new_fe:
		new_bkpts.append(fe)
	for mn in new_mn:
		new_bkpts.append(mn)
	for mg2 in new_mg2:
		new_bkpts.append(mg2)
	for mg1 in new_mg1:
		new_bkpts.append(mg1)
	for o in other:
		new_bkpts.append(o)
	print new_bkpts
	new_bkpts.sort()
	print new_bkpts
	i=input("ss")
		
bkpts_files=glob.glob('Data/breakpoints/test_bkpts/*flm')
num=len(bkpts_files)
breakpoints=[]
for i in range(len(bkpts_files)):
	bkpts=Table.read(bkpts_files[i],format='ascii')
	breakpoints.append(bkpts["col1"])
	
for a in range(num):
	mc=10
	for i in range(mc):
		bps=rand_break(breakpoints[a])
		print bps