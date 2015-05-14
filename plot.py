import matplotlib.pyplot as plt
import numpy as np
import operator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math

colors = {'92a': 'black', '11by': 'blue', '11ek': 'brown', '11fe': 'green', '11iv': 'magenta', '12cg': 'orange', '13dy': 'purple', '14j': 'black'}
markers = {'92a': '+', '11by': '^', '11ek': 'p', '11fe': 'o', '11iv': '*', '12cg': '>', '13dy': 'D', '14j': 's'}
sn = ['92a', '11by', '11ek', '11fe', '11iv', '12cg', '13dy', '14j']
font = {'family':'serif' , 'color' :'black' , 'size':24,} 

def plot(ion, iname):
	fig = plt.figure(dpi = 100, figsize = [10, 8], facecolor = 'w')
	ax1 = fig.add_subplot(111)
	
	max = -999
	min = 999
	for s in sn:
		set = ion[s]
		sorted_set = sorted(set.items(), key = operator.itemgetter(0))
		set = sorted_set
		epochs = []
		ews = []
		hi = []
		lo = []
		for i in set:
			epochs.append(i[0])
			ews.append(i[1][0])
			lo.append(i[1][1][0])
			hi.append(i[1][1][1])
		"""
		epochs = set.keys()
		ews = []
		hi = []
		lo = []
		for e in epochs:
			ews.append(set[e][0])
			lo.append(set[e][1][0])
			hi.append(set[e][1][1])
		"""
			
		weights = []
		ew = []
		for i in range(len(hi)):
			if ((lo[i] != 0) | (hi[i] != 0)):
				avg = (lo[i] + hi[i]) / 2
				var = (1/avg**2)
				weights.append(var)
				ew.append(ews[i])
	
		try:	
			average = np.average(ew, weights = weights)	
			deltas = ew - average
		except:
			print "No weights"
		
		depoch = epochs
		deltaews = []
		dlo = lo
		dhi = hi
		
		j = 0
		for i in range(len(ews)):	
			if((lo[i] == 0) & (hi[i] == 0)):
				deltaews.append(ews[i])
			else:	
				deltaews.append(deltas[j])
				j += 1

		label = s
		if (s == '92a'):
			label = '92A'
		if (s == '14j'):
			label = '14J'
		plot = False
		for l in dlo:
			if l != 0:
				plot = True
		if plot:
			for i in range(len(dhi)):
				if (deltaews[i] < 6):
					h = dhi[i] + deltaews[i] + .5
					l = deltaews[i] - dlo[i] - .45
					if (h > max):
						max = h
					if (l < min):
						min = l
			plt.scatter(depoch, deltaews, color = colors[s], marker = markers[s], label = label)
			plt.errorbar(depoch, deltaews, yerr = [dlo, dhi], color = colors[s], linewidth = 2)
			
			for i in range(len(dlo)):
				if ((dlo[i] == 0) & (dhi[i] == 0)):
					plt.arrow(depoch[i], deltaews[i], 0, -.3, head_width = .3, head_length = .1, fc = colors[s], ec = colors[s])
			

	plt.ylabel('$\Delta$EW ($\AA$)', fontdict = font)
	plt.xlabel('Rest-frame Days Relative to Maximum Brightness', fontdict = font)
	plt.hlines(0, -20, 50, linestyles = 'dashed')
	plt.xlim(-20, 45)
	plt.ylim(math.floor(min * 10) / 10, math.ceil(max * 10) / 10)
	legend = plt.legend(loc = 'upper right', shadow = False, numpoints = 1, prop = {'family': 'serif', 'size':24,})
	majorLocator = MultipleLocator(5)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(1)
	ax1.xaxis.set_major_locator(majorLocator)
	ax1.xaxis.set_major_formatter(majorFormatter)
	ax1.xaxis.set_minor_locator(minorLocator)
	majorLocator = MultipleLocator(1)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(.25)
	ax1.yaxis.set_major_locator(majorLocator)
	ax1.yaxis.set_major_formatter(majorFormatter)
	ax1.yaxis.set_minor_locator(minorLocator)
	name = iname[:2] + " " + iname[2:]
	plt.text(-15, (math.ceil(max * 10) / 10) - .5, name, fontsize = 18, color = 'k', family = 'serif')
	plt.savefig(iname + '_delta.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
	fig = plt.figure(dpi = 100, figsize = [10, 8], facecolor = 'w')
	ax2 = fig.add_subplot(111)
	max = -999
	min = 999
	for s in sn:
		set = ion[s]
		sorted_set = sorted(set.items(), key = operator.itemgetter(0))
		set = sorted_set
		epochs = []
		ews = []
		hi = []
		lo = []
		for i in set:
			epochs.append(i[0])
			ews.append(i[1][0])
			lo.append(i[1][1][0])
			hi.append(i[1][1][1])
		
		for i in range(len(hi)):
			if (ews[i] < 6):
				h = hi[i] + ews[i] + .5
				l = ews[i] - lo[i] - .45
				if (h > max):
					max = h
				if (l < min):
					min = l	
		
		label = s
		if (s == '92a'):
			label = '92A'
		if (s == '14j'):
			label = '14J'
		plt.scatter(epochs, ews, color = colors[s], marker = markers[s], label = label)
		if (lo[0] != 0):
			plt.errorbar(epochs, ews, yerr = [lo, hi], color = colors[s], linewidth = 2)
		else:
			plt.errorbar(epochs, ews, yerr = [lo, hi], color = colors[s], linestyle = 'None', linewidth = 2)
		
		for i in range(len(lo)):
			if ((lo[i] == 0) & (hi[i] == 0)):
				plt.arrow(epochs[i], ews[i], 0, -.3, head_width = .3, head_length = .1, fc = colors[s], ec = colors[s])
				
	plt.ylabel('EW ($\AA$)', fontdict = font)
	plt.xlabel('Rest-frame Days Relative to Maximum Brightness', fontdict = font)
	plt.hlines(0, -20, 50, linestyles = 'dashed')
	plt.xlim(-20, 45)
	plt.ylim(math.floor(min * 10) / 10, math.ceil(max * 10) / 10)
	legend = plt.legend(loc = 'upper right', shadow = False, numpoints = 1, prop = {'family': 'serif', 'size':24,})
	majorLocator = MultipleLocator(5)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(1)
	ax2.xaxis.set_major_locator(majorLocator)
	ax2.xaxis.set_major_formatter(majorFormatter)
	ax2.xaxis.set_minor_locator(minorLocator)
	majorLocator = MultipleLocator(1)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(.25)
	ax2.yaxis.set_major_locator(majorLocator)
	ax2.yaxis.set_major_formatter(majorFormatter)
	ax2.yaxis.set_minor_locator(minorLocator)
	plt.savefig(iname + '.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
def plotMW(ion):
	set = ion
	sorted_set = sorted(set.items(), key = operator.itemgetter(0))
	set = sorted_set
	epochs = []
	ews = []
	hi = []
	lo = []
	for i in set:
		epochs.append(i[0])
		ews.append(i[1][0])
		lo.append(i[1][1][0])
		hi.append(i[1][1][1])
		
	weights = []
	ew = []
	for i in range(len(hi)):
		if ((lo[i] != 0) | (hi[i] != 0)):
			avg = (lo[i] + hi[i]) / 2
			var = (1/avg**2)
			weights.append(var)
			ew.append(ews[i])
			
	try:		
		deltas = ew - np.average(ew, weights = weights)
	except:
		print "No weights"
	
	deltaews = []
	
	j = 0
	for i in range(len(ews)):
		if((lo[i] == 0) & (hi[i] == 0)):
			deltaews.append(ews[i])
		else:
			deltaews.append(deltas[j])
			j += 1
			
	plt.scatter(epochs, deltaews, color = 'black', marker = '<', label = "MW")
	plt.errorbar(epochs, deltaews, yerr = [lo, hi], color = 'black')
	legend = plt.legend(bbox_to_anchor = (.85, 1), loc = 2, shadow = False, numpoints = 1, prop = {'family': 'serif'})

	