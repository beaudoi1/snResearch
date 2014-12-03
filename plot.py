import matplotlib.pyplot as plt
import numpy as np
import operator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

colors = {'92a': 'black', '11by': 'blue', '11ek': 'brown', '11fe': 'green', '11iv': 'magenta', '12cg': 'orange', '13dy': 'purple', '14j': 'black'}
markers = {'92a': '+', '11by': '^', '11ek': 'p', '11fe': 'o', '11iv': '*', '12cg': '>', '13dy': 'D', '14j': 's'}
sn = ['92a', '11by', '11ek', '11fe', '11iv', '12cg', '13dy', '14j']
font = {'family':'serif' , 'color' :'black' , 'size':10,} 

def plot(ion, iname):
	fig = plt.figure(dpi = 100, figsize = [12, 6], facecolor = 'w')
	ax = fig.add_subplot(111)
	
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
			deltas = ew - np.average(ew, weights = weights)
		except:
			print "No weights"
		
		depoch = []
		deltaews = []
		dlo = []
		dhi = []
		
		j = 0
		for i in range(len(ews)):
			if((lo[i] == 0) & (hi[i] == 0)):
				continue
			else:
				depoch.append(epochs[i])
				deltaews.append(deltas[j])
				dlo.append(lo[i])
				dhi.append(hi[i])
				j += 1

			
		print s	
		print len(depoch), len(deltaews), len(dlo), len(dhi)
		if (len(dlo) != 0):
			label = s
			if (s == '92a'):
				label = '92A'
			if (s == '14j'):
				label = '14J'
			plt.scatter(depoch, deltaews, color = colors[s], marker = markers[s], label = label)
			plt.errorbar(depoch, deltaews, yerr = [dlo, dhi], color = colors[s], linewidth = 2)

	plt.ylabel('$\Delta$EQW ($\AA$)', fontdict = font)
	plt.xlabel('Rest-frame Days Relative to Maximum Brightness', fontdict = font)
	plt.hlines(0, -20, 50, linestyles = 'dashed')
	plt.xlim(-15, 45)
	plt.ylim(-1, 6)
	legend = plt.legend(bbox_to_anchor = (.85, 1), loc = 2, shadow = False, numpoints = 1, prop = {'family': 'serif'})
	majorLocator = MultipleLocator(5)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(1)
	ax.xaxis.set_major_locator(majorLocator)
	ax.xaxis.set_major_formatter(majorFormatter)
	ax.xaxis.set_minor_locator(minorLocator)
	majorLocator = MultipleLocator(1)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(.25)
	ax.yaxis.set_major_locator(majorLocator)
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.yaxis.set_minor_locator(minorLocator)
	plt.savefig(iname + '_delta.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
	fig = plt.figure(dpi = 100, figsize = [12, 6], facecolor = 'w')
	ax = fig.add_subplot(111)
	
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
		
		label = s
		if (s == '92a'):
			label = '92A'
		if (s == '14j'):
			label = '14J'
		plt.scatter(epochs, ews, color = colors[s], marker = markers[s], label = s)
		if (lo[0] != 0):
			plt.errorbar(epochs, ews, yerr = [lo, hi], color = colors[s], linewidth = 2)
		else:
			plt.errorbar(epochs, ews, yerr = [lo, hi], color = colors[s], linestyle = 'None', linewidth = 2)
		
		for i in range(len(lo)):
			if ((lo[i] == 0) & (hi[i] == 0)):
				plt.arrow(epochs[i], ews[i], 0, -.4, head_width = .3, head_length = .15, fc = colors[s], ec = colors[s])
				
	plt.ylabel('EQW ($\AA$)', fontdict = font)
	plt.xlabel('Rest-frame Days Relative to Maximum Brightness', fontdict = font)
	plt.hlines(0, -20, 50, linestyles = 'dashed')
	plt.xlim(-15, 45)
	plt.ylim(-1, 6)
	legend = plt.legend(bbox_to_anchor = (.85, 1), loc = 2, shadow = False, numpoints = 1, prop = {'family': 'serif'})
	majorLocator = MultipleLocator(5)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(1)
	ax.xaxis.set_major_locator(majorLocator)
	ax.xaxis.set_major_formatter(majorFormatter)
	ax.xaxis.set_minor_locator(minorLocator)
	majorLocator = MultipleLocator(1)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(.25)
	ax.yaxis.set_major_locator(majorLocator)
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.yaxis.set_minor_locator(minorLocator)
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

	