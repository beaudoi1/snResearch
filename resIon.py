import matplotlib.pyplot as plt
import numpy as np
import operator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math

colors = {'92a': 'black', '11by': 'blue', '11ek': 'brown', '11fe': 'green', '11iv': 'magenta', '12cg': 'orange', '13dy': 'purple', '14j': 'black'}
markers = {'92a': '+', '11by': '<', '11ek': 'p', '11fe': 'o', '11iv': '*', '12cg': '>', '13dy': 'D', '14j': 's'}
sn = ['92a', '11by', '11ek', '11fe', '11iv', '12cg', '13dy', '14j']
font = {'family':'serif' , 'color' :'black' , 'size':24,} 
rwaves = {'FeIIa': 'Fe II $\lambda$2344.21', 'FeIIb': 'Fe II $\lambda$2374.4601', 'FeIIc': 'Fe II $\lambda$2382.76', 'FeIId': 'Fe II $\lambda$2586.65', 'FeIIe': 'Fe II $\lambda$2600.17', 'MgIIa': 'Mg II $\lambda$2796.35', 'MgIIb': 'Mg II $\lambda$2803.53', 'MgI': 'Mg I $\lambda$2852.96'}

def plotDelta(ions, inames):
	fig = plt.figure(dpi = 100, figsize = [12, 8 * len(ions)], facecolor = 'w')
	count = 1
	for ion in ions:
		sub = str(len(ions)) + '1' + str(count)
		sub = int(sub)
		ax = fig.add_subplot(sub)
		
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
				plt.scatter(depoch, deltaews, color = colors[s], marker = markers[s], s = 50, label = label)
				plt.errorbar(depoch, deltaews, yerr = [dlo, dhi], color = colors[s], linewidth = 2)
				
				for i in range(len(dlo)):
					if ((dlo[i] == 0) & (dhi[i] == 0)):
						plt.arrow(depoch[i], deltaews[i], 0, -.6, head_width = .5, head_length = .25, fc = colors[s], ec = colors[s])
		
		if (count == math.ceil(float(len(ions)) / 2)):
			plt.ylabel('$\Delta$EW ($\AA$)', fontdict = font)
		plt.hlines(0, -20, 50, linestyles = 'dashed')
		if (count != len(ions)):
			for x in plt.gca().axes.get_xticklabels():
				x.set_visible(False)
				x.set_fontsize(0.0)
		plt.xlim(-20, 47)
		plt.ylim(math.floor(min * 10) / 10, math.ceil(max * 10) / 10)
		if (count == 1):	
			legend = plt.legend(loc = 'upper right', ncol = 2, columnspacing = .1, handletextpad = .1, borderpad = .2, shadow = False, scatterpoints = 1, prop = {'family': 'serif', 'size':18,})
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
		for tick in ax.yaxis.get_major_ticks():
			tick.label.set_fontsize(16)
		name = rwaves[inames[count - 1]]
		plt.text(-18.5, (math.ceil(max * 10) / 10) - 1, name, fontsize = 18, color = 'k', family = 'serif')
	
		count += 1
	plt.gca().axes.get_xaxis().set_visible(True)
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	plt.xlabel('Rest-frame Days Relative to Maximum Brightness', fontdict = font)
	plt.subplots_adjust(hspace = 0, right = 0.97, left = 0.08, top = 0.97, bottom = 0.1)
	plt.savefig('delta_res_' + inames[0] + '.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	
def plotEw(ions, inames):
	fig = plt.figure(dpi = 100, figsize = [12, 8 * len(ions)], facecolor = 'w')
	count = 1
	for ion in ions:
		sub = str(len(ions)) + '1' + str(count)
		sub = int(sub)
		ax = fig.add_subplot(sub)
		
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
			"""

			label = s
			if (s == '92a'):
				label = '92A'
			if (s == '14j'):
				label = '14J'
			
			"""
			plot = False
			for l in lo:
				if l != 0:
					plot = True
			"""		
			plot = True
				
			if plot:
				for i in range(len(hi)):
					if (ews[i] < 6):
						h = hi[i] + ews[i] + .5
						l = ews[i] - lo[i] - .9
						if (h > max):
							max = h
						if (l < min):
							min = l
				plt.scatter(epochs, ews, color = colors[s], marker = markers[s], s = 50, label = label)
				plt.errorbar(epochs, ews, yerr = [lo, hi], color = colors[s], linewidth = 2)
				
				for i in range(len(lo)):
					if ((lo[i] == 0) & (hi[i] == 0)):
						plt.arrow(epochs[i], ews[i], 0, -.6, head_width = .5, head_length = .25, fc = colors[s], ec = colors[s])
		
		if (count == math.ceil(float(len(ions)) / 2)):
			plt.ylabel('EW ($\AA$)', fontdict = font)
		#plt.hlines(0, -20, 50, linestyles = 'dashed')
		if (count != len(ions)):
			for x in plt.gca().axes.get_xticklabels():
				x.set_visible(False)
				x.set_fontsize(0.0)
		plt.xlim(-20, 47)
		plt.ylim(math.floor(min * 10) / 10, math.ceil(max * 10) / 10)
		if (count == 1):	
			legend = plt.legend(loc = 'upper right', ncol = 2, columnspacing = .1, handletextpad = .1, borderpad = .2, shadow = False, scatterpoints = 1, prop = {'family': 'serif', 'size':18,})
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
		for tick in ax.yaxis.get_major_ticks():
			tick.label.set_fontsize(16)
		name = rwaves[inames[count - 1]]
		plt.text(-18.5, (math.ceil(max * 10) / 10) - 1, name, fontsize = 18, color = 'k', family = 'serif')
	
		count += 1
	plt.gca().axes.get_xaxis().set_visible(True)
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	plt.xlabel('Rest-frame Days Relative to Maximum Brightness', fontdict = font)
	plt.subplots_adjust(hspace = 0, right = 0.97, left = 0.08, top = 0.97, bottom = 0.1)
	plt.savefig('res_' + inames[0] + '.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	plt.show()
	