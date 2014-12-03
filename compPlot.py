import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import operator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

sn = ['92a', '11by', '11ek', '11fe', '11iv', '12cg', '13dy', '14j']
colors = {'92a': 'black', '11by': 'blue', '11ek': 'brown', '11fe': 'green', '11iv': 'magenta', '12cg': 'orange', '13dy': 'purple', '14j': 'black'}
markers = {'92a': '+', '11by': '^', '11ek': 'p', '11fe': 'o', '11iv': '*', '12cg': '>', '13dy': 'D', '14j': 's'}
font = {'family':'serif' , 'color' :'black' , 'size':10,} 

data = {}
for s in sn:
	temp = Table.read('norm_flux/' + s + '/eqw/com.flm', format = 'ascii')
	ions = temp['ions']
	try:
		ews = temp['EQWS']
	except:
		ews = np.zeros(len(ions))
	try:
		errs = temp['EQW_Errs']
	except:
		errs = np.zeros(len(ions))
	try:
		lims = temp['Limits']
	except:
		lims = np.zeros(len(ions))
		
	d = {}
	for i in range(len(ions)):
		d.update({ions[i]: [ews[i], errs[i], lims[i]]})
	
	data.update({s: d})

ions = data['13dy'].keys()
ions.sort()
plots = ions
temp = {}
j = 0
for i in ions:
	temp.update({i:j})
	j += 1
ions = temp

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
			
		ew = data[s][iname][0]
		if (ew == 0):
			ew = data[s][iname][2]
		err = data[s][iname][1]
		
		plt.hlines(ew, -15, 45, color = colors[s])
		plt.vlines(15, ew - err, ew + err, color = colors[s])
		
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


"""

fig = plt.figure(dpi = 100, figsize = [12, 6], facecolor = 'w')
ax = fig.add_subplot(111)

for s in sn:
	d = data[s]
	keys = d.keys()
	for k in range(len(keys)):
		if keys[k] == 'MgI':
			keys[k] = 'MgIa'
	keys.sort()
	
	ews = []
	errs = []
	lims = []
	for k in keys:
		if ((s != '11by') & (s != '13dy')):
			if k == 'MgIa':
				k = 'MgI'
		ews.append(d[k][0])
		errs.append(d[k][1])
		lims.append(d[k][2])
	
	ion = []
	for k in keys:
		ion.append(ions[k])
		
	label = s
	if (s == '92a'):
		label = '92A'
	if (s == '14j'):
		label = '14J'
		
	plt.scatter(ion, ews, color = colors[s], marker = markers[s], label = label)
	plt.errorbar(ion, ews, yerr = [errs, errs], color = colors[s])
	plt.scatter(ion, lims, color = colors[s], marker = markers[s])
	

labels = [item.get_text() for item in ax.get_xticklabels()]
j = 1
for k in keys:
	labels[j] = k
	j += 1
ax.set_xticklabels(labels)
legend = plt.legend(bbox_to_anchor = (.85, 1), loc = 2, shadow = False, numpoints = 1, prop = {'family': 'serif'})
plt.ylabel('EQW ($\AA$)', fontdict = font)
plt.xlabel('Ions', fontdict = font)
plt.show()
"""