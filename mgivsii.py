import orgData as o
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

sn = ['92a', '11by', '11ek', '11fe', '11iv', '12cg', '13dy', '14j']
colors = {'92a': 'black', '11by': 'blue', '11ek': 'cyan', '11fe': 'green', '11iv': 'red', '12cg': 'orange', '13dy': 'purple', '14j': 'grey'}

fig = plt.figure(dpi = 100, figsize = [12,6], facecolor = 'w')
ax = fig.add_subplot(111)

for s in sn:
	a = o.mgiia[s]
	b = o.mgiib[s]
	i = o.mgi[s]
	
	epochs = a.keys()
	
	mgii = []
	mgiierr = []
	mgi = []
	mgierr = []
	
	for e in epochs:
		mgii.append(a[e][0] + b[e][0])
		mgiierr.append([np.sqrt(a[e][1][0]**2 + b[e][1][0]**2), np.sqrt(a[e][1][1]**2 + b[e][1][1]**2)])
		mgi.append(i[e][0])
		mgierr.append([i[e][1][0], i[e][1][1]])
		
	x = []
	y = []
	xlo = []
	xhi = []
	ylo = []
	yhi = []
	xweights = []
	yweights = []
	for j in range(len(mgii)):
		x.append(mgii[j])
		y.append(mgi[j])
		xlo.append(mgiierr[j][0])
		xhi.append(mgiierr[j][1])
		ylo.append(mgierr[j][0])
		yhi.append(mgierr[j][1])
		try:
			xweights.append(1/(np.abs(xlo[j] - xhi[j]) / 2)**2)
		except:
			xweights.append(0)
		try:
			yweights.append(1/(np.abs(ylo[j] - yhi[j]) / 2)**2)
		except:
			yweights.append(0)
		
	if ((xweights[0] != 'inf') & (yweights[0] != 0)):
		xavg = np.average(x, weights = xweights)
		yavg = np.average(y, weights = yweights)		
	else:
		xavg = np.average(x)
		yavg = np.average(y)
	
	
	plt.scatter(x, y, color = colors[s], marker = '.', label = s)	
	if((xlo[0] != 0) & (ylo[0] != 0)):
		plt.errorbar(x, y, xerr = [xlo, xhi], yerr = [ylo,yhi], color = colors[s], linestyle = 'None')

	plt.scatter(xavg, yavg, color = colors[s], marker = 'o')
	
plt.ylabel('MgI')
plt.xlabel('MgII')
legend = plt.legend(bbox_to_anchor = (.85, .6), loc = 2, shadow = False, numpoints = 1, prop = {'family': 'serif'})
majorLocator = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(0.25)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(.25)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
plt.savefig('mgi_vs_mgii.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
plt.show()
