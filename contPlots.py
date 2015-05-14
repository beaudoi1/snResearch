import spAnalyze as a
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

a.pickSN('11fe')
spectra = a.loadData('sn2011fe/Data')
for s in spectra:
	if s.visit == 5:
		sp = s
s = sp
param = s.param
baseline_x = []
baseline_y = []
plt.figure()
plt.plot(s.wave, s.norm_flux)
try:
	a.loadInc(s)
except:
	a.setInc(s)
plt.close()
include = s.inc
gaussian = lambda p, x: 1 + p[0] * np.exp(-(x - p[1])**2 / (2 * p[2]**2))

w = s.wave
f = s.norm_flux

start = np.where(w < include[0])
end = np.where(w > include[len(include) - 1])

for s in range(len(start[0])):
	baseline_x.append(w[start[0][s]])
	baseline_y.append(1)
	
s = sp	
	
i = 0
p = 0

while (i < len(include)):
	inc = np.where((w > include[i]) & (w < include[i + 1]))
	wave = w[inc]
	flux = f[inc]
	x = np.array(wave)
	y = np.array(flux)
	hix = np.linspace(min(x), max(x), len(x) * 10)
	gaus = gaussian([param[p], param[p + 1], param[p + 2]], hix)
	for j in range(len(gaus)):
		baseline_x.append(hix[j])
		baseline_y.append(gaus[j])
	if (i < len(include) - 2):
		ex = np.where((w > include[i + 1]) & (w < include[i + 2]))
		ex_wave = w[ex]
		for k in range(len(ex_wave)):
			baseline_x.append(ex_wave[k])
			baseline_y.append(1)
	i += 2
	p += 3

for e in range(len(end[0])):
	baseline_x.append(w[end[0][e]])
	baseline_y.append(1)
	
font = {'family':'serif' , 'color' :'black' , 'size':24,} 
	
fig = plt.figure(dpi = 100, figsize = [12, 12], facecolor = 'w')
ax = fig.add_subplot(211)
scale = 4e13
plt.plot(s.wave, s.flux * scale, label = 'SN 2011fe', color = 'k', linewidth = 2, linestyle = 'steps')
plt.plot(s.wave, s.sflux * scale, label = 'Continuum', color = 'b', linewidth = 2)
plt.text(2778.8, .776293, 'Mg II', fontsize = 18, color = 'k', family = 'serif')
plt.text(2830.6, 4.98202, 'Mg I', fontsize = 18, color = 'k', family = 'serif')
plt.text(2460.1, -1.2154, 'Fe II', fontsize = 18, color = 'k', family = 'serif')
plt.hlines(-.212375, 2345.2, 2602, lw = 2)
plt.vlines(2345.2, -.212375, .729047, lw = 2)
plt.vlines(2375.64, -.212375, .944216, lw = 2)
plt.vlines(2383.78, -.212375, .68294, lw = 2)
plt.vlines(2587.71, -.212375, 1.29915, lw = 2)
plt.vlines(2602, -.212375, .779957, lw = 2)
plt.hlines(1.61664, 2797.74, 2805.17, lw = 2)
plt.vlines(2797.74, 1.61664, 2.24505, lw = 2)
plt.vlines(2805.17, 1.61664, 2.56688, lw = 2)
plt.vlines(2853.65, 5.78679, 6.75, lw = 2)
plt.ylabel('Scaled f$_\lambda$', fontdict = font)
plt.xlim(2290, 2880)
plt.ylim(-1.5, 11)
legend = plt.legend(loc = 'upper left', shadow = False, prop = {'family':'serif', 'size':24,})
plt.gca().axes.get_xaxis().set_visible(False)
majorLocator = MultipleLocator(100)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(20)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(2)
#majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(.25)
ax.yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
for tick in ax.yaxis.get_major_ticks():
	tick.label.set_fontsize(16)
ax = fig.add_subplot(212)
plt.plot(s.wave, s.norm_flux, color = 'k', linestyle = 'steps', linewidth = 2)
for b in range(len(baseline_x)):
	baseline_x[b] -= .5 
plt.plot(baseline_x, baseline_y, color = 'b', linestyle = 'steps', linewidth = 2)
plt.ylabel('Normalized Flux', fontdict = font)
plt.xlabel('Observed Wavelength ($\AA$)', fontdict = font)
plt.xlim(2290, 2880)
plt.ylim(-.1, 1.2)
majorLocator = MultipleLocator(100)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(20)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(.2)
#majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(.025)
ax.yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
plt.subplots_adjust(hspace = 0, right = .97, left = 0.08, top = 0.97, bottom = 0.1)
#inc1 = np.where((s.wave > 2780) & (s.wave < 2880))
#bx = np.array(baseline_x)
#by = np.array(baseline_y)
#inc2 = np.where((bx > 2780) & (bx < 2880))
a = plt.axes([.17, .14, .35, .13], axisbg = 'y')
#a.set_xticklabels([])
#a.set_yticklabels([])
font = {'family':'serif' , 'color' :'black' , 'size':24,} 
#plt.xlabel('Observed Wavelength ($\AA$)', fontdict = {'family':'serif' , 'color' :'black' , 'size':10,})
#plt.ylabel('Normalized Flux', fontdict = {'family':'serif' , 'color' :'black' , 'size':10,})
majorLocator = MultipleLocator(25)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(5)
a.xaxis.set_major_locator(majorLocator)
a.xaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(.2)
minorLocator = MultipleLocator(.025)
a.yaxis.set_major_locator(majorLocator)
a.yaxis.set_minor_locator(minorLocator)
plt.plot(s.wave, s.norm_flux, color = 'k', linestyle = 'steps', linewidth = 2)
plt.plot(baseline_x, baseline_y, color = 'b', linestyle = 'steps', linewidth = 2)
a.set_axis_bgcolor('white')
plt.xlim(2780, 2880)
plt.ylim(.2, 1.2)
for tick in ax.xaxis.get_major_ticks():
	tick.label.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
	tick.label.set_fontsize(16)
#plt.setp(a, plt.xlim(2780, 2880))
plt.savefig('cont_fit.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
plt.show()
"""
plt.xlim(2740, 2910)
plt.ylim(.25, 1.25)
legend = plt.legend(loc = 'upper left', shadow = False, prop = {'family':'serif'})
majorLocator = MultipleLocator(50)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(10)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(.1)
#majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(.02)
ax.yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
plt.savefig('cont_fit.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
plt.show()



plt.savefig('cont_fit.png', dpi = 800, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
plt.show()
"""