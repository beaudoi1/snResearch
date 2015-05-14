sp.res = sp.norm_flux - sp.comb

font = {'family' : 'serif','color'  : 'black','weight' : 'bold','size' : 10,} 

i = 1
for s in spectra:
	fig = plt.figure(num = i, dpi = 100, figsize = [10,8], facecolor = 'w')
	plt.title(s.name, fontdict = font)
	plt.gca().axes.get_xaxis().set_visible(False)
	plt.gca().axes.get_yaxis().set_visible(False)
	ax1 = fig.add_subplot(2, 1, 1)
	ax1.plot(s.wave, s.norm_flux, color = 'k')
	plt.xlim(2739, 2889)
	plt.ylim(.2, 1.4)
	ax1.set_xticklabels([])
	plt.gca().axes.get_yaxis().set_visible(True)
	plt.ylabel('Norm Flux', fontdict = font)
	plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))

	ax2 = fig.add_subplot(2, 1, 2)
	ax2.plot(s.wave, s.res, color = 'b')
	plt.xlim(2739, 2889)
	plt.ylim(-1, 1)
	plt.gca().axes.get_yaxis().set_visible(True)
	plt.ylabel('Residual', fontdict = font)

	plt.gca().axes.get_yaxis().set_visible(True)
	plt.xlabel('Wavelength', fontdict = font)
	plt.subplots_adjust(hspace = 0, right = .8)
	plt.savefig('test_plots/norm_vs_res/norm_vs_res_' + str(i) + '.png', dpi = 600, facecolor = 'w', edgecolor = 'w', pad_inches = .1)
	i += 1