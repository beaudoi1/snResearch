#sp has attributes wave, norm_flux, 
sp = a.spectrum()

gaussian = lambda p, w, x: p[0]*np.exp(-(x - w*(1 + p[2]))**2 / (2*p[1]**2))
	
norm_flux = sp.norm_flux

sigma = np.linspace(2.5, 4, 16)
# [2.5, 2.6, 2.7, 2.8, ..., 4]
height = np.linspace(-.02, -.6, 30)
# [-.02, -.04, 0.06, ..., -.6]
red = np.linspace(.01, .03, 9)
# [.01, .0125, .015, ..., .03]

# reference wavelength to move with redshift
rest = 2800

ew_meas = []
ew_guess = []
ew_err = []
for s in sigma:
	for h in height:
		for r in red:
			# makes a new gaussian with a set of params
			gauss = gaussian([h, s, r], rest, sp.wave)
			
			# finds wavelength range for gaussian
			xindecies = np.where(gauss < -.001)
			xindecies = xindecies[0]
			xmin = sp.wave[min(xindecies)] - 1
			xmax = sp.wave[max(xindecies)] + 1
			
			# subtracts Gaussian line from normalized flux
			flux = sp.norm_flux + gauss
			# Finds EW value for ew_in
			ew_guess.append(np.abs(integrate.quad(lambda x: h * np.exp(-(x-rest*(1+r))**2/(2*s**2)), min(sp.wave), max(sp.wave))[0]))
			
			# Find EW value for ew_out
			sp.norm_flux = flux
			sp.inc = [xmin, xmax]
			guesses = [h, rest*(1+r), s]
			findEqw(sp, False)
			
			ew_meas.append(sp.eqws[0])
			ew_err.append(sp.eqw_errs[0])
			sp.norm_flux = norm_flux

# gets rid of ew_in = 0 (ew_err also is 0 for these values)
guess = []
meas = []
err = []
for i in range(len(ew_guess)):
	if (ew_guess[i] > .5):
		guess.append(ew_guess[i])
		meas.append(ew_meas[i])
		err.append(ew_err[i])

plt.figure()
plt.scatter(guess, meas)
fit = np.polyfit(guess, meas, 1)
x = np.linspace(0, 7, 8)
y = fit[0] * x + fit[1]
plt.plot(x, y, linewidth = 3, color = 'r')
plt.show()

# SNR: ew_out / ew_err
snr = []
for i in range(len(meas)):
	snr.append(meas[i]/err[i])
	
plt.figure()
plt.scatter(guess, snr)
plt.show()