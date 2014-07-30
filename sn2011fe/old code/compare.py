from astropy.table import Table
import glob
import matplotlib.pyplot as plt

files = glob.glob("Data/residuals/*.flm")

waves = []
fluxes = []
for i in range(len(files)):
	spectrum = Table.read(files[i],format='ascii')
	waves.append(spectrum["col1"])
	fluxes.append(spectrum["col2"])
	
plt.figure(1)
for i in range(len(waves)):
	plt.plot(waves[i],fluxes[i], label=files[i])

legend=plt.legend(loc='upper left', shadow=True)
plt.ylim(0,1.5)
plt.show()