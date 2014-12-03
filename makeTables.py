import orgData as o
from astropy.table import Table

sn = ['92a', '11by', '11ek', '11fe', '11iv', '12cg', '13dy', '14j']
names = {'92a': 'SN~1992A'}
names.update({'11by': 'SN~2011by'})
names.update({'11ek': 'SN~2011ek'})
names.update({'11fe': 'SN~2011fe'})
names.update({'11iv': 'SN~2011iv'})
names.update({'12cg': 'SN~2012cg'})
names.update({'13dy': 'SN~2013dy'})
names.update({'14j': 'SN~2014J'})
dates = {'92a': [19920124.21, 19920304.31, 19921105.59]}
dates.update({'11by': [20110430.56, 20110509.38]})
dates.update({'11ek': [20110812.56, 20110819.45]})
dates.update({'11fe': [20110828.17, 20110831.27, 20110903.43, 20110907.42, 20110910.44, 20110913.68, 20110919.63, 20111001.27, 20111007.32, 20111021.2]})
dates.update({'11iv': [20111211.10, 20111215.94, 20111220.72, 20111224.58, 20111228.44, 20120101.36, 20120109.21]})
dates.update({'12cg': [20120604.50, 20120618.40]})
dates.update({'13dy': [20130721.50, 20130725.61, 20130727.34, 20130729.32, 20130801.59, 20130805.51, 20130809.11, 20130811.15, 20130815.01, 20130817.94]})
dates.update({'14j': [20140126.60, 20140128.44, 20140130.49, 20140201.61, 20140204.73, 20140208.53, 20140210.44, 20140213.29, 20140216.41, 20140226.07]})

data = []
data.append(o.feiia)
data.append(o.feiib)
data.append(o.feiic)
data.append(o.feiid)
data.append(o.feiie)
data.append(o.mgiia)
data.append(o.mgiib)
data.append(o.mgi)

strings = []
b = '$ & $'

for s in sn:
	strings.append('\multicolumn{10}{c}{' + names[s] + '}\\')
	date = dates[s]
	epochs = data[0][s].keys()
	epochs.sort()
	i = 0
	for e in epochs:
		temp = '$' + str(date[i]) + b
		temp += str(round(e,2))
		for d in data:
			ews = d[s][e]
			temp += b + str(round(ews[0], 2))
			if(ews[1][0] != 0):
				temp += '^{+' + str(round(ews[1][1], 2)) + '}_{-' + str(round(ews[1][0], 2)) + '}'
		temp += '$\\'
		strings.append(temp)
		i += 1
	
f = open('data_tables.tex', 'w')
for s in strings:
	f.write(s)
	f.write('\n')
f.close()

values = []
for s in sn:
	epochs = data[0][s].keys()
	epochs.sort()
	for e in epochs:
		for d in data:
			temp = [s, e, d[s][e][0], d[s][e][1][0], d[s][e][1][1]]
			values.append(temp)
			
sns = []
ews = []
lo = []
hi = []
for v in vales:
	sns.append(v[0])
	ews.append(v[1])
	lo.append(v[2])
	hi.append(v[3])
	
t = Table([sns, ews, lo, hi], names = ('SN', 'ews', 'err_lo', 'err_hi',))
t.write('data.flm', format = 'ascii')