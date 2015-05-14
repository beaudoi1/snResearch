import spAnalyze as a

feiia = {}
feiib = {}
feiic = {}
feiid = {}
feiie = {}
mgiia = {}
mgiib = {}
mgi = {}


a.pickSN('92a')
s92a = a.loadData('sn1992a/Data')
temp = {}
for s in s92a:
	temp.update({s.epoch: [s.sig3[0], [0,0]]})
feiia.update({'92a': temp})
temp = {}
for s in s92a:
	temp.update({s.epoch: [s.sig3[1], [0,0]]})
feiib.update({'92a': temp})
temp = {}
for s in s92a:
	temp.update({s.epoch: [s.sig3[2], [0,0]]})
feiic.update({'92a': temp})
temp = {}
for s in s92a:
	temp.update({s.epoch: [s.sig3[3], [0,0]]})
feiid.update({'92a': temp})
temp = {}
for s in s92a:
	temp.update({s.epoch: [s.sig3[4], [0,0]]})
feiie.update({'92a': temp})
temp = {}
for s in s92a:
	temp.update({s.epoch: [s.sig3[5], [0,0]]})
mgiia.update({'92a': temp})
temp = {}
for s in s92a:
	temp.update({s.epoch: [s.sig3[6], [0,0]]})
mgiib.update({'92a': temp})
temp = {}
for s in s92a:
	temp.update({s.epoch: [s.sig3[7], [0,0]]})
mgi.update({'92a': temp})


a.pickSN('11by')
s11by = a.loadData('sn2011by/Data')
temp = {}
for s in s11by:
	temp.update({s.epoch: [s.eqws[0], [s.err_lo[0],s.err_hi[0]]]})
feiia.update({'11by': temp})
temp = {}
for s in s11by:
	temp.update({s.epoch: [s.eqws[1], [s.err_lo[1],s.err_hi[1]]]})
feiib.update({'11by': temp})
temp = {}
for s in s11by:
	temp.update({s.epoch: [s.eqws[2], [s.err_lo[2],s.err_hi[2]]]})
feiic.update({'11by': temp})
temp = {}
for s in s11by:
	temp.update({s.epoch: [s.eqws[3], [s.err_lo[3],s.err_hi[3]]]})
feiid.update({'11by': temp})
temp = {}
for s in s11by:
	temp.update({s.epoch: [s.eqws[4], [s.err_lo[4],s.err_hi[4]]]})
feiie.update({'11by': temp})
temp = {}
for s in s11by:
	temp.update({s.epoch: [s.eqws[7], [s.err_lo[7],s.err_hi[7]]]})
mgiia.update({'11by': temp})
temp = {}
for s in s11by:
	temp.update({s.epoch: [s.eqws[8], [s.err_lo[8],s.err_hi[8]]]})
mgiib.update({'11by': temp})
temp = {}
for s in s11by:
	temp.update({s.epoch: [s.eqws[10], [s.err_lo[10],s.err_hi[10]]]})
mgi.update({'11by': temp})


a.pickSN('11ek')
s11ek = a.loadData('sn2011ek/Data')
temp = {}
for s in s11ek:
	temp.update({s.epoch: [s.sig3[0], [0,0]]})
feiia.update({'11ek': temp})
temp = {}
for s in s11ek:
	temp.update({s.epoch: [s.sig3[1], [0,0]]})
feiib.update({'11ek': temp})
temp = {}
for s in s11ek:
	temp.update({s.epoch: [s.sig3[2], [0,0]]})
feiic.update({'11ek': temp})
temp = {}
for s in s11ek:
	temp.update({s.epoch: [s.sig3[3], [0,0]]})
feiid.update({'11ek': temp})
temp = {}
for s in s11ek:
	temp.update({s.epoch: [s.sig3[4], [0,0]]})
feiie.update({'11ek': temp})
temp = {}
for s in s11ek:
	temp.update({s.epoch: [s.sig3[5], [0,0]]})
mgiia.update({'11ek': temp})
temp = {}
for s in s11ek:
	temp.update({s.epoch: [s.sig3[6], [0,0]]})
mgiib.update({'11ek': temp})
temp = {}
for s in s11ek:
	temp.update({s.epoch: [s.sig3[7], [0,0]]})
mgi.update({'11ek': temp})


a.pickSN('11fe')
s11fe = a.loadData('sn2011fe/Data')
temp = {}
for s in s11fe:
	temp.update({s.epoch: [s.eqws[0], [s.err_lo[0],s.err_hi[0]]]})
feiia.update({'11fe': temp})
temp = {}
for s in s11fe:
	temp.update({s.epoch: [s.eqws[1], [s.err_lo[1],s.err_hi[1]]]})
feiib.update({'11fe': temp})
temp = {}
for s in s11fe:
	temp.update({s.epoch: [s.eqws[2], [s.err_lo[2],s.err_hi[2]]]})
feiic.update({'11fe': temp})
temp = {}
for s in s11fe:
	temp.update({s.epoch: [s.eqws[3], [s.err_lo[3],s.err_hi[3]]]})
feiid.update({'11fe': temp})
temp = {}
for s in s11fe:
	temp.update({s.epoch: [s.eqws[4], [s.err_lo[4],s.err_hi[4]]]})
feiie.update({'11fe': temp})
temp = {}
for s in s11fe:
	temp.update({s.epoch: [s.eqws[5], [s.err_lo[5],s.err_hi[5]]]})
mgiia.update({'11fe': temp})
temp = {}
for s in s11fe:
	temp.update({s.epoch: [s.eqws[6], [s.err_lo[6],s.err_hi[6]]]})
mgiib.update({'11fe': temp})
temp = {}
for s in s11fe:
	temp.update({s.epoch: [s.eqws[7], [s.err_lo[7],s.err_hi[7]]]})
mgi.update({'11fe': temp})


a.pickSN('11iv')
s11iv = a.loadData('sn2011iv/Data')
temp = {}
for s in s11iv:
	temp.update({s.epoch: [s.sig3[0], [0,0]]})
feiia.update({'11iv': temp})
temp = {}
for s in s11iv:
	temp.update({s.epoch: [s.sig3[1], [0,0]]})
feiib.update({'11iv': temp})
temp = {}
for s in s11iv:
	temp.update({s.epoch: [s.sig3[2], [0,0]]})
feiic.update({'11iv': temp})
temp = {}
for s in s11iv:
	temp.update({s.epoch: [s.sig3[3], [0,0]]})
feiid.update({'11iv': temp})
temp = {}
for s in s11iv:
	temp.update({s.epoch: [s.sig3[4], [0,0]]})
feiie.update({'11iv': temp})
temp = {}
for s in s11iv:
	temp.update({s.epoch: [s.sig3[5], [0,0]]})
mgiia.update({'11iv': temp})
temp = {}
for s in s11iv:
	temp.update({s.epoch: [s.sig3[6], [0,0]]})
mgiib.update({'11iv': temp})
temp = {}
for s in s11iv:
	temp.update({s.epoch: [s.sig3[7], [0,0]]})
mgi.update({'11iv': temp})


a.pickSN('12cg')
s12cg = a.loadData('sn2012cg/Data')
temp = {}
for s in s12cg:
	temp.update({s.epoch: [s.sig3[0], [0,0]]})
feiia.update({'12cg': temp})
temp = {}
for s in s12cg:
	temp.update({s.epoch: [s.sig3[1], [0,0]]})
feiib.update({'12cg': temp})
temp = {}
for s in s12cg:
	temp.update({s.epoch: [s.sig3[2], [0,0]]})
feiic.update({'12cg': temp})
temp = {}
for s in s12cg:
	temp.update({s.epoch: [s.sig3[3], [0,0]]})
feiid.update({'12cg': temp})
temp = {}
for s in s12cg:
	temp.update({s.epoch: [s.sig3[4], [0,0]]})
feiie.update({'12cg': temp})
temp = {}
for s in s12cg:
	temp.update({s.epoch: [s.eqws[7], [s.err_lo[7],s.err_hi[7]]]})
mgiia.update({'12cg': temp})
temp = {}
for s in s12cg:
	temp.update({s.epoch: [s.eqws[8], [s.err_lo[8],s.err_hi[8]]]})
mgiib.update({'12cg': temp})
temp = {}
for s in s12cg:
	temp.update({s.epoch: [s.eqws[9], [s.err_lo[9],s.err_hi[9]]]})
mgi.update({'12cg': temp})


a.pickSN('13dy')
s13dy = a.loadData('sn2013dy/Data')
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[0], [s.err_lo[0],s.err_hi[0]]]})
feiia.update({'13dy': temp})
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[1], [s.err_lo[1],s.err_hi[1]]]})
feiib.update({'13dy': temp})
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[2], [s.err_lo[2],s.err_hi[2]]]})
feiic.update({'13dy': temp})
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[3], [s.err_lo[3],s.err_hi[3]]]})
feiid.update({'13dy': temp})
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[4], [s.err_lo[4],s.err_hi[4]]]})
feiie.update({'13dy': temp})
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[5], [s.err_lo[5],s.err_hi[5]]]})
mgiic = temp
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[6], [s.err_lo[6],s.err_hi[6]]]})
mgiid = temp
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[7], [s.err_lo[7],s.err_hi[7]]]})
mgiia.update({'13dy': temp})
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[8], [s.err_lo[8],s.err_hi[8]]]})
mgiib.update({'13dy': temp})
temp = {}
for s in s13dy:
	temp.update({s.epoch: [s.eqws[10], [s.err_lo[10],s.err_hi[10]]]})
mgi.update({'13dy': temp})


a.pickSN('14j')
s14j = a.loadData('sn2014j/Data')
temp = {}
for s in s14j:
	temp.update({s.epoch: [s.sig3[0], [0,0]]})
feiia.update({'14j': temp})
temp = {}
for s in s14j:
	temp.update({s.epoch: [s.sig3[1], [0,0]]})
feiib.update({'14j': temp})
temp = {}
for s in s14j:
	temp.update({s.epoch: [s.sig3[2], [0,0]]})
feiic.update({'14j': temp})
temp = {}
for s in s14j:
	temp.update({s.epoch: [s.sig3[3], [0,0]]})
feiid.update({'14j': temp})
temp = {}
for s in s14j:
	temp.update({s.epoch: [s.sig3[4], [0,0]]})
feiie.update({'14j': temp})
temp = {}
for s in s14j:
	temp.update({s.epoch: [s.eqws[5], [s.err_lo[5],s.err_hi[5]]]})
mgiia.update({'14j': temp})
temp = {}
for s in s14j:
	temp.update({s.epoch: [s.eqws[6], [s.err_lo[6],s.err_hi[6]]]})
mgiib.update({'14j': temp})
temp = {}
for s in s14j:
	temp.update({s.epoch: [s.eqws[7], [s.err_lo[7],s.err_hi[7]]]})
mgi.update({'14j': temp})

chisqs = {}
chisqs.update({'11by': a.chiSq(s11by)})
chisqs.update({'11fe': a.chiSq(s11fe)})
chisqs.update({'12cg': a.chiSq(s12cg)})
chisqs.update({'13dy': a.chiSq(s13dy)})
chisqs.update({'14j': a.chiSq(s14j)})

spectra = [s92a, s11by, s11ek, s11fe, s11iv, s12cg, s13dy, s14j]