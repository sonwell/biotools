import numpy as np
import matplotlib.pyplot as plt
import math, os

def smoothed(unsmoothed, factor):
	l  = len(unsmoothed)
	r  = [math.exp(-(i-l)**2/float(factor)) for i in xrange(1,l)]
	w  = np.array(r + [1] + r[::-1])
	c  = (w[l-i-1:2*l-i-1] for i in xrange(l))
	wm = np.array([x/sum(x) for x in c])
	return np.dot(unsmoothed, wm.transpose())

def plot(snpdata, directory, bottom=True, side=True, legend=True, filename='untitled.pdf', **kwargs):
	sep = os.sep
	try: os.mkdir(directory)
	except: pass

	upperbound = 0.05
	lowerbound = -upperbound/6

	snt = smoothed(snpdata['nt']['var'], 21)
	lnt =      len(snpdata['nt']['var'])
	saa = smoothed(snpdata['aa']['var'], 21)
	laa =      len(snpdata['aa']['var'])

	xnt = np.arange(lnt)*(0.0+laa)/lnt+1
	xaa = np.arange(laa)              +1

	x = 0.09 if side   else 0.02
	y = 0.09 if bottom else 0.04

	fig = plt.figure(1, facecolor='w', edgecolor='w')
	ax = fig.add_axes([x,y,0.98-x,0.98-y], xlim=[0, laa*1.06], ylim=[lowerbound, upperbound])
	ax.minorticks_on()
	ax.tick_params(axis='x', which='minor', length=3)

	for loc, spine in ax.spines.iteritems():
		if loc in ['right', 'top']:
			spine.set_color('none')
			continue
		if loc == 'bottom':
			ax.xaxis.set_ticks_position('bottom')
			if bottom:
				ax.set_xlabel("Amino acids")
			continue
		if loc == 'left':
			if side:
				ax.set_ylabel("Sequence variance")
				ax.yaxis.set_ticks_position('left')
			else:
				spine.set_color('none')
				ax.tick_params('y', which='both', color='none', labelcolor='none')

	ax.hlines(ax.get_yticks(), 0, laa*1.06, color='0.75', linestyle='dashed')

	ax.hlines(0, 0, laa*1.06, color='k', linestyle='solid')

	nt_lines = plt.plot(xnt, snt, color='#0000ff', linestyle='solid')
	aa_lines = plt.plot(xaa, saa, color='#00ff00', linestyle='solid')

	starts = snpdata['aa']['starts']
	ends   = snpdata['aa']['ends']
	counts = snpdata['aa']['count']
	scale  = laa/max(ends)
	ys = (np.arange(len(starts))+1)*lowerbound/3

	ax.hlines(ys, starts, ends, colors='k', lw=4, linestyle='solid')
	for i, c in zip(xrange(len(counts)), counts):
		ax.text(laa + 10, lowerbound/3*(i+1.25), c)
	if legend:
		fig.legend((nt_lines, aa_lines), ('Nucleotide', 'Amino acid'), 'upper right')

	fig.savefig(directory+filename)

	print '=============', filename, '============='
	print 'Average variance: '
	print '\t', sum(snpdata['nt']['var'])/lnt, 'per base pair'
	print '\t', sum(snpdata['aa']['var'])/laa, 'per amino acid'
