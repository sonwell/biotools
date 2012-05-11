#!/usr/bin/python
import sys, os
import biotools.sequence as sequ
import biotools.analysis.options as options

def SaySNPs(input):
	'''SaySNPs(file)
Takes a clustalw alignment and will return a dictionary of data
relevent to plotting the sequence variance for the sequences in the
given clustalw alignment. These data are:
* var	   = the measure of sequence variation,
* starts = the starting positions for each gene model in amino acids,
* ends   = the ending positions for each gene model in amino acids, and
* count  = the number of sequences with a particular gene model.
The values given in starts, ends, and counts are sorted to that the nth
element in starts corresponds to the nth value in ends and the nth value 
in counts.'''

	SNPcatalogue = []
	lengths = {}
	for seq in sequ.open(input, 'r'):
		key = (seq.start, seq.end)
		lengths[key] = lengths.get(key,0) + 1
		for (i,c) in zip(xrange(key[1]-key[0]+1),seq):
			if i >= len(SNPcatalogue): SNPcatalogue.append({})
			if c != " ": SNPcatalogue[i][c] = SNPcatalogue[i].get(c,0) + 1

	SNPcalc = []
	for s in SNPcatalogue:
		tot = 0.0+sum(s.values())
		cnt = 0.0+len(s)
		SNPcalc.append(1.0-sum((s[c]/tot)**2 for c in s))
	llist = list(lengths.keys())
	llist.sort()

	return {'var':    SNPcalc, 
					'starts': [s for s,e in llist], 
					'ends':   [e for s,e in llist],
					'count': [lengths[k] for k in llist]}

def var(strain, fmt):
	'''var(strain, fmt)
Returns plotdata and metadata for plotting later on in the pipeline.'''
	plotdata = {'nt': SaySNPs( options.DIRECTORY + fmt % {'strain': strain, 'type': 'nt'}),
							'aa': SaySNPs( options.DIRECTORY + fmt % {'strain': strain, 'type': 'aa'})}
	metadata = {'strain': strain, 'filename': strain + '.pdf'}
	
	return plotdata, metadata
