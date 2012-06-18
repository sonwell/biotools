#!/usr/bin/python
import biotools.IO               as io
import biotools.sequence         as sequ
import biotools.annotation       as anno
import biotools.BLAST            as BLAST
import biotools.analysis.options as options
from biotools.align      import OptimalCTether as align
from biotools.translate  import *
from biotools.complement import *

try: import Queue as queue
except: import queue

import threading
import random
import sys
import os

PIPING = True

def ORFGenerator(sequ):
	'''ORFGenerator(sequ)
Scans both strands of the give sequence and yields the longest subsequence that starts with a start codon and contains no stop codon other than the final codon.'''

	comp   = complement(sequ[::-1])
	seq    = sequ.seq
	cseq   = comp.seq
	slen   = len(sequ)
	starts = [-1, 0, 1, -1, 0, 1]  # locations of start codons in each frame
	stops  = [ 0, 1, 2,  0, 1, 2]  # locations of stop  codons in each frame

	for i in xrange(slen-2):
		fcodon, rcodon = seq[i:i+3], cseq[i:i+3]
		if fcodon in options.STOP_CODONS:
			if starts[i%3+3] >= stops[i%3+3] and \
			   i - starts[i%3+3] >= options.MIN_ORFLEN:
				yield sequ[starts[i%3+3]:i+3]
			stops[i%3+3] = i+3
		elif fcodon in options.START_CODONS:
			if starts[i%3+3] < stops[i%3+3]:
				starts[i%3+3] = i
		if rcodon in options.STOP_CODONS:
			if starts[i%3+0] >= stops[i%3+0] and \
			   i - starts[i%3+0] >= options.MIN_ORFLEN:
				yield comp[starts[i%3+0]:i+3]
			stops[i%3+0] = i+3
		elif rcodon in options.START_CODONS:
			if starts[i%3+0] < stops[i%3+0]:
				starts[i%3+0] = i

	raise StopIteration
	
def _genepredict_target(qin,qout,orfs,subj):
	while 1:
		try: res = qin.get(PIPING, 1)
		except queue.Empty:
			if not PIPING:
				break
			else:
				continue

		qname, sname = res['query']['name'],  res['subject']['name']
		start, end   = res['query']['start'], res['query']['end']

		if subj[sname].type == 'nucl': subject = translate(subj[sname])
		else: subject = subj[sname]
		length = len(subject)*3

		while qname:
			try:
				o = orfs[qname]
				break
			except KeyError:
				qname = qname[:-1]
		if not qname:
			qin.task_done()
			continue

		alignments = []
		for orf in o:
			if _genepredict_in_range(orf, start, end, res['frame']):
				orf = orf[:-3]
				query = translate(orf)
				options.debug("Aligning %33s v. %33s." % (qname, sname))

				alignments.append((orf,sname,align(subject.seq,query.seq)))
			
		max_match = (options.MIN_IDENTITY,None)
		for orf, refname, alignment in alignments:
			hitlen = alignment['sublength']
			identity = float(alignment['identities'])/alignment['length']
			if identity >= max_match[0]:
				max_match = (identity,(orf[-3*hitlen:], sname, alignment))
		if max_match[1]:
			seq, name, _ = max_match[1]
			odl = subject.defline.split('[')[0].strip()
			src = seq.original.name
			start, end, strand = seq.start, seq.end, seq.step
			defline = '%s[source=%s] [start=%d] [end=%d] [strand=%d]' % \
				(odl + (' ' if odl else ''), src, start, end, strand)
				
			qout.put(sequ.Sequence(name.strip(), seq.seq, defline=defline,
				original=seq.original, type=seq.type,
				start=seq.start, end=seq.end, step=seq.step))

		qin.task_done()

def _genepredict_in_range(seq, start, end, frame):
	ss, se = sorted((seq.start, seq.end))
	os, oe = sorted((start, end))
	frame  = int(frame)

	return ss <= os and se >= oe and \
	     (se%3 == oe%3 or ss%3 == oe%3) and \
	     ((frame < 0 and seq.step < 0) or \
	      (frame > 0 and seq.step > 0))

class ThreadQueue(queue.Queue):
	def __init__(self, qout, orfs, subj):
		queue.Queue.__init__(self)
		self.threadcount = 0
		self.qout = qout
		self.orfs = orfs
		self.subj = subj

	def put(self, item):
		options.lock.acquire()
		queue.Queue.put(self, item)
		if self.threadcount < options.NUM_THREADS - 1:
			thread = threading.Thread(target=_genepredict_target,
				args=(self, self.qout, self.orfs, self.subj))
			thread.start()
			self.threadcount += 1
		options.lock.release()

def GeneFromBLAST(db, sequences, pref, names):
	'''GeneFromBLAST(database, sequences, prefix)
BLASTs database against sequences, and for those results that pass the length and percent identity requirements, attempt to locate the full gene that corresponds to that BLAST hit. Genes that are found are saved in the subdirectory sequences under the given directory, divided depending on whether the sequnece is amino acid or nucleotide.'''
	global PIPING

	sep = os.sep
	wd = options.DIRECTORY + 'sequences' + sep
	try: os.mkdir(options.DIRECTORY)
	except: pass
	try: os.mkdir(wd)
	except: pass
	try: os.mkdir(wd + 'nt' + sep)
	except: pass
	try: os.mkdir(wd + 'aa' + sep)
	except: pass

	subj = dict((s.name, s) for s in io.open(db, 'r'))
	options.debug("Database sequences loaded from file %s." % db)

	try:
		orfs = dict((s.name, [orf for orf in ORFGenerator(s)]) \
			for s in io.open(sequences, 'r'))
		options.debug("ORFs loaded from file %s." % sequences)
	except IOError:
		options.debug("No file \"" + sequences + ",\" skipping.")
		return

	q_outputs = queue.Queue()
	q_inputs  = ThreadQueue(q_outputs, orfs, subj)

	blastopts = {
		'evalue': options.MAX_EVALUE,
		'num_threads': options.NUM_THREADS
	}
	
	for res in BLAST.run(db, sequences, **blastopts): 
		if float(res['expect']) > options.MAX_EVALUE:
			continue

		sbjl  = len(subj[res['subject']['name']])
		ident = float(res['identities'].split('(')[1][:-2]) / 100
		lerr  = float(res['subject']['length'])/sbjl
		
		if ident >= options.MIN_IDENTITY:
			if lerr >= (1.0-options.LENGTH_ERR):
				q_inputs.put(res)
	PIPING = False
	options.debug("BLAST done.")
	
	_genepredict_target(q_inputs,q_outputs,orfs,subj)
	q_inputs.join()

	seqs = {}
	nuc_file = io.open(wd + pref + '.fasta', 'w')
	while not q_outputs.empty():
		try:
			seq = q_outputs.get(False)
			if seq.seq not in seqs:
				seqs[seq.seq] = set()
			seqs[seq.seq].add(seq)
			nuc_file.write(seq)
		except: break
	nuc_file.close()

	gh = io.open(wd + pref + '.gff3', 'w' )
	names.append(wd + pref + '.fasta')

	sa = sequ.annotation
	for id in seqs:
		gh.write(sa(seqs[id].copy().pop(), pref, 'gene', 
			homologs=','.join(s.name for s in seqs[id])))
	gh.close()

def run(subject, query, prefix, names):
	GeneFromBLAST(subject, query, prefix, names)

