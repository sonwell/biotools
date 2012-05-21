import biotools.analysis.predict as pr
import biotools.IO               as io

import os

directory = '/Volumes/Stasklab Server/Cassava Group Tools/Contigs/'

for sequences in (f for f in os.listdir(directory) if f.endswith('.fa')):
	orfs = [orf for orf in pr.ORFGenerator(s) for s in io.open(sequences, 'r')]

	fh = io.open(sequences[:-3] + ' orfs.fa', 'w')

	for orf in orfs:
		orf.name = orf.original.name
		orf.defline = '[start=%d] [end=%d] [strand=%d]' % (orf.start, orf.end, orf.step)
		fh.write(orf)

	fh.close()
