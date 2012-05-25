#!/opt/local/bin/python2.7

import biotools.analysis.predict as genepredict
import biotools.analysis.options as options
import threading, sys, os
import Queue as queue

def _run_genepredict(q, infile, names):
	while 1:
		try: strainf = q.get(False)
		except: break
		strain = strainf.split(os.sep)[-1]
		pos = strain.rfind('.')
		if pos > 1 or (pos == 1 and strain[0] != '.'):
			strain = strain[:pos]

		options.debug("Predicting for %s." % strain)

		try: genepredict.run(infile, strainf, strain, names)
		except RuntimeError: pass
		q.task_done()

def run(infile, strains):
	'''run()
Run several instances of genepredict.run at once.'''

	q = queue.Queue()
	for strain in strains:
		q.put(strain)

	filenames = []
	for i in range(options.NUM_PROCESSES-1):
		curr = threading.Thread(target=_run_genepredict, args=(q, infile, filenames))
		curr.start()

	_run_genepredict(q, infile, filenames)
	q.join()

	return filenames

if __name__ == "__main__":
	import sys
	try: run(sys.argv[1:])
	except: pass
