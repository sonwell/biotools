from optparse  import OptionParser
from threading import Lock
import os

LENGTH_ERR    = 0.2
MIN_IDENTITY  = 0.45
MAX_EVALUE    = 1e-30
MIN_ORFLEN    = 300
NUM_THREADS   = 16
NUM_PROCESSES = 2
DIRECTORY     = '.' + os.sep
PLOTTER       = 'biotools.analysis.plot'

START_CODONS = ["ATG"]
STOP_CODONS  = ["TAG", "TAA", "TGA"]
args = tuple()

predicting  = True
clustering  = True
renaming    = True
calculating = True
reporting   = True
plotting    = True
verbose     = False

lock = Lock()

parser = OptionParser(usage = "Usage: %prog [options] <database> <sequences ...>")
parser.add_option("-S", "--start", action="append", dest="start",
	default=START_CODONS, help="define a start codon [default: %s]" % 
	' '.join( "-S " + s for s in START_CODONS), type="string")
parser.add_option("-E", "--stop", action="append", dest="stop",
	default=STOP_CODONS, help="define a stop codon [default: %s]" %
	' '.join( "-E " + s for s in STOP_CODONS), type="string")
parser.add_option("-j", "--threads", action="store",	dest="threads",
	default=NUM_THREADS, help="number of threads [default: %default]", type="int")
parser.add_option("-p", "--processes", action="store", dest="processes",
	default=NUM_PROCESSES, help="number of parallel processes to run [default: %default]", type="int")
parser.add_option("-e", "--evalue", action="store",	dest="evalue",	 
	default=MAX_EVALUE, help="maximum e-value [default: %default]", type="float")
parser.add_option("-I", "--identity", action="store",	dest="identity", 
	default=MIN_IDENTITY, help="minimum percent identity [default: %default]", type="float")
parser.add_option("-L", "--length",	action="store",	dest="fraction",	 
	default=LENGTH_ERR, help="allowable relative error in hit length (|hit-ref|/ref) [default: %default]", type="float")
parser.add_option("-O", "--orflen", action="store", dest="orflen", metavar="bases",
  default=MIN_ORFLEN, help="minimum allowable length for ORFs [default: %default]", type="int")
parser.add_option("-d", "--directory", action="store", dest="directory",	 
	default=DIRECTORY, help="minimum percent length [default: current]", type="string")
parser.add_option("-P", "--plotter", action="store", dest="plotter",
	default=PLOTTER, help="plotting module [default: %default]", type="string")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
	default=verbose, help="print debug messages [default: False]")

parser.add_option("--no-plots", action="store_false", dest="plotting",
	default=plotting, help="suppress the drawing of plots [default: False]")
parser.add_option("--no-predict", action="store_false", dest="predicting",
	default=predicting, help="don't predict genes, instead treat the input files as predicted genes [default: False]")
parser.add_option("--no-cluster", action="store_false", dest="clustering",
	default=True, help="don't cluster the sequences, instead treat the input files as alignments [default: False]")
parser.add_option("--no-rename", action="store_false", dest="renaming",
	default=True, help="don't rename the fasta and clustal files [default: False]")
parser.add_option("--no-reports", action="store_false", dest="reporting",
	default=True, help="don't generate files for variance data [default: False]")
parser.add_option("--no-calculation", action="store_false", dest="calculating",
	default=True, help="don't calculate sequence variance [default: False]")


def debug(msg):
	if verbose:
		lock.acquire(True)
		print msg
		lock.release()

def parse(pargs):
	'''parse(pargs)
Parses pargs and sets the various variables to be accessible to other modules.

These variables are:
* LENGTH_ERR
* MIN_IDENTITY
* MAX_EVALUE
* NUM_THREADS
* NUM_PROCESSES
* START_CODONS
* START_CODONS
* DIRECTORY
* PLOTTER
* args'''
	global LENGTH_ERR,    \
	       MIN_IDENTITY,  \
	       MAX_EVALUE,    \
	       MIN_ORFLEN,    \
	       NUM_THREADS,   \
	       NUM_PROCESSES, \
	       START_CODONS,  \
	       STOP_CODONS,   \
	       DIRECTORY,     \
	       PLOTTER,       \
	       args,          \
	       predicting,    \
	       clustering,    \
	       renaming,      \
	       calculating,   \
	       reporting,     \
	       plotting,      \
	       verbose

	opts, largs   = parser.parse_args(pargs)

	if opts.directory[-1] != os.sep:
		opts.directory += os.sep
	try:    os.makedirs(opts.directory)
	except: pass

	if '*' in opts.start:
		DNA = 'ATCG'
		opts.start = set(i+j+k for i in DNA \
		                       for j in DNA \
		                       for k in DNA) - \
		             set(opts.stop)

	if not (0 <= opts.fraction <= 1):
		raise RuntimeError, "Allowable length error must be between 0 and 1, inclusive."

	LENGTH_ERR    = opts.fraction
	MIN_IDENTITY  = opts.identity
	MAX_EVALUE    = opts.evalue
	MIN_ORFLEN    = opts.orflen
	NUM_THREADS   = opts.threads
	NUM_PROCESSES = opts.processes
	STOP_CODONS   = opts.stop
	START_CODONS  = opts.start
	DIRECTORY     = opts.directory
	PLOTTER				= opts.plotter

	predicting    = opts.predicting
	clustering    = opts.clustering
	renaming      = opts.renaming
	calculating   = opts.calculating
	reporting     = opts.reporting
	plotting      = opts.plotting

	verbose       = opts.verbose

	args = largs

def help():
	'''help()
Prints the usage.'''
	parser.print_help()
