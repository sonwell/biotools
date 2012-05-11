from optparse import OptionParser
import os

LENGTH_ERR    = 0.2
MIN_IDENTITY  = 0.45
MAX_EVALUE    = 1e-30
NUM_THREADS   = 16
NUM_PROCESSES = 2
DIRECTORY     = '.' + os.sep
PLOTTER       = 'biotools.analysis.plot'

START_CODONS = ["ATG"]
STOP_CODONS  = ["TAG", "TAA", "TGA"]
args = tuple()

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
parser.add_option("-d", "--directory", action="store", dest="directory",	 
	default=DIRECTORY, help="minimum percent length [default: current]", type="string")
parser.add_option("-P", "--plotter", action="store", dest="plotter",
	default=PLOTTER, help="plotting module [default: %default]", type="string")

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
	       NUM_THREADS,   \
	       NUM_PROCESSES, \
	       START_CODONS,  \
	       STOP_CODONS,   \
	       DIRECTORY,     \
	       PLOTTER,       \
	       args

	opts, largs   = parser.parse_args(pargs)

	if opts.directory[-1] != os.sep:
		opts.directory += os.sep
	try:    os.makedirs(opts.directory)
	except: pass

	if not (0 <= opts.fraction <= 1):
		raise RuntimeError, "Allowable length error must be between 0 and 1, inclusive."

	LENGTH_ERR    = opts.fraction
	MIN_IDENTITY  = opts.identity
	MAX_EVALUE    = opts.evalue
	NUM_THREADS   = opts.threads
	NUM_PROCESSES = opts.processes
	STOP_CODONS   = opts.stop
	START_CODONS  = opts.start
	DIRECTORY     = opts.directory
	PLOTTER				= opts.plotter
	args = largs

def help():
	'''help()
Prints the usage.'''
	parser.print_help()
