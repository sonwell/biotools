#!/usr/bin/python
from biotools.iomanager import IOManager
import __builtin__

def _parseAttrs(attr,token='='):
  '''_parseAttrs(attr, token) /internal/
Creates a dictionary from the atrributes (9th column) of a gff file. By default, token is '=', which is the separator used in gff version 3.

In other words, attr "a=b;c=d;" and token '=' will yield the dictionary {'a':'b','c':'d'}. The other separator (';') cannot be changed.

This function is not to be called on its own.'''

  attributes = {}
  attrs = [a.strip() for a in attr.strip().split(';')]
  for attribute in attrs:
    pos = attribute.find(token)
    if pos > -1:
      var,val = attribute[:pos],attribute[pos+1:]
      attributes[var] = attributes.get(var,[]) + [val]
	
  for key in attributes:
    attributes[key] = ','.join(attributes[key])
  return attributes

class Annotation(object):
  '''class Annotation
An object to help with reading and writing GFF files.'''
  unknowns = 0
	
  def __init__(self,ref,src,type,start,end,score,strand,phase,attr,name_token='ID',gff_token='='):
    '''Annotation(ref, src, type, start, end, score,
           strand, phase, attr, name_token, gff_token)
Constructs an Annotation object with the necessary values. The parameters are passed in the same order as the columns from a GFF (version 3) file and the name_token and gff_token parameters are the defaults for a gff version 3 file from phytozome. Just write (e.g.)
	Annotation(*line.split('\\t')) #(splitting on tabs),
and the rest of the work will be done for you. Other sources may require changes to name_tokens and gff_token.

Instantiating an Annotation will generate for it an id of the form SEQNAME_TYPE[START:END], where SEQNAME is the name of the sequence (column 1) from the GFF file, and type is like 'gene' or 'CDS'. If no SEQNAME is provided, then 'X' be used in its place, and if no identifier can be found in the attributes, the Annotation will generate an identifier for itself in the form of unknown #.'''

    start,end = int(start),int(end)
    self.strand = strand
    self.type   = type
    self.source = src
    self.seq    = ref
    self.start  = min(start,end)
    self.end    = max(end,start)
    self.attr   = _parseAttrs(attr,gff_token)
    self.phase  = phase
    self.score  = score
    self.ntoken = name_token
    self.id     = ((self.seq or 'X') + '_' + self.type + "[%d:%d]" % (self.start,self.end))
    try:
      self.name   = self.attr[name_token]
    except:
      Annotation.unknowns += 1
      self.name = "unknown %d" % Annotation.unknowns
    self.parent = None
    self.children = []

  '''Some things that you can do to Annotation objects:
* len(annotation)  => gets the length of the annotation (end-start+1)
* dictionary[annotation] => store annotations as keys of a dictionary or as elements in a set
* annA == annB     => compare two Annotations, they are the same if they have the same id.
* print annotation => prints the annotation as a line of a GFF version 3 file.
'''
		
  def __len__(self):
		return max(self.start,self.end)-min(self.end,self.start)+1
    
  def __hash__(self):
    return self.id.__hash__()
    
  def __eq__(self,other):
    try:
      return self.id == other.id
    except:
      return False

  def __str__(self):
    return '\t'.join((self.seq,self.source, \
      self.type,str(self.start),str(self.end),self.score, \
      self.strand,str(self.phase), \
      ';'.join(k+'='+self.attr[k] for k in self.attr)))

def _io_methods():
	def read_gff(fh):
		for line in fh:
			if line[0] != '#':
				yield Annotation(*line.split('\t'))
		raise StopIteration

	def write_gff(fh, a):
		fh.write(str(a) + '\n')

	def probe_gff(fh):
		for line in fh:
			line = line.strip()
			if line:
				bits = line.split()
				if bits[0] == '##gff-version':
					return {'anntype': 'gff', 
									'version': float(bits[1])}
				return False
		return False

	nil = lambda x: None

	def whook_gff(fh):
		fh.write('##gff-version 3\n')

	return {
		'gff': {'rhook': nil,       'read': read_gff,
	          'whook': whook_gff, 'write': write_gff,
	                              'probe': probe_gff}
	}

class IOBase(object):
	'''class IOBase(object)
Generic IO class for sequence files.'''
	methods = IOManager(_io_methods())
		
	def __init__(self, name, mode):
		'''IOBase(name, mode)
Opens file name with mode mode. This function will attempt to guess at the filetype by 1. looking at the file extension and failing that, will 2. read the first few lines to determine the file type.

Recoginized file extensions include fa, fsa, fas, fasta, fastc, fastq, clustalw, clustal, aln.'''

		self.handle = __builtin__.open(name, mode)
		self.method = self.methods.default
		self.anntype = None
		self.haswritten = False

		suffixes = {'gff': 'gff', 'gff3':  'gff'}

		p = name.rfind('.')
		if p > -1:
			ext = name[p+1:]
			if ext in suffixes:
				self.seqtype = suffixes[ext]
				self.format(self.seqtype)
				return
		try:
			for method in methods:
				ret = method['probe'](__builtin__.open(name, 'r'))
				if ret:
					for key in ret:
						object.__setitem__(self, key, ret[key])
					break
		except IOError:
			return

	def format(self, fmt):
		'''format(fmt)
Forces a file to be parsed as a particular format. By default, the values for fmt can be fasta, fastq, or fastc.'''

		if fmt in self.methods:
			method = self.methods[fmt]
			ret = method['probe'](__builtin__.open(name, 'r'))
			if ret:
				self.method = method
				for key in ret:
					object.__setitem__(self, key, ret[key])
				return self
			else: raise ValueError, "File cannot be parsed as type %s." % fmt
		self.method = self.methods.default
		return self

	def write(self, sequence):
		'''write(sequence):
Writes sequence as the correct format to the file.'''
		if not self.haswritten:
			self.method['whook'](self.handle)
			self.haswritten = True
		self.method['write'](self.handle, sequence)

	def read(self, count = None):
		'''read(n=None)
If n is provided, the next (up to) n entries are parsed and returned. Otherwise, all remaining entries are parsed and returned.'''
		if count is None: return [s for s in self]
		return [s for s, i in zip(iter(self), xrange(int(count)))]

	def __iter__(self):
		self.method['rhook'](self.handle)
		return self.method['read'](self.handle)

	def next(self):
		'''next()
Reads a single entry in the file and returns it.'''
		try: return self.read(1)[0]
		except ValueError: raise StopIteration

	def close(self):
		'''close()
Close the file handle.'''
		self.handle.close()

class Reader(IOBase):
	'''class Reader(IOBase)
A class that wraps IOBase and restricts the ability to write.'''
	def __init__(self, filename, mode='r'):
		IOBase.__init__(self, filename, mode)

	def write(self, *args):
		raise AttributeError, "Annotation reader cannot write."

class Writer(IOBase):
	'''class Writer(IOBase)
A class that wraps IOBase and restricts the ability to read.'''
	def __init__(self, filename, mode='w'):
		IOBase.__init__(self, filename, mode)
			
	def read(self, *args):
		raise AttributeError, "Annotation writer cannot read."

def open(name, mode='r'):
	'''open(filename, mode='r')
Open a file for parsing or creation. Returns either a Reader or Writer object, depending on the open mode.'''
	if mode == 'r': return Reader(name, mode)
	if mode == 'w': return Writer(name, mode)
	if mode == 'a': return Writer(name, mode)
