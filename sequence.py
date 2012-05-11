#!/usr/bin/python
from biotools.annotation import Annotation
from biotools.iomanager import IOManager
import __builtin__

def chop(seq,length = 70):
  '''chop( sequence, length = 70 )
Yields a chunk of a sequence of no more than length characters,
it is meant to be used to print fasta files.'''

  while seq:
    try: piece,seq = seq[:length],seq[length:]
    except: piece,seq = seq,''
    yield piece
  raise StopIteration

class Sequence(object):
  '''class Sequence
A wrapper class for sequences.'''

  def __init__(self,name,seq,**kwargs):
    '''Sequence( name, sequence, ... )
Instantiates a Sequence object with a given sequence. Some other useful
parameters that the Sequence constructor can handle are:
  * qual        => the quality scores (an array of integers) of the sequence,
  * type        => the type of the sequence, either prot or nucl,
  * start       => the starting position of the sequence within a supersequence,
  * end         => the ending position of the sequnece within a supersequence,
  * step        => the 'step' of the sequence, usually +1 for top-strand 
                   sequences, and -1 for bottom-strand sequences, but can handle 
                   other values as well,
  * original    => the original Sequence object from which this one derives,
  * defline     => the definition line for this sequnce from a fasta file.
If one of these are not given, they will default to the most logical value
that can be determined from the other values and sequence (e.g., 
if end < start, then step is probably -1).'''

    self.name = name
    self.seq = seq.upper()
    self.qual = kwargs.get('qual',None)
    self.type = kwargs.get('type','prot' if set(seq)-set('ATCGNYR') else 'nucl')
    self.start = kwargs.get('start',1)
    self.end = kwargs.get('end',len(seq))
    self.step = kwargs.get('step', -1 if self.start > self.end else 1)
    self.original = kwargs.get('original',self)
    self.defline = kwargs.get('defline','')

  def __getitem__(self,key):
    '''sequence[i] or sequence[start:end] or sequence[start:end:step] constructs a new Sequence that is a subsequence as described by the way this function is called. This function will automatically fill in the name, sequence, start, stop, step, original, and type of the subsequence. It also tries to fill in the annotations, but annotations are handled pretty poorly right now, so it's probably best not to worry about those, but it will work if you really want to.'''

    try:
      start, stop, step = key.indices(len(self.seq))
    except AttributeError: start, stop, step = key, key+1, 1
    order = 1 if step >= 0 else -1
    sorder = 1 if self.step > 0 else -1
    dstart, dend = (self.step * start, self.step * (stop-1))[::order]
    pstart, pend = (self.start, self.end)[::sorder]
    newstart, newend = pstart + dstart, pstart + dend
    sub = Sequence("subsequence(%s,%d,%d,%d)" % (self.name, start, stop, step),
                   self.seq[start: stop: step], qual = self.qual and self.qual[start: stop: step],
                   original = self.original, type = self.type, start = newstart, end = newend, step = step * self.step)
    return sub

  '''Some other things you can do with a Sequence object:
* len(sequence)    => gives the length of the sequence.
* sequence.upper() => technically does nothing, but returns the sequence with a capitalized version of the sequence (which is done on instantiation, now).
* for character in sequence: => allows you to loop over each character in the sequence.
* dictionary[sequence] => allows sequences to be used as keys for dictionaries and allows you to have sequences in sets. This relies on the test seqA == seqB, described next.
* seqA == seqB     => compare two sequences. The sequences are the same if they have the same sequence AND name. Therefore, two sequences with different names are treated as separate items in a set and separate keys in a dictionary. If you need to match only the sequence, use seqA.seq == seqB.seq.
* print sequence   => print a fasta / fastq (depending on whether there are any quality scores) representation of the sequence. Sequence objects in any other data structure (e.g., list, dictionary) are printed as (e.g., <Sequence 0x000000>). If you want to change that, you can do:
	def __repr__(self): return self.__str__()'''

  def __len__(self):       return len(self.seq)
		
  def upper(self):         return self
	
  def __iter__(self):
    for i in xrange(len(self)): yield self.seq[i:i+1]
    raise StopIteration
		
  def __hash__(self):      return self.seq.__hash__();
		
  def __eq__(self,other):
    try:                   return self.seq == other.seq and self.name == other.name
    except AttributeError: return (self.seq == other)
			
  def __str__(self):
    if self.qual:          return '@%s\n%s\n+\n%s' % (self.name, self.seq, ''.join(chr(ord('A')-1+q) for q in self.qual))
    else:                  return '>%s %s\n%s' % (self.name, self.defline, '\n'.join(chop(self.seq,70)))

def _io_methods():
	def clean_alignment(x):
		i = 0
		for c in x:
			if c != "-": break
			i += 1

		j = 0
		for c in reversed(x):
			if c != "-": break
			j += 1

		return (' '*i + x[(i or None):(-j or None)] + ' '*j,i,len(x)-j)

	def read_fasta(fh):
		name, defline, seq = '', '', ''
		for line in fh:
			line = line.strip()
			if line[0] == '>':
				if name or seq:
					yield Sequence(name, seq, defline=defline)
				seq = ''
				name = line[1:].split()[0]
				defline = line[1+len(name):].strip()
				continue
			seq += line
		if name or seq:
			yield Sequence(name, seq, defline=defline)
		fh.close()
		raise StopIteration

	def read_fastq(fh):
		while 1:
			try:
				line = fh.next().strip()
				if line[0] == '@':
					name = line[1:].split()[0]
					defline = line[1+len(name):].strip()
					seq = f.next().strip()
					fh.next()
					qual = [ord(c)-self.phred for c in fh.next().strip()]
					yield Sequence(name, seq, qual=qual, defline=defline)
			except:
				raise StopIteration
			finally:
				fh.close()

	def read_fastc(fh):
		c = set()
		for s in read_fasta(fh):
			c.add(s)
			if s.seq:
				for r in c: r.seq = s.seq
				yield c
				c = set()
		if c: yield c
		raise StopIteration

	def read_clustalw(fh):
		seqs = {}
		for line in fh:
			st = line.strip()
			if st:
				bits = st.split()
				if len(bits) != 2: continue
				if bits[0] not in seqs:
					seqs[bits[0]] = ''
				seqs[bits[0]] += bits[1]

		for k in seqs:
			seq, start, end = clean_alignment(seqs[k])
			yield Sequence(k, seq, start=start, end=end)
		raise StopIteration

	def write_fasta(fh, s):
		fh.write('>%s %s\n' % (s.name, s.defline) + \
		         '\n'.join(chop(s.seq, 70)) + '\n')

	def write_fastq(fh, s):
		fh.write('@%s %s\n%s\n+\n%s\n' % (s.name, s.defline, s.seq, \
		         ''.join(q+chr('A')-1 for q in s.qual)) + '\n')

	def write_fastc(fh, c):
		c = set(c)
		fh.write('\n'.join('>%s %s' % (s.name, s.defline) for s in c) + \
		         '\n' + '\n'.join(chop(c.pop().seq,70)) + '\n')

	def probe_fasta(fh):
		for line in fh:
			st = line.strip()
			if st:
				fh.close()
				if st[0] == '>':
					return {'seqtype': 'fasta'}
				return False
		fh.close()
		return False

	def probe_fastq(fh):
		for line in fh:
			st = line.strip()
			if st:
				fh.close()
				if st[0] == '@':
					fh.next(); fh.next() # ignore sequence, duplicate name line
					qual = [ord(c) for c in fh.next().strip()]
					phred = 32 if min(qual)<ord('A') else 64
					qual = [q-phred for q in qual]
					return {'seqtype': 'fastq', 'phred': phred, 'qual': qual}
				return False
		fh.close()
		return False

	def probe_fastc(fh):
		return False

	def probe_clustalw(fh):
		if fh.next().strip().split()[0].lower() == 'CLUSTAL':
			return {'seqtype': 'clustalw'}
		return False

	nil = lambda *x: None
	def pop(fh):
		try: fh.next()
		except StopIteration: pass

	return {
		'fasta': {'rhook': nil, 'read':  read_fasta,
		          'whook': nil, 'write': write_fasta,
		                        'probe': probe_fasta},
		'fastq': {'rhook': nil, 'read':  read_fastq,
		          'whook': nil, 'write': write_fastq,
		                        'probe': probe_fastq},
		'fastc': {'rhook': nil, 'read':  read_fastc,
		          'whook': nil, 'write': write_fastc,
		                        'probe': probe_fastc},
		'clustalw': {'rhook': pop, 'read': read_clustalw,
		             'whook': nil, 'write': nil, 'probe': probe_clustalw}
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
		self.phred  = 0
		self.seqtype = None
		self.haswritten = False

		suffixes = {'fsa': 'fasta', 'fa':  'fasta',
		            'fs': 'fasta',  'fas': 'fasta',
								'fna': 'fasta', 'gff3': 'gff',
								'fastc': 'fastc', 'fasta': 'fasta',
								'gff': 'gff', 'clustalw': 'clustalw',
		            'clustal': 'clustalw',
		            'aln': 'clustalw', 'fastq': 'fastq'}

		p = name.rfind('.')
		if p > -1:
			ext = name[p+1:]
			if ext in suffixes:
				try:
					self.format(self.seqtype)
					return
				except ValueError:
					pass
		try:
			for method in methods:
				try:
					self.format(method)
					return
				except ValueError:
					pass
		except IOError:
			raise

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

	def read(self, n=None):
		'''read(n=None)
If n is provided, the next (up to) n entries are parsed and returned. Otherwise, all remaining entries are parsed and returned.'''
		if count is None: return [s for s in self]
		return [s for s, i in zip(iter(self), xrange(int(n)))]

	def __iter__(self):
		self.method['rhook'](self.handle)
		return self.method['read'](self.handle)

	def next(self):
		'''next()
Reads a single entry in the file and returns it.'''
		try: return self.read(1)[0]
		except (StopIteration, ValueError):
			raise StopIteration

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
		raise AttributeError, "Sequence reader cannot write."

class Writer(IOBase):
	'''class Writer(IOBase)
A class that wraps IOBase and restricts the ability to read.'''
	def __init__(self, filename, mode='w'):
		IOBase.__init__(self, filename, mode)
			
	def read(self, *args):
		raise AttributeError, "Sequence writer cannot read."

def open(filename, mode='r'):
	'''open(filename, mode='r')
Open a file for parsing or creation. Returns either a Reader or Writer object, depending on the open mode.'''
	if mode == 'r': return Reader(filename)
	if mode == 'w': return Writer(filename)
	if mode == 'a': return Writer(filename, mode='a')
  
def dict(filename):
  '''dict(filename)
Transforms the output of open into a dictionary with the key as the sequence name and the value the corresponding Sequnece.'''
  return __builtin__.dict((s.name,s) for s in open(filename, 'r'))
  
def annotation(seq,source,type,**kwargs):
  '''annotation(sequence, source, type, **kwargs)
Creates an Annotation object for the given sequence from a source (e.g., "phytozome7.0") of a particular type (e.g., "gene").'''
  try: sname = source.name
  except AttributeError: sname = source
  
  start,end = min(seq.start,seq.end),max(seq.start,seq.end)
  strand = '+' if seq.step == 1 else '-'
  attrs = ';'.join('%s=%s'%(key,str(kwargs[key])) for key in kwargs)
  
  return Annotation(seq.original.name,sname,type,start,end,'.',strand,start%3,attrs)
