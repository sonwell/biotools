from biotools.IO.manager import IOManager
from biotools.sequence import Sequence
from biotools.annotation import Annotation
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
					return {'type': 'gff', 
					        'version': float(bits[1])}
				return False
		return {'type': 'gff',
		        'version': 3}

	nil = lambda x: None

	def whook_gff(fh):
		fh.write('##gff-version 3\n')

	def read_fasta(fh):
		name, defline, seq = '', '', ''
		for line in fh:
			line = line.strip()
			if not line: continue
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
					return {'type': 'fasta'}
				return False
		fh.close()
		return {'type': 'fasta'}

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
					return {'type': 'fastq', 'phred': phred}
				return False
		return {'type': 'fastq', 'phread': 64}

	def probe_fastc(fh):
		return False

	def probe_clustalw(fh):
		for line in fh:
			st = line.strip()
			if st:
				if st.split()[0].lower() == 'CLUSTAL':
					return {'type': 'clustalw'}
				return False
		return {'type': 'clustalw'}

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
		             'whook': nil, 'write': nil,
		                           'probe': probe_clustalw},
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

		self.file   = name
		self.handle = __builtin__.open(name, mode)
		self.method = self.methods.default
		self.type   = None

		suffixes = {'fsa': 'fasta', 'fa':  'fasta',
		            'fs': 'fasta',  'fas': 'fasta',
		            'fna': 'fasta', 'fastc': 'fastc',
		            'fasta': 'fasta', 'clustalw': 'clustalw',
		            'clustal': 'clustalw', 'aln': 'clustalw',
		            'fastq': 'fastq', 'gff': 'gff',
		            'gff3': 'gff'}

		p = name.rfind('.')
		if p > -1:
			ext = name[p+1:]
			if ext in suffixes:
				try:
					self.format(suffixes[ext])
					return
				except ValueError:
					pass
		try:
			print name
			for method in IOBase.methods:
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
			ret = method['probe'](__builtin__.open(self.file, 'r'))
			if ret:
				self.method = method
				for key in ret:
					object.__setattr__(self, key, ret[key])
				return self
			else: raise ValueError, "File cannot be parsed as type %s." % fmt
		self.method = self.methods.default
		return self

	def close(self):
		'''close()
Close the file handle.'''
		self.handle.close()

class Reader(IOBase):
	'''class Reader(IOBase)
A class that wraps IOBase and restricts the ability to write.'''
	def __init__(self, filename, mode='r'):
		IOBase.__init__(self, filename, mode)

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

class Writer(IOBase):
	'''class Writer(IOBase)
A class that wraps IOBase and restricts the ability to read.'''
	def __init__(self, filename, mode='w'):
		IOBase.__init__(self, filename, mode)
		self.haswritten = False
			
	def write(self, sequence):
		'''write(sequence):
Writes sequence as the correct format to the file.'''
		if not self.haswritten:
			self.method['whook'](self.handle)
			self.haswritten = True
		self.method['write'](self.handle, sequence)

def open(filename, mode='r'):
	'''open(filename, mode='r')
Open a file for parsing or creation. Returns either a Reader or Writer object, depending on the open mode.'''
	if mode == 'r': return Reader(filename)
	if mode == 'w': return Writer(filename)
	if mode == 'a': return Writer(filename, mode='a')
  
