'''
A module for reading and writing to sequence and annotation files. Currently
supported file types are: FASTA, FASTQ, CLUSTAL alignments, and GFF3 files.
'''

from biotools.IO import fasta, fastq, gff, clustal
from biotools.IO.manager import IOManager
try:
    import __builtin__
except ImportError:
    import builtins as __builtin__


def get_methods():
    methods = {}
    nil = lambda *x: iter([])
    for module in [fasta, fastq, gff, clustal]:
        modname = module.__name__.split('.')[-1]
        methods[modname] = {}
        for method in ['read', 'write', 'rhook', 'whook', 'probe']:
            methods[modname][method] = module.__dict__.get(method, nil)
    return methods


class IOBase(object):
    '''
    Generic IO class for sequence files.
    '''
    methods = IOManager(get_methods())

    def __init__(self, name, mode):
        '''
        Opens file name with mode mode. This function will attempt to guess at
        the filetype by 1. looking at the file extension and failing that,
        will 2. read the first few lines to determine the file type.

        Recoginized file extensions include fa, fsa, fas, fasta, fastq,
        clustalw, clustal, aln.
        '''

        self.file = name
        self.handle = __builtin__.open(name, mode)
        self.method = self.methods.default
        self.type = None

        self.suffixes = {
            'fsa': 'fasta',
            'fa':  'fasta',
            'fs': 'fasta',
            'fas': 'fasta',
            'fna': 'fasta',
            'fasta': 'fasta',
            'clustalw': 'clustal',
            'clustal': 'clustal',
            'aln': 'clustal',
            'fastq': 'fastq',
            'gff': 'gff',
            'gff3': 'gff'
        }

        p = name.rfind('.')
        if p > -1:
            ext = name[p + 1:]
            if ext in self.suffixes:
                try:
                    self.format(self.suffixes[ext])
                    return
                except ValueError:
                    pass
        try:
            for method in IOBase.methods:
                try:
                    self.format(method)
                    return
                except ValueError:
                    pass
        except IOError:
            raise

    def format(self, fmt):
        '''
        Forces a file to be parsed as a particular format. By default, the
        values for fmt can be any recognized format.
        '''
        if fmt in self.methods:
            method = self.methods[fmt]
            ret = method['probe'](__builtin__.open(self.file, 'r'))
            if ret:
                self.method = method
                for key in ret:
                    object.__setattr__(self, key, ret[key])
                return self
            else:
                raise ValueError("File cannot be parsed as type %s." % fmt)
        self.method = self.methods.default
        return self

    def close(self):
        '''
        Close the file handle.
        '''
        self.handle.close()


class Reader(IOBase):
    '''
    A class that wraps IOBase and restricts the ability to write.
    '''
    def __init__(self, filename, mode='r'):
        IOBase.__init__(self, filename, mode)
        self.method['rhook'](self.handle)
        self.iter = self.method['read'](self.handle)

    def read(self, n=None):
        '''
        If `n` is provided, the next (up to) `n` entries are parsed and 
        returned. Otherwise, all remaining entries are parsed and returned.
        '''
        if n is None:
            return [s for s in self]
        return [self.iter.next() for i in xrange(int(n))]

    def __iter__(self):
        return self.iter

    def next(self):
        '''
        Reads a single entry in the file and returns it.
        '''
        try:
            return self.read(1)[0]
        except (StopIteration, ValueError, IndexError):
            raise StopIteration()


class Writer(IOBase):
    '''
    A class that wraps IOBase and restricts the ability to read.
    '''

    def __init__(self, filename, mode='w'):
        IOBase.__init__(self, filename, mode)
        self.haswritten = False

    def write(self, sequence):
        '''
        Writes sequence as the correct format to the file.
        '''
        if not self.haswritten:
            self.method['whook'](self.handle)
            self.haswritten = True
        self.method['write'](self.handle, sequence)


def open(filename, mode='r'):
    '''
    Open a file for parsing or creation. Returns either a Reader or Writer
    object, depending on the open mode.
    '''
    if mode == 'r':
        return Reader(filename)
    if mode == 'w':
        return Writer(filename, mode='w')
    if mode == 'a':
        return Writer(filename, mode='a')
