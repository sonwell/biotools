'''
This module is home to the IOManager class, which manages the various input
and output formats (specifically, FASTA, FASTQ, CLUSTAL alignments, and GFF
files, currently).
'''


class IOManager(object):
    '''
    A class used by the `IOBase` class to manage the various input and output
    methods for the different file types. Additional file types can be added 
    to the manager by using 

    ```python
    manager[format] = methods
    ```

    From the above example, `methods` is a dictionary with keys `rhook`, 
    `read`, `whook`, `write`, and `probe`. Each of the values must be callable
    object:
    * `rhook` => takes a file handle opened for reading; called before reading
        of the file has begun,
    * `whook` => takes a file handle opened for writing; called before writing
        to the file has begun,
    * `read`  => takes a file handle opened for reading; should be a generator
        that yields entries,
    * `write` => takes a file handle opened for writing and a single entry;
        writes the entry to the file,
    * `probe` => takes a file handle opened for reading; returns a dictionary
        of attributes to be applied to the `IOBase` instance.

    This class behaves similarly to a dictionary, except that the get method
    will default to the default method (which does nothing) if no truthy
    second parameter is passed.
    '''

    def __init__(self, methods=None):
        '''
        Instantiates an `IOManager` with methods, where the keys of methods are
        the formats and the values are dictionaries with `rhook`, `whook`, 
        `read`, `write`, and `probe` callables.
        '''
        self.methods = methods or {}
        self.protected = set(self.methods.keys())

        nil = lambda *x: iter([])
        self.default = {
            'rhook': nil,
            'read':  nil,
            'whook': nil,
            'write': nil,
            'probe': nil
        }

    def __contains__(self, key):
        return key in self.methods

    def __getitem__(self, key):
        return self.methods.get(key, self.default)

    def get(self, key, default=None):
        '''
        Try to get a set of methods via format (e.g., 'fasta') or fall-back
        to the default methods (which do nothing).
        '''
        try:
            return self.methods[key]
        except KeyError:
            return default or self.default

    def __setitem__(self, key, value):
        if key not in self.protected:
            self.methods[key] = value
        return value

    def __iter__(self):
        return iter(self.methods)
