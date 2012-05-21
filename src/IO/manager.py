class IOManager(object):
	'''class IOManager(object)
IOManager is a class used by the IOBase classes in sequence.py and annotation.py to manage the various input and output methods for the different file types. Additional file types can be added to the manager by using
	manager[format] = methods.

From the above example, methods is a dictionary with keys rhook, read, whook, write, and probe. Each of the values must be callable object:
* rhook => takes a file handle, opened for reading; called before reading of the file has begun,
* whook => takes a file handle, opened for writing; called before writing to the file has begun,
* read  => takes a file handle, opened for reading; should be a generator that yields entries,
* write => takes a file handle, opened for writing and a single entry; writes the entry to the file,
* probe => takes a file handle, opened for reading; returns a dictionary of attributes to be applied to the IOBase instance.

This class behaves similarly to a dicitonary, except that the get method will default to the default method (which does nothing) if no truthy second parameter is passed.'''
	def __init__(self, methods=None):
		'''IOManager(methods=None)
Instantiates an IOManager with methods, where the keys of methods are the formats and the values are dictionaries with rhook, whook, read, write, and probe callables.'''
		self.methods   = methods or {}
		self.protected = set(self.methods.keys())

		nil = lambda *x: None
		self.default = {'rhook': nil, 'read':  nil,
		                'whook': nil, 'write': nil,
										'probe': nil}

	def __contains__(self, key):
		return key in self.methods

	def __getitem__(self, key):
		return self.methods.get(key, self.default)

	def get(self, key, default=None):
		try: return self.methods[key]
		except: return default or self.default

	def __setitem__(self, key, value):
		if key not in self.protected:
			self.methods[key] = value
		return value

	def __iter__(self):
		return iter(self.methods)
