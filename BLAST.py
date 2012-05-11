#!/usr/bin/python
import biotools.sequence as sequ
import subprocess
import os, sys

def run(db,sfile,**kwargs):
	'''BLAST(database, query, **params)
This function takes a database and a query and runs the appropriate type of
BLAST on them. The database can be an existing BLAST database or a fasta/fastq
file. If it is a sequence file, this function will look in the places where BLAST
would look for an existing database created from that file and use that instead.
If there is no such database, this function will make one for you and then use the
newly created database with BLAST. 

Optional named arguments can currently only be evalue.'''

	sep = os.sep
	cmds = {'p':{'p':'blastp','n':'tblastn'},'n':{'n':'blastn','p':'blastx'}}

	seq = sequ.open(sfile, 'r').next()
	qtype = 'p' if set(seq.seq)-set('ATCGNYR') else 'n'

	rcloc = ''
	for loc in (".:~:"+ (os.getenv("NCBI") or "")).split(':'):
		if loc and loc[-1] == sep: loc += sep
		try:
			for line in (l.strip() for l in open(loc+'.ncbirc','r')):
				pos = line.find('=')
				if pos >= 0 and line[:pos].strip() == "BLASTDB":
					rcloc = line[pos+1:].strip()
		except IOError: pass

	dbtype = None
	bdbenv = os.getenv("BLASTDB")
	dblocations = (":." + ((':' + bdbenv) if bdbenv else '') + ((':' +rcloc) if rcloc else '')).split(':')
	for loc in dblocations:
		if loc and loc[-1] != sep: loc += sep
		try:
			open(loc + db+'.pin','r')
			dbtype = 'p'
			break
		except IOError:
			try:
				open(loc + db+'.nin','r')
				dbtype = 'n'
				break
			except IOError: pass
				
	if not dbtype:
		pos = db.rfind(".")
		if pos >= 0 and db[pos+1:] in ["txt","fasta","fa","fas"]:
			for seq in sequ.open(db, 'r'):
				if set(seq.seq) - set('ATCGNRY'): dbtype = 'p'
				else: dbtype = 'n'
				break
			if not dbtype: raise IOError, "Database not found: " + db
			ndb = None
			dbdir = sep.join(db.split(sep)[:-1]) or '.'
			for file in os.listdir(dbdir):
				dpos = file.rfind('.')
				if dpos >= 0 and file[dpos+1:] == dbtype + 'in':
					fh = open(dbdir + sep + file,'r')
					c = ord(fh.read(12)[-1])
					fname = fh.read(c)
					if fname in set([db,'"%s"'%db]):
						ndb = dbdir + sep + file[:dpos]
						break
			if not ndb:
				ndb = '_'.join(db[:pos].split())
				try:
					ignore = open('/dev/null','w')
					mbdb	 = 'makeblastdb'
				except:
					ignore = open('nul', 'w')
					mbdb	 = 'makeblastdb.exe'
				dt = 'nucl' if dbtype == 'n' else 'prot'
				subprocess.call([mbdb,"-in",'"%s"' % db,"-out",ndb,"-dbtype",dt],stdout=ignore)
				db = ndb
			else: db = ndb
		else: raise IOError, "Database not found: " + db
	allowed = set(["evalue", "gapopen", "gapextend", "num_threads"]) & set( kwargs.keys() )
	cmd = cmds[qtype][dbtype]
	if sys.platform in ('win32', 'cygwin'): cmd += '.exe'
	return Result(subprocess.check_output([cmds[qtype][dbtype],"-db",db,"-query",sfile] + \
				 [arg for pair in [["-"+k,str(kwargs[k])] for k in allowed] for arg in pair]).split('\n'))
		
class Result(object):
	'''Class BLASTResult
A class which take the raw output from BLAST and generates dictionaries from
the data from BLAST. This data includes the alignment, percetn identity, gaps,
e-value, score, length of subject, length of query, and start and stop positions
for both sequences.

This class should be used in a for loop like so:
for res in BLASTResult(file_or_data): pass

The class instance has a single other property, headers, which are the lines in
BLAST results before the BLAST hits (e.g., citation info, etc.).'''

	def __init__(self,file):
		self.file = file
		self.headers = []
		
	def __iter__(self):
		name, qname, counts, newseq = '', '', 0, False
		sequence, naming = [], False
		try: handle = open(self.file,'r')
		except TypeError: handle = iter(self.file)
		while 1:
			try: line = (lambda x: x[:-1] if x and x[-1] == '\n' else x)(handle.next())
			except StopIteration:
				if name:
					obj = _parse(sequence,name,qname)
					if obj: yield obj
				raise StopIteration
			if not line:
				counts += 1
				if counts == 2: newseq = True
				continue
			else: counts = 0
			if line[:6] == 'Query=':
				qname = line[6:]
				sequence = []
				newseq = False
				name = ''
				naming = True
			elif naming:
				if line[:7] == 'Length=':
					naming = False
					qname = qname.strip().split()[0]
					sequence.append(line)
				else:
					qname += line[:-1]
			elif line[0] == '>':
				newseq = False
				if name:
					obj = _parse(sequence,name,qname)
					if obj: yield obj
				else:
					self.headers = sequence[:]
				name = line[1:].split()[0]
				sequence = []
			elif newseq:
				newseq = False
				if name:
					obj = _parse(sequence,name,qname)
					if obj: yield obj
				sequence = [line]
			else: sequence.append(line)
			
	def __getattribute__(self,name):
		if name == 'headers':
			if not object.__getattribute__(self,'headers'):
				for a in self: pass
		return object.__getattribute__(self,name)
	
def _parse(sequence,name,qname):
	'''parse(data, subjectName, queryName)
Parses a single BLAST result and returns a dictionary of the information.'''

	obj = {'query':{'name':qname,'start':None,'end':None,'sequence':''},'subject':{'name':name,'start':None,'end':None,'sequence':''}}
	subheaders = True
	for line in sequence:
		if line[:5] == 'Query': curr = 'query'
		elif line[:5] == 'Sbjct': curr = 'subject'
		elif subheaders:
			for name,value in ((lambda x,y='': (x.split('(')[0].strip().lower(),y.strip()))(*i.split('=')) for i in line.split(',') if line if len(i.split('=')) == 2):
				obj[name] = value
			continue
		else: continue # garbage: empty or alignment lines
		subheaders = False
		bits = line.split()
		obj[curr]['start'] = obj[curr]['start'] or bits[1] # it's a string for now ('0' is not falsy)
		obj[curr]['end']	 = bits[3]
		obj[curr]['sequence'] += ''.join(bits[2].split('-'))
	for i in ('query','subject'):
		for j in ('start','end'):
			if obj[i][j]: obj[i][j] = int(obj[i][j])
			else: return None
		obj[i]['length'] = abs(obj[i]['end']-obj[i]['start'])+1
	return obj
		
if __name__ == "__main__":
	import sys
	for result in BLAST(sys.argv[1],sys.argv[2]):
		print result or "BUGGY BUGGY BOO"
