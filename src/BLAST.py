#!/usr/bin/python
import biotools.IO as io
import subprocess
import os, sys

def run(db, sfile, mega_blast=False, **kwargs):
	'''BLAST(database, query, **params)
This function takes a database and a query and runs the appropriate type of BLAST on them. The database can be an existing BLAST database or a fasta/fastq file. If it is a sequence file, this function will look in the places where BLAST would look for an existing database created from that file and use that instead. If there is no such database, this function will make one for you and then use the newly created database with BLAST. 

Optional named arguments can currently only be evalue or num_threads.'''

	sep = os.sep
	cmds = {'p':{'p':'blastp','n':'tblastn'},'n':{'n':'blastn','p':'blastx'}}

	seq = io.open(sfile, 'r').next()
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
			for seq in io.open(db, 'r'):
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
	pn = ["-db", "-query"]
	if mega_blast:
		cmd = "megablast"
		pn = ["-d", "-i"]
		allowed = ["e", "a"]

	if sys.platform in ('win32', 'cygwin'): cmd += '.exe'
	return Result(subprocess.check_output([cmd,pn[0],db,pn[1],sfile] + \
				 [arg for pair in [["-"+k,str(kwargs[k])] for k in allowed] for arg in pair]).split('\n'))
		
class Result(object):
	'''Class BLASTResult
A class which take the raw output from BLAST and generates dictionaries from the data from BLAST. This data includes the alignment, percetn identity, gaps, e-value, score, length of subject, length of query, and start and stop positions for both sequences.

This class should be used in a for loop like so: 
    for res in Result(file_or_data): pass

The class instance has a single other property, headers, which are the lines in BLAST results before the BLAST hits (e.g., citation info, etc.).'''

	def __init__(self,file):
		self.file = file
		self.headers = []
		
	def __iter__(self):
		try: ipt = open(self.file, 'r')
		except: ipt = self.file
		mode = 0
		headers = []
		curr = None
		length = 0
		
		def sh(sn, qn, l):
			return {
				'subject': {
					'name':  sn,
					'start': None,
					'end':   None,
					'sequence': ''
				},
				'query': {
					'name':  qn.split()[0],
					'start': None,
					'end':   None,
					'sequence': ''
				},
				'length':  l
			}
	
		def ra(sh):
			for res in ('subject', 'query'):
				sh[res]['start']  = int(sh[res]['start'])
				sh[res]['end']    = int(sh[res]['end'])
				sh[res]['length'] = abs(sh[res]['end'] - \
				                    sh[res]['start'] + 1)
			return sh		

		for line in ipt:
			line = line.strip()
			if not line:
				if mode == 4:
					mode = 5
				continue
	
			if mode == 0:
				if line[:6] == 'Query=':
					mode  = 1
					qname = line[6:].strip()
					self.headers = headers
				else:
					headers.append(line)
		
			elif mode == 1:
				if line[0] == '>':	
					mode = 3
					subheaders = sh(line[1:], qname, length)
				elif line[:6] == 'Length' or \
						(line[0] == '(' and line.endswith('letters)')):
					if line[:7] == 'Length=':
						length = int(''.join(line[7:].strip().split(',')))
					else:
						length = int(''.join(line[1:-8].strip().split(',')))
					mode  = 2
				elif line[:6] == 'Query=':
					qname = line[6:].strip()
				else:
					qname += line
		
			elif mode == 2:
				if line[0] == '>':
					mode = 3
					subheaders = sh(line[1:], qname, length)
				elif line[:6] == 'Query=':
					qname = line[6:].strip()
					mode = 1
	
			elif mode == 3:
				if line[:5] == 'Score':
					subheaders['subject']['name'] = subheaders['subject']['name'].split()[0]
					for pairs in (a.strip() for a in line.split(',')):
						l, r = tuple(a.strip() for a in pairs.split('=')[:2])
						subheaders[l.lower()] = r
					mode = 4
				else:
					subheaders['subject']['name'] += line
		
			elif mode == 4:
				for pairs in (a.strip() for a in line.split(',')):
					l, r = tuple(a.strip() for a in pairs.split('=')[:2])
					subheaders[l.lower()] = r
		
			elif mode == 5:
				if   line[:6] == 'Query=':
					mode  = 1
					qname = line[6:].strip().split()[0]
					yield ra(subheaders)
					continue
				elif line[0] == '>':
					yield ra(subheaders)
					subheaders = sh(line[1:], qname, length)
					mode = 3
					continue
				elif line[:5] == 'Score':
					yield ra(subheaders)
					subheaders = sh(subheaders['subject']['name'], qname, length)
					for pairs in (a.strip() for a in line.split(',')):
						l, r = tuple(a.strip() for a in pairs.split('=')[:2])
						subheaders[l.lower()] = r
					mode = 4
					continue
				elif line[:5] == 'Sbjct': curr = 'subject'
				elif line[:5] == 'Query': curr = 'query'
				else: continue
		
				_, start, seq, end = line.split()
				subheaders[curr]['start']     = subheaders[curr]['start'] or start
				subheaders[curr]['end']       = end
				subheaders[curr]['sequence'] += ''.join(seq.split('-')).upper()
	
		yield ra(subheaders)
		raise StopIteration

if __name__ == "__main__":
	import sys
	for result in BLAST(sys.argv[1],sys.argv[2]):
		print result or "BUGGY BUGGY BOO"
