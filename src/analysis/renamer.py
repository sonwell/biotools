import os, sys
import biotools.IO as io
import options

def rename(direc, db):
	'''rename(directory, database)
This isn't really for bioinformatics, this is more for the pipeline, to rename the files generated by cluster.py with a little human interaction.'''

	sep = os.sep
	nt_dir = direc + 'nt' + sep
	aa_dir = direc + 'aa' + sep

	names = []
	files = [f for f in os.listdir(nt_dir) if f[-6:] == ".fasta"]
	seqdb = dict((s.name.split("|")[1],s) for s in io.open(db, 'r'))
	for f in files:
		seq = io.open(nt_dir+f, 'r').next()
		ids = seq.defline.split(', ')
		print "File\033[33;1m", f, "\033[0mis described by the following effectors:"
		try:
			for id in ids:
				seqdb[id]
				print "*", seqdb[id].defline.split('[')[0]
		except: 
			print "* (none)"
			continue
		pre = raw_input("\033[33;1mWhat should we call this file (or hit enter to skip)? \033[0m")
		fpre = f[:f.find('.')]

		if pre != "":
			count = 0
			while True:
				rpre = pre + ((" (%d)"%count) if count > 0 else "")
				try:
					fh = open(nt_dir+rpre+".fasta", 'r')
					fh.close()
					count += 1
					continue
				except:# IOError:
					try: 
						print "Renaming "+fpre+".fasta to "+rpre+".fasta"
						os.rename(nt_dir+fpre+".fasta", nt_dir+rpre+".fasta")
						os.rename(aa_dir+fpre+".fasta", aa_dir+rpre+".fasta")
						print "Renaming "+fpre+".clustalw to "+rpre+".clustalw"
						os.rename(nt_dir+fpre+".clustalw", nt_dir+rpre+".clustalw")
						os.rename(aa_dir+fpre+".clustalw", aa_dir+rpre+".clustalw")
						names.append(rpre)
					except: pass
					break
	return names