#`biotools`
  **Needs documentation**

##`biotools.BLAST`
  **Needs documentation**

###`biotools.BLAST.Result(self, file)`
  
  A class which take the raw output from BLAST and generates dictionaries
  from the data from BLAST. This data includes the alignment, percent
  identity, gaps, e-value, score, length of subject, length of query, and
  start and stop positions for both sequences.
  
  ```python
  This class should be used in a for loop like so:
      for res in Result(file_or_data): pass
  ```
  
  The class instance has a single other property, `headers`, which are the
  lines in BLAST results before the BLAST hits (e.g., citation info, etc.).

###`biotools.BLAST.run(db, sfile, mega_blast=False, **kwargs)`
  
  Takes a database and a query and runs the appropriate type of BLAST on 
  them. The database can be an existing BLAST database or a fasta/fastq 
  file. If it is a sequence file, this function will look in the places 
  where BLAST would look for an existing database created from that file and 
  use that instead. If there is no such database, this function will make 
  one for you and then use the newly created database with BLAST.
  
  Optional named arguments can currently only be `evalue`, `num_threads`, 
  `gapopen`, or `gapextend`. The correspond to the BLAST options of the same 
  name.

##`biotools.IO`
  **Needs documentation**

###`biotools.IO.IOBase(self, name, mode)`
  
  Generic IO class for sequence files.

####`biotools.IO.IOBase.close(self)`
  
  Close the file handle.

####`biotools.IO.IOBase.format(self, fmt)`
  
  Forces a file to be parsed as a particular format. By default, the
  values for fmt can be any recognized format.

###`biotools.IO.Reader(self, filename, mode='r')`
  
  A class that wraps IOBase and restricts the ability to write.

####`biotools.IO.Reader.next(self)`
  
  Reads a single entry in the file and returns it.

####`biotools.IO.Reader.read(self, n=None)`
  
  If n is provided, the next (up to) n entries are parsed and returned.
  Otherwise, all remaining entries are parsed and returned.

###`biotools.IO.Writer(self, filename, mode='w')`
  
  A class that wraps IOBase and restricts the ability to read.

####`biotools.IO.Writer.write(self, sequence)`
  
  Writes sequence as the correct format to the file.

###`biotools.IO.clustal`
  **Needs documentation**

####`biotools.IO.clustal.clean_alignment(x)`
  **Needs documentation**

####`biotools.IO.clustal.probe(fh)`
  **Needs documentation**

####`biotools.IO.clustal.read(fh)`
  **Needs documentation**

####`biotools.IO.clustal.rhook(fh)`
  **Needs documentation**

###`biotools.IO.fasta`
  **Needs documentation**

####`biotools.IO.fasta.probe(fh)`
  **Needs documentation**

####`biotools.IO.fasta.read(fh)`
  **Needs documentation**

####`biotools.IO.fasta.write(fh, s)`
  **Needs documentation**

###`biotools.IO.fastq`
  **Needs documentation**

####`biotools.IO.fastq.probe(fh)`
  **Needs documentation**

####`biotools.IO.fastq.read(fh)`
  **Needs documentation**

####`biotools.IO.fastq.write(fh, s)`
  **Needs documentation**

###`biotools.IO.get_methods()`
  **Needs documentation**

###`biotools.IO.gff`
  **Needs documentation**

####`biotools.IO.gff.probe(fh)`
  **Needs documentation**

####`biotools.IO.gff.read(fh)`
  **Needs documentation**

####`biotools.IO.gff.whook(fh)`
  **Needs documentation**

####`biotools.IO.gff.write(fh, a)`
  **Needs documentation**

###`biotools.IO.manager`
  **Needs documentation**

####`biotools.IO.manager.IOManager(self, methods=None)`
  
  IOManager is a class used by the IOBase classes in sequence.py and
  annotation.py to manage the various input and output methods for the
  different file types. Additional file types can be added to the manager
  by using manager[format] = methods.
  
  From the above example, methods is a dictionary with keys rhook, read,
  whook, write, and probe. Each of the values must be callable object:
  * rhook => takes a file handle, opened for reading; called before reading
      of the file has begun,
  * whook => takes a file handle, opened for writing; called before writing
      to the file has begun,
  * read  => takes a file handle, opened for reading; should be a generator
      that yields entries,
  * write => takes a file handle, opened for writing and a single entry;
      writes the entry to the file,
  * probe => takes a file handle, opened for reading; returns a dictionary
      of attributes to be applied to the IOBase instance.
  
  This class behaves similarly to a dicitonary, except that the get method
  will default to the default method (which does nothing) if no truthy
  second parameter is passed.

#####`biotools.IO.manager.IOManager.get(self, key, default=None)`
  **Needs documentation**

###`biotools.IO.open(filename, mode='r')`
  
  Open a file for parsing or creation. Returns either a Reader or Writer
  object, depending on the open mode.

##`biotools.align`
  **Needs documentation**

###`biotools.align.OptimalCTether(reference, translation, gp=1, c=10)`
  
  This function will take two sequences: a `reference` sequence any other
  protein sequence (`translation`; usually, this is an open reading frame
  that has been translated). Needleman-Wunsch alignment will be performed
  and the substring of translation with the highest identity that begins
  with a start codon [default: `['ATG']`] is reported.
  
  This function returns a dictionary of relevent information from the
  alignment; specifically, the alignments itself [keys: `query`, `subject`],
  the score [key: `score`], the length of the alignment [key: `length`], the
  length of the substring of translation used [key: `sublength`], the number
  of identities [key: `identities`], the theoretical perfect score for that
  alignment [key: `perfect`], and the number of gaps [key: `gaps`].

##`biotools.analysis`
  **Needs documentation**

###`biotools.analysis.cluster`
  **Needs documentation**

####`biotools.analysis.cluster._run_clustal(q, clusters, direc, names)`
  **Needs documentation**

####`biotools.analysis.cluster.run(direc, inputs)`
  
  Takes a collection of files generated by gene prediction, creates clusters
  based off of the genes that have homology to those predicted genes, and
  creates new fasta files in the clusters sub directory under the given
  directory and separated according to whether they are nucleotide or amino
  acid sequnces. These new fasta files are then used to create clustalw
  alignments of the genes if more than 1 sequence exists in the fasta file.

###`biotools.analysis.loaddata`
  
  This is a pretty simple JSON-like parser. Specifically, it can load Python-like
  object, list, and other literals, i.e., the sort of stuff you'd get it you
  dumped the the string representation of some data into a file.
  
  The real difference is that you must specify a variable name, e.g.:
  
      ```
      my_stuff = { ... }
      ```
  
  These variable names don't need to be on a newline or anything like that, you
  should be able to omit any and all whitespace. The result of a successful 
  parse is a dictionary:
  
      ```
      {'my_stuff': { ... }}
      ```
  
  This function really only works for `None`, `True`, `False`, numbers, strings, 
  dictionaries, and lists.

####`biotools.analysis.loaddata.parse(ipt)`
  **Needs documentation**

###`biotools.analysis.options`
  **Needs documentation**

####`biotools.analysis.options.debug(msg)`
  **Needs documentation**

####`biotools.analysis.options.help()`
  
  help()
  Prints the usage.

####`biotools.analysis.options.parse(pargs)`
  
  Parses pargs and sets the various variables to be accessible to other
  modules.
  
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
  * args

###`biotools.analysis.plot`
  **Needs documentation**

####`biotools.analysis.plot.axes(bottom, side, bound, fig, **kwargs)`
  **Needs documentation**

####`biotools.analysis.plot.draw(x, y, ax, color, **kwargs)`
  **Needs documentation**

####`biotools.analysis.plot.models(starts, ends, counts, bound, ax, **kwargs)`
  **Needs documentation**

####`biotools.analysis.plot.plot(plotdata, directory, bottom=True, side=True, legend=True, save=True, filename='untitled.pdf', upperbound=0.05, factor=21, fig=<matplotlib.figure.Figure object at 0x102237790>, **kwargs)`
  **Needs documentation**

####`biotools.analysis.plot.report(ntvar, aavar, lnt, laa)`
  **Needs documentation**

####`biotools.analysis.plot.smoothed(unsmoothed, factor)`
  **Needs documentation**

###`biotools.analysis.predict`
  **Needs documentation**

####`biotools.analysis.predict.GeneFromBLAST(db, sequences, pref, names)`
  
  BLASTs database against sequences, and for those results that pass the
  length and percent identity requirements, attempt to locate the full gene
  that corresponds to that BLAST hit. Genes that are found are saved in the
  subdirectory sequences under the given directory, divided depending on
  whether the sequnece is amino acid or nucleotide.

####`biotools.analysis.predict.ORFGenerator(sequ)`
  
  Scans both strands of the given sequence and yields the longest subsequence
  that starts with a start codon and contains no stop codon other than the
  final codon.

####`biotools.analysis.predict.ThreadQueue`
  **Needs documentation**

####`biotools.analysis.predict._genepredict_in_range(seq, start, end, frame)`
  **Needs documentation**

####`biotools.analysis.predict._genepredict_target(qin, qout, orfs, subj)`
  **Needs documentation**

####`biotools.analysis.predict.run(subject, query, prefix, names)`
  **Needs documentation**

###`biotools.analysis.renamer`
  **Needs documentation**

####`biotools.analysis.renamer.rename(direc, db, files)`
  
  This isn't really for bioinformatics, this is more for the pipeline, to
  rename the files generated by cluster.py with a little human interaction.

###`biotools.analysis.report`
  **Needs documentation**

####`biotools.analysis.report.plot(plotdata, directory, bottom=True, side=True, legend=True, save=True, filename='untitled.pdf', upperbound=0.05, factor=21, fig=<matplotlib.figure.Figure object at 0x1022a0b90>, **kwargs)`
  **Needs documentation**

####`biotools.analysis.report.report(plotdata, **kwargs)`
  **Needs documentation**

###`biotools.analysis.run`
  **Needs documentation**

####`biotools.analysis.run._run_genepredict(q, infile, names)`
  **Needs documentation**

####`biotools.analysis.run.run(infile, strains)`
  
  Run several instances of genepredict.run at once.

###`biotools.analysis.variance`
  **Needs documentation**

####`biotools.analysis.variance.SaySNPs(input)`
  
  Takes a clustalw alignment and will return a dictionary of data
  relevent to plotting the sequence variance for the sequences in the
  given clustalw alignment. These data are:
  * `var`: the measure of sequence variation,
  * `starts`: the starting positions for each gene model in amino acids,
  * `ends`: the ending positions for each gene model in amino acids, and
  * `count`: the number of sequences with a particular gene model.
  The values given in `starts`, `ends`, and `counts` are sorted to that the 
  nth element in starts corresponds to the nth value in ends and the nth 
  value in counts.

####`biotools.analysis.variance.var(strain, fmt)`
  
  Returns plot data and metadata for plotting later on in the pipeline.

##`biotools.annotation`
  **Needs documentation**

###`biotools.annotation.Annotation(self, ref, src, type, start, end, score, strand, phase, attr, name_token='ID', gff_token='=')`
  
  An object to help with reading and writing GFF files.

###`biotools.annotation._parseAttrs(attr, token='=')`
  
  Creates a dictionary from the atrributes (9th column) of a gff file. By
  default, token is `=`, which is the separator used in gff version 3.
  
  In other words, `attr "a=b;c=d;"` and token `=` will yield the dictionary
  `{'a':'b','c':'d'}`. The other separator (`;`) cannot be changed.
  
  This function is not to be called on its own.

##`biotools.blosum62`
  
  This is a pretty uninteresting file, the BLOSUM62 matrix I think I ripped from
  wikipedia, and just hacked it to pieces to make it do what I want. Lookup
  looks something like `blosum62['A','C']` and would give the value corresponding
  to the replacement of alanine with cystine.
  
  Basically, you don't need to be using this file, it is just to help in
  align.py.

###`biotools.blosum62.blosum62(self)`
  **Needs documentation**

##`biotools.clustal`
  **Needs documentation**

###`biotools.clustal.run(infile, outfile, **kwargs)`
  **Needs documentation**

##`biotools.complement`
  **Needs documentation**

###`biotools.complement.complement(s)`
  
  Creates the complement of a sequence, which can then be reversed by using
  `seq[::-1]`, if it needs to be reversed. This function accepts either
  `Sequence`s or strings.

##`biotools.sequence`
  **Needs documentation**

###`biotools.sequence.Sequence(self, name, seq, **kwargs)`
  
  A wrapper class for sequences.

####`biotools.sequence.Sequence.upper(self)`
  **Needs documentation**

###`biotools.sequence.annotation(seq, source, type, **kwargs)`
  
  Creates an `Annotation` object for the given sequence from a source
  (e.g., "phytozome7.0") of a particular type (e.g., "gene").

###`biotools.sequence.chop(seq, length=70)`
  
  Yields a chunk of a sequence of no more than `length` characters,
  it is meant to be used to print fasta files.

##`biotools.translate`
  **Needs documentation**

###`biotools.translate.translate(sequence)`
  
  Translate a nucleotide using the standard genetic code. The sequence
  parameter can be either a string or a `Sequence` object. Stop codons are
  denoted with an asterisk (*).

