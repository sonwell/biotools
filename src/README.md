`biotools`
==========

A bunch of bioinformatics tools.

`biotools.align`
----------------

`biotools.align.OptimalCTether(reference, translation)`\: uses a hybrid between the Needleman-Wunsch and Smith-Waterman algorithms to find the optimal sub-string of `translation` that begins with a start codon (see also: `biotools.analysis.options`) aligning to `reference` by requiring global alignment at the C terminus and local alignment at the N terminus. 

`OptimalCTether` returns a dictionary of relevent information from the alignment; specifically, the alignments itself [keys: `query`, `subject`], the score [key: `score`], the length of the alignment [key: `length`], the length of the substring of translation used [key: `sublength`], the number of identities [key: `identities`], the theoretical perfect score for that alignment [key: `perfect`], and the number of gaps [key: `gaps`].

`biotools.annotation`
---------------------

`_parseAttrs(attr, token='=')`\: creates a dictionary from the attributes (9th column) of a gff file. By default, token is '=', which is the separator used in gff version 3.

`Annotation(self, ref, src, type, start, end, score, strand, phase, attr, name_token='ID', 'gff_token')`\: constructs an `Annotation` object with the necessary values. The parameters are passed in the same order as the columns from a GFF (version 3) file and the `name_token` and `gff_token` parameters are the defaults for a gff version 3 file from phytozome. Just write (e.g.) `Annotation(*line.split('\t')) #(splitting on tabs)`, and the rest of the work will be done for you. Other sources may require changes to `name_tokens` and `gff_token`.

Instantiating an Annotation will generate for it an id of the form *SEQNAME*\_*TYPE*[*START*\:*END*], where *SEQNAME* is the name of the sequence (column 1) from the GFF file, and *TYPE* is like 'gene' or 'CDS'. If no *SEQNAME* is provided, then 'X' be used in its place, and if no identifier can be found in the attributes, the `Annotation` will generate an identifier for itself in the form of '*unknown #*'.

An instance of this class has the following attributes:
* `strand`
* `type`
* `source`
* `seq`
* `start`
* `end`
* `attr`
* `phase`
* `score`
* `ntoken`
* `id`
* `name`
* `parent`
* `children`

It length (`end-start+1`) can be determined by calling `len(instance)`, they are hashable, they can be tested for equivalence, and casting them to a string renders them as a line of a GFF (version 3 file).

`biotools.BLAST`
----------------

`run(db, sfile, **kwargs)`\: takes a database and a query and runs the appropriate type of BLAST on them. The database can be an existing BLAST database or a fasta/fastq file. If it is a sequence file, this function will look in the places whe    re BLAST would look for an existing database created from that file and use that instead. If there is no such database, this function will make one for you and then use the newly created database with BLAST.

Optional named arguments can currently only be `evalue` or `num_threads`.

`Result(self, file)`\: A class which take the raw output from BLAST and generates dictionaries from the data from BLAST. This data includes the alignment, percent identity, gaps, e-value, score, length of subject, length of query, and start and stop positions for both sequences.

This class should be used in a for loop like so: 
    for res in Result(file_or_data): pass

The class instance has a single other property, `headers`, which are the lines in BLAST results before the BLAST hits (e.g., citation info, etc.).

`biotools.blosum62`
-------------------

`blosum62(self)`\: creates an object that can be subscripted with a tuple of 2 amino acids, and will return the score for substituting the first with the second. This file is mostly just to help in align.py.

`biotools.complement`
---------------------

`complement(s)`\: generates the complement of a sequence. This function accepts either Sequences (see `biotools.sequence`) or strings.

`biotools.IO`
-------------

`IOBase(self, name, mode)`\: opens the file `name` with mode `mode`. An `IOBase` will attempt to guess at the filetpye by 1. looking at the file extension and failing that, will 2. read the first few lines to determine the file type.

Recognized file extension include fa, fsa, fas, fasta, fastc, fastq, clusalw, clustalw, aln, gff, gff3.

`chop(seq, length=70)`\: yields a chunk of a sequence of no more than `length` characters, it is meant to be used to print fasta files.

`IOBase` has a class property, `methods`, which is an `IOManager` (see `biotools.iomanager`).

An `IOBase` instance has the following attributes:
* `file`
* `handle`
* `method`
* `type`
and others that are specific to the particular file type.

`IOBase.format(self, fmt)`\: tries to parse the file as a file of type `fmt`. If it cannot, it will not read the file.

`IOBase.close(self)`\: closes the handle.

`Reader(filename, mode='r')`\: extends `IOBase` to allow the reading of files. Can be used in the iterator protocol (that is, in a `for` loop).

`Reader.read(self, n=None)`\: if `n` is provided, the next (up to) `n` entries are parsed and returned. Otherwise, all remaining entries are parsed and returned.

`Reader.next(self)`\: reads a single entry in the file and returns it.

`Writer(filename, mode='r')`\: extends `IOBase` to allow the writing of files.

`Writer.write(self, sequence)`\: appropriately formats a sequence and writes it to the disk.

`open(filename, mode='r')`\: open a file for parsing or writing; returns either a Reader or Writer, depending on the open mode.

`biotools.IO.manager`
---------------------

`IOManager(self, methods=None)`\: a class used by the `IOBase` class to manage the various input and output methods for the different file types. Additional file types can be added to the manager by using `manager[format] = methods`.

From the above example, methods is a dictionary with keys `rhook`, `read`, `whook`, `write`, and `probe`. Each of the values must be callable object:
* `rhook`\: takes a file handle, opened for reading; called before reading of the file has begun,
* `whook`\: takes a file handle, opened for writing; called before writing to the file has begun,
* `read`\: takes a file handle, opened for reading; should be a generator that yields entries,
* `write`\: takes a file handle, opened for writing and a single entry; writes the entry to the file,
* `probe`\: takes a file handle, opened for reading; returns a dictionary of attributes to be applied to the `IOBase` instance.

This class behaves similarly to a dictionary, except that the `get` method will default to the default method (which does nothing) if no truthy second parameter is passed.

`biotools.sequence`
-------------------

`Sequence(self, name, seq, **kwargs)`\: instantiates a `Sequence` object with sequence `seq`. 

A sequence object has attributes
* `name`
* `seq`
* `qual`
* `type`
* `start`
* `end`
* `step`
* `original`
* `defline`

Some other useful parameters that the Sequence constructor can handle are:
* `qual`\: the quality scores (an array of integers) of the sequence,
* `type`\: the type of the sequence, either 'prot' or 'nucl',
* `start`\: the starting position of the sequence within a supersequence,
* `end`\: the ending position of the sequnece within a supersequence,
* `step`\: the 'step' of the sequence, usually +1 for top-strand sequences, and -1 for bottom-strand sequences, but can handle other values as well,
* `original`\: the original Sequence object from which this one derives,
* `defline`\: the definition line for this sequence from a fasta file.
If one of these are not given, they will default to the most logical value that can be determined from the other values and sequence (e.g., if `end` < `start`, then step is probably -1).

A subsequence can be constructed using the `seq[start:stop:step]` syntax. The length can be calculated using `len`. Instances have a method `upper`, but it doesn't do anything since all sequences are made upper-case upon instantiation. Sequences can be hashed and tested for equivalence and casted to a string (which will render either a fasta entry or a fastq entry depending on the presence of quality scores).

The bases or residues of the sequence can be iterated over by writing, e.g., `for base in seq: pass`.

`annotation(seq, source, type, **kwarg)`\: Creates an `Annotation` (see `biotools.annotation`) object for the given sequence from a source (e.g., "phytozome7.0") of a particular type (e.g., "gene").

`biotools.translate`
--------------------

`translate(x)`\: Translate a nucleotide sequence using the standard genetic code. The sequence parameter can be either a string or a `Sequence` (see `biotools.sequence`) object. Stop codons are denoted with an asterisk (\*).

`biotools.analysis`
===================

Modules used in Bart, Rebecca *et al.* _PNAS Plus_.

`biotools.analysis.cluster`
---------------------------

`_run_clustal(q, clusters, direc)`\: takes a queue of clusters and a dictionary of all of the clusters from each of the inputs from `run` and generates fasta files to be used in a ClustalW alignment.

`run(direc, inputs)`\: generates a set of clusters of sequences based on the presence of these sequences in each of the `inputs`, which are fastc files. Sequence IDs are obtained for each cluster in each input. The resulting clusters are those where for each input, each of the sequences in the cluster are also in the same input cluster or are absent from the input. That is,

    >seq1 input1
    >seq2 input1
    ...

    >seq1 input2
    ...
    >seq2 input2
    ...

will generate clusters `{seq1}, {seq2}` because `input2` has distinct sequences for `seq1` and `seq2` and thus cannot be placed in the same cluster as with `input1`. Instead, we split up `seq1` and `seq2` from `input1` to match `input2`. 


`biotools.analysis.options`
---------------------------

Globals:
* `LENGTH_ERR`\: maximum allowable relative error for hit length (default: 0.2)
* `MAX_EVALUE`\: maximum e-value (default: 1e-30)
* `MIN_IDENTITY`\: minimum identity (default: 0.45)
* `NUM_THREADS`\: number of threads (default: 16)
* `NUM_PROCESSES`\: number of parallel processes (default: 2)
* `PLOTTER`\: plotting module name (default: biotools.analysis.plot)
* `DIRECTORY`\: directory to save results (default: current directory)

* `START_CODONS`\: start codons (default: ATG)
* `STOP_CODONS`\: stop codons (default: TAG, TAA, TGA)
* `args`\: additional remaining arguments 

`parse(pargs)`\: parse pargs and set the global variables.

`biotools.analysis.plot`
------------------------

`plot(snpdata, directory, bottom=True, side=True, legend=True, save=True, filename='untitled.pdf', fig=matplotlib.pyplot.figure(None, facecolor='w', edgecolor='w'), upperbound=0.05, **kwargs)`\: uses matplotlib to plot the sequence variance (`snpdata`) into a pdf in `directory` with under the name `filename`. This code can be used to generate multiple subfigures in the same figure using the booleans `bottom`, `side`, and `lengend` to draw or hide the x-axis, y-axis, and legend, respectively.

`axes(bottom, side, bound, fig, **kwargs)`\: set up the axes for plotting; will hide the x-axis label if `bottom` is `False`, the y-axis if `side` is `False` on with `bound = [min x, max x, min y, max y]` figure `fig`.

`draw(x, y, color, ax, **kwargs)`\: draw data `x` vs. `y` with color `color` on axes `ax`.

`models(starts, ends, counts, bound, ax, **kwargs)`\: draws gene models over ranges between the values in `starts` and the corresponding values in `ends` and writes the number of genes with that model from `counts` on the axes `ax`.

`report(ntvar, aavar, lnt, laa)`\: prints out the of average nucleotide variance (`ntvar`) per nucleotide (`lnt`) and the average amino acid variance (`aavar`) per nucleotide (`laa`).

`biotools.analysis.predict`
---------------------------

`ORFGenerator(sequence)`\: scans both strands of the given sequence for ORFs, and yields them one by one.

`_genepredict_target(qin, qout, orfs, subj)`\: obtains BLAST results from `qin` and searches `orfs` for ORFs that are in range (see below) and determines the best-aliging ORF subsequence (see `biotools.align`) and adds it to `qout` if it meets the minimum requirements (see `biotools.analysis.options`).

`_genepredict_in_range(seq, start, end)`\: returns True if the given sequence spans the range (`start`, `end`).

`GeneFromBLAST(db, sequences, pref)`\: begins BLAST between `db` and `sequences`. Results are sent to `_genepredict_target`, and good hits are saved to disk as protein, nucleotide, and gff files.

`run(*args)`\: simple wrapper for GeneFromBLAST that verifies the correct number of arguments are passed.

`biotools.analysis.renamer`
---------------------------

`rename(direc, db)`\: searches through `direc` for files with a fastc extension and renames them according to user input. The user is given descriptions of the sequences contained within the fastc file from sequences in the `db` file.

`biotools.analysis.run`
-----------------------

`_run_genepredict(q, infile)`\: runs a single instance of `biotools.analysis.predict`.

`run()`\: runs several (see `biotools.analysis.options`) instances of `_run_genepredict` at once.

`biotools.analysis.variance`
----------------------------

`var(strain, fmt)`\: returns plotdata and metadata for plotting later on in the pipeline.

`SaySNPs(input)`\: takes a clustalw alignment and will return a dictionary of data relevent to plotting the sequence variance for the sequences in the
given clustalw alignment. These data are:
* `var`\: the measure of sequence variation,
* `starts`\: the starting positions for each gene model in amino acids,
* `ends`\: the ending positions for each gene model in amino acids, and
* `counts`\: the number of sequences with a particular gene model.
The values given in `starts`, `ends`, and `counts` are sorted to that the nth element in `starts` corresponds to the nth value in `ends` and the nth value in `counts`.

`prok-geneseek`
===============

`prok-geneseek` is a command-line utility that is installed onto your system when you install biotools via [pip](http://www.pip-installer.org). It glues the pieces of the analysis into a single tool.

    Usage: prok-geneseek [options] <database> <sequences ...>

    Options:
      -h, --help            show this help message and exit
      -S START, --start=START
                            define a start codon [default: -S ATG]
      -E STOP, --stop=STOP  define a stop codon [default: -E TAG -E TAA -E TGA]
      -j THREADS, --threads=THREADS
                            number of threads [default: 16]
      -p PROCESSES, --processes=PROCESSES
                            number of parallel processes to run [default: 2]
      -e EVALUE, --evalue=EVALUE
                            maximum e-value [default: 1e-30]
      -I IDENTITY, --identity=IDENTITY
                            minimum percent identity [default: 0.45]
      -L FRACTION, --length=FRACTION
                            allowable relative error in hit length (|hit-ref|/ref)
                            [default: 0.2]
      -d DIRECTORY, --directory=DIRECTORY
                            minimum percent length [default: current]
      -P PLOTTER, --plotter=PLOTTER
                            plotting module [default: biotools.analysis.plot]
