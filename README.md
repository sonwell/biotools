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

`Annotation(self, ref, src, type, start, end, score, strand, phase, attr, name\_token='ID', 'gff\_token')`\: constructs an `Annotation` object with the necessary values. The parameters are passed in the same order as the columns from a GFF (version 3) file and the `name\_token` and `gff\_token` parameters are the defaults for a gff version 3 file from phytozome. Just write (e.g.) `Annotation(\*line.split('\\t')) #(splitting on tabs)`, and the rest of the work will be done for you. Other sources may require changes to `name\_tokens` and `gff\_token`.

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

`IOBase` has a class property, `methods`, which is an `IOManager` (see `biotools.iomanager`).

An `IOBase` instance has the following attributes:
* file
* handle
* method
* type
and others that are specific to the particular file type.

`IOBase.format(self, fmt)`\: tries to parse the file as a file of type `fmt`. If it cannot, it will not read the file.

`IOBase.close(self)`\: closes the handle.

`Reader(filename, mode='r')`\: extends `IOBase` to allow the reading of files. Can be used in the iterator protocol (that is, in a `for` loop).

`Reader.read(self, n=None)`\: if `n` is provided, the next (up to) `n` entries are parsed and returned. Otherwise, all remaining entries are parsed and returned.

`Reader.next(self)`\: reads a single entry in the file and returns it.

`Writer(filename, mode='r')`\: extends `IOBase` to allow the writing of files.

`Writer.write(self, sequence)`\: appropriately formats a sequence and writes it to the disk.

`open(filename, mode='r')`\: open a file for parsing or writing; returns either a Reader or Writer, depending on the open mode.

`biotools.iomanager`
--------------------

`IOManager(self, methods=None)`\: a class used by the `IOBase` class to manage the various input and output methods for the different file types. Additional file types can be added to the manager by using `manager[format] = methods`.

From the above example, methods is a dictionary with keys rhook, read, whook    , write, and probe. Each of the values must be callable object:
* rhook\: takes a file handle, opened for reading; called before reading of the file has begun,
* whook\: takes a file handle, opened for writing; called before writing to the file has begun,
* read\: takes a file handle, opened for reading; should be a generator that yields entries,
* write\: takes a file handle, opened for writing and a single entry; writes the entry to the file,
* probe\: takes a file handle, opened for reading; returns a dictionary of attributes to be applied to the `IOBase` instance.

This class behaves similarly to a dictionary, except that the `get` method will default to the default method (which does nothing) if no truthy second parameter is passed.

`biotools.sequence`
-------------------

`chop(seq, length=70)`\: yields a chunk of a sequence of no more than `length` characters, it is meant to be used to print fasta files.

`Sequence(self, name, seq, **kwargs)`\: instantiates a `Sequence` object with sequence `seq`. 

A sequence object has attributes
* name
* seq
* qual
* type
* start
* end
* step
* original
* defline

Some other useful parameters that the Sequence constructor can handle are:
* qual        => the quality scores (an array of integers) of the sequence,
* type        => the type of the sequence, either 'prot' or 'nucl',
* start       => the starting position of the sequence within a supersequence,
* end         => the ending position of the sequnece within a supersequence,
* step        => the 'step' of the sequence, usually +1 for top-strand sequences, and -1 for bottom-strand sequences, but can handle other values as well,
* original    => the original Sequence object from which this one derives,
* defline     => the definition line for this sequence from a fasta file.
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

`biotools.analysis.options`
---------------------------

`biotools.analysis.plot`
------------------------

`biotools.analysis.predict`
---------------------------

`biotools.analysis.renamer`
---------------------------

`biotools.analysis.run`
-----------------------

`biotools.analysis.variance`
----------------------------
