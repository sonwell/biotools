`biotools`
==========

A bunch of bioinformatics tools.

`biotools.align`
----------------

`biotools.align.OptimalCTether(reference, translation)`\: uses a hybrid between the Needleman-Wunsch and Smith-Waterman algorithms to find the optimal sub-string of `translation` that begins with a start codon (see also: `biotools.analysis.options`) aligning to `reference` by requiring global alignment at the C terminus and local alignment at the N terminus. 

`OptimalCTether` returns a dictionary of relevent information from the alignment; specifically, the alignments itself [keys: `query`, `subject`], the score [key: `score`], the length of the alignment [key: `length`], the length of the substring of translation used [key: `sublength`], the number of identities [key: `identities`], the theoretical perfect score for that alignment [key: `perfect`], and the number of gaps [key: `gaps`].

`biotools.annotation`
---------------------

`biotools.BLAST`
----------------

`biotools.blosum62`
-------------------

`biotools.complement`
---------------------

`biotools.IO`
-------------

`biotools.iomanager`
--------------------

`biotools.sequence`
-------------------

`biotools.translate`
--------------------

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