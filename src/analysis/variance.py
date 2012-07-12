#!/usr/bin/env python

import biotools.IO as io
import biotools.analysis.options as options


def SaySNPs(input):
    '''
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
    '''

    catalogue = []
    lengths = {}
    for seq in io.open(input, 'r'):
        key = (seq.start, seq.end)
        lengths[key] = lengths.get(key, 0) + 1
        for (i, c) in zip(xrange(key[1] - key[0] + 1), seq):
            if i >= len(catalogue):
                catalogue.append({})
            if c != " ":
                catalogue[i][c] = catalogue[i].get(c, 0) + 1

    calc = []
    for s in catalogue:
        tot = float(sum(s.values()))
        cnt = float(len(s))
        calc.append(1.0 - sum((s[c] / tot) ** 2 for c in s))
    llist = sorted(list(lengths.keys()))

    return {
        'var': calc,
        'starts': [s for s, e in llist],
        'ends': [e for s, e in llist],
        'count': [lengths[k] for k in llist]
    }


def var(strain, fmt):
    '''
    Returns plot data and metadata for plotting later on in the pipeline.
    '''
    plotdata = {
        'nt': SaySNPs(options.DIRECTORY + fmt % {
            'strain': strain,
            'type': 'nt'
        }),
        'aa': SaySNPs(options.DIRECTORY + fmt % {
            'strain': strain,
            'type': 'aa'
        })
    }
    metadata = {'strain': strain, 'filename': strain + '.pdf'}

    return plotdata, metadata
