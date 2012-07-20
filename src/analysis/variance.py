#!/usr/bin/env python

import biotools.IO as io
import biotools.analysis.options as options
from biotools.translate import translate
from os import sep


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


def var(files):
    '''
    Returns plot data and metadata for plotting later on in the pipeline.
    '''
    sort = {}
    for f in files:
        seqs = [s for s in io.open(f)]
        type = set(s.type for s in seqs)
        if len(type) > 1:
            type = set(['prot'])
        fid = (type.pop(), f)
        seqs = [''.join(s.seq.split('-')).strip() for s in seqs]
        seqs = [translate(s) if fid[0] == 'nucl' else s for s in seqs]
        sset = frozenset(seqs)
        srtr = (len(seqs), sset)
        sort[srtr] = sort.get(srtr, set()) | set([fid])

    couples = []
    for partners in sort.values():
        trim = lambda x: '.'.join(x.split('.')[:-1]) \
                         if f.endswith('.clustalw') or \
                            f.endswith('.clustal') or \
                            f.endswith('.aln') else x
        names = ', '.join(set(trim(f.split(sep)[-1]) for type, f in partners))
        pair = {}
        for type, f in partners:
            if len(pair) == 2:
                break
            if type in pair:
                continue
            pair[type] = f
        if 0 < len(pair) < 2:
            raise TypeError("Unmatched clustal alignment(s): " + 
                            ", ".join(f for type, f in partners))
        if len(pair) == 0:
          continue
        couples.append((pair['nucl'], pair['prot'], names))

    for nt, aa, strain in couples:
        plotdata = {
            'nt': SaySNPs(nt),
            'aa': SaySNPs(aa)
        }
        metadata = {'strain': strain, 'filename': strain + '.pdf'}

        yield {'plotdata': plotdata, 'metadata': metadata}
    raise StopIteration
