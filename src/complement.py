#!/usr/bin/env python
from biotools.sequence import isprot

_ref = {
    'DNA': {
        'A': 'T', 'T': 'A', 'a': 't', 't': 'a', 
        'C': 'G', 'G': 'C', 'c': 'g', 'g': 'c',
        'R': 'Y', 'Y': 'R', 'r': 'y', 'y': 'r',
        ' ': ' ', '-': '-'
    },
    'RNA': {
        'A': 'U', 'U': 'A', 'a': 'u', 'u': 'a',
        'C': 'G', 'G': 'C', 'c': 'g', 'g': 'c',
        'R': 'Y', 'Y': 'R', 'r': 'y', 'y': 'r',
        ' ': ' ', '-': '-'
    }
}


def complement(s):
    '''
    Creates the complement of a sequence, which can then be reversed by using
    `seq[::-1]`, if it needs to be reversed. This function accepts either
    `Sequence`s or strings.
    '''

    if isprot(s):
        return s
    has_u = ('U' in s or 'u' in s)
    has_t = ('T' in s or 't' in s)
    if has_u and not has_t:
        repl = _ref['RNA']
    elif has_t and not has_u:
        repl = _ref['DNA']
    else:
        repl = _ref['DNA']
    value = ''.join(repl.get(c, 'N') for c in s)
    try:
        return s.__class__("complement(%s)" % s.name, value,
                           original=s.original, start=s.start,
                           end=s.end, step=s.step, qual=s.qual)
    except (AttributeError, TypeError):
        return s.__class__(value)


if __name__ == '__main__':
    assert complement('ATCGTAGCTGATCGAT') == 'TAGCATCGACTAGCTA'
    assert complement('AUCGUAGCUGAUCGAU') == 'UAGCAUCGACUAGCUA'
    print(complement('AUCgu--cuGAUCGAU'))
