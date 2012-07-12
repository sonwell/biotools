#!/usr/bin/env python
import sys

_ref = {
    'DNA': {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R'},
    'RNA': {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R'}
}


def complement(s):
    '''
    Creates the complement of a sequence, which can then be reversed by using
    `seq[::-1]`, if it needs to be reversed. This function accepts either
    `Sequence`s or strings.
    '''

    if set(s) - set('ATUCGNRY'):
        return s
    if 'U' in s and 'T' not in s:
        repl = _ref['RNA']
    elif 'T' in s and 'U' not in s:
        repl = _ref['DNA']
    else:
        repl = _ref['DNA']
    value = ''.join(repl.get(c, 'N') for c in s.upper())
    try:
        return s.__class__("complement(%s)" % s.name, value,
                           original=s.original, start=s.start,
                           end=s.end, step=s.step, qual=s.qual)
    except (AttributeError, TypeError):
        return s.__class__(value)


if __name__ == '__main__':
    assert complement('ATCGTAGCTGATCGAT') == 'TAGCATCGACTAGCTA'
    assert complement('AUCGUAGCUGAUCGAU') == 'UAGCAUCGACUAGCUA'
    assert complement('ATCGUAGCUGAUCGAU') == 'ATCGUAGCUGAUCGAU'
