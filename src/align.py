#!/usr/bin/env python

'''
This module is used to align sequences. Currently, there is only a single
alignment algorithm implementented; it is a hybrid between Needleman-Wunsch
and Smith-Waterman and is used to find the subsequence within a larger sequence
that best aligns to a reference.
'''

from biotools.translate import translate
import biotools.analysis.options as options

DIAG_MARK, VGAP_MARK, HGAP_MARK = 3, 2, 1
bl = {
 '*': {'*': 0, 'A': -9, 'C': -9, 'E': -9, 'D': -9, 'G': -9, 'F': -9, 'I': -9,
       'H': -9, 'K': -9, 'M': -9, 'L': -9, 'N': -9, 'Q': -9, 'P': -9, 'S': -9,
       'R': -9, 'T': -9, 'W': -9, 'V': -9, 'Y': -9, 'X': 0},
 'A': {'*': -9, 'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1,
       'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': 1,
       'R': -1, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': 0},
 'C': {'*': -9, 'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1,
       'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1,
       'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2, 'X': 0},
 'E': {'*': -9, 'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3,
       'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0,
       'R': 0, 'T': 0, 'W': -3, 'V': -3, 'Y': -2, 'X': 0},
 'D': {'*': -9, 'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3,
       'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0,
       'R': -2, 'T': 1, 'W': -4, 'V': -3, 'Y': -3, 'X': 0},
 'G': {'*': -9, 'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4,
       'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': -2, 'Q': -2, 'P': -2, 'S': 0,
       'R': -2, 'T': 1, 'W': -2, 'V': 0, 'Y': -3, 'X': 0},
 'F': {'*': -9, 'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0,
       'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2,
       'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3, 'X': 0},
 'I': {'*': -9, 'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4,
       'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2,
       'R': -3, 'T': -2, 'W': -3, 'V': 1, 'Y': -1, 'X': 0},
 'H': {'*': -9, 'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': -2, 'F': -1, 'I': -3,
       'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1,
       'R': 0, 'T': 0, 'W': -2, 'V': -2, 'Y': 2, 'X': 0},
 'K': {'*': -9, 'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3,
       'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0,
       'R': 2, 'T': 0, 'W': -3, 'V': -3, 'Y': -2, 'X': 0},
 'M': {'*': -9, 'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1,
       'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1,
       'R': -1, 'T': -1, 'W': -1, 'V': -2, 'Y': -1, 'X': 0},
 'L': {'*': -9, 'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2,
       'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2,
       'R': -2, 'T': -2, 'W': -2, 'V': 3, 'Y': -1, 'X': 0},
 'N': {'*': -9, 'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3,
       'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1,
       'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2, 'X': 0},
 'Q': {'*': -9, 'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3,
       'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0,
       'R': 1, 'T': 0, 'W': -2, 'V': -2, 'Y': -1, 'X': 0},
 'P': {'*': -9, 'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3,
       'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -1, 'Q': -1, 'P': 7, 'S': -1,
       'R': -2, 'T': 1, 'W': -4, 'V': -2, 'Y': -3, 'X': 0},
 'S': {'*': -9, 'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2,
       'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4,
       'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2, 'X': 0},
 'R': {'*': -9, 'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3,
       'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1,
       'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2, 'X': 0},
 'T': {'*': -9, 'A': -1, 'C': -1, 'E': 0, 'D': 1, 'G': 1, 'F': -2, 'I': -2,
       'H': 0, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1,
       'R': -1, 'T': 4, 'W': -3, 'V': -2, 'Y': -2, 'X': 0},
 'W': {'*': -9, 'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3,
       'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3,
       'R': -3, 'T': -3, 'W': 11, 'V': -3, 'Y': 2, 'X': 0},
 'V': {'*': -9, 'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3,
       'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2,
       'R': -3, 'T': -2, 'W': -3, 'V': 4, 'Y': -1, 'X': 0},
 'Y': {'*': -9, 'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1,
       'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2,
       'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7, 'X': 0},
 'X': {'*': 0, 'A': 0, 'C': 0, 'E': 0, 'D': 0, 'G': 0, 'F': 0, 'I': 0,
       'H': 0, 'K': 0, 'M': 0, 'L': 0, 'N': 0, 'Q': 0, 'P': 0, 'S': 0,
       'R': 0, 'T': 0, 'W': 0, 'V': 0, 'Y': 0, 'X': 0}
}


def OptimalCTether(reference, translation, extend=1, create=10):
    '''
    This function will take two sequences: a `reference` sequence and  another
    protein sequence (`translation`; usually, this is an open reading frame
    that has been translated). Needleman-Wunsch alignment will be performed
    and the substring of translation with the highest identity that begins
    with a start codon [default: `['ATG']`] is reported.

    This function returns a dictionary of relevent information from the
    alignment; specifically, the alignments itself [keys: `query`, `subject`],
    the score [key: `score`], the length of the alignment [key: `length`], the
    length of the substring of translation used [key: `sublength`], the number
    of identities [key: `identities`], and the number of gaps [key: `gaps`].
    '''

    starts = set(translate(s) for s in options.START_CODONS)
    v, w = reference, translation

    try:
        v = v.seq
    except AttributeError:
        pass
    try:
        w = w.seq
    except AttributeError:
        pass
    if not starts & set(w):
        raise ValueError("Open reading frame does not contain a start codon.")

    v, w = v[::-1], w[::-1]
    lv, lw = len(v), len(w)
    rv, rw = range(lv + 1), range(lw + 1)
    gpc = [[create * int(not (i | j)) for i in rw] for j in rv]
    mat = [[-(i + j) * extend - create * (not (i | j) and w[0] != v[0])
           for i in rw] for j in rv]
    pnt = [[VGAP_MARK if i > j else HGAP_MARK if j > i else DIAG_MARK
           for i in rw] for j in rv]
    ids = [[0 for i in rw] for j in rv]
    optimal = [None, 0, 0]
    for i in range(lv):
        for j in range(lw):
            vals = [[mat[i][j] + bl[v[i]][w[j]], DIAG_MARK],
                    [mat[i + 1][j] - extend - gpc[i + 1][j], VGAP_MARK],
                    [mat[i][j + 1] - extend - gpc[i][j + 1], HGAP_MARK]]
            mat[i + 1][j + 1], pnt[i + 1][j + 1] = max(vals)
            gpc[i + 1][j + 1] = create * int(pnt[i + 1][j + 1] == DIAG_MARK)
            if (optimal[0] is None or mat[i + 1][j + 1] > optimal[0]) and \
                    abs(lv - i) / float(lv) <= options.LENGTH_ERR and \
                    w[j] in starts:
                optimal = [mat[i + 1][j + 1], i + 1, j + 1]

    i, j = optimal[1], optimal[2]
    seq, ids = ['', ''], 0
    gapcount, length, sublen = 0, 0, 0
    methods = {
        VGAP_MARK:
            lambda s, i, j, l, g, n:
                (['-' + s[0], w[j - 1] + s[1]], i, j - 1, l + 1, g + 1, n),
        DIAG_MARK:
            lambda s, i, j, l, g, n:
                ([v[i - 1] + s[0], w[j - 1] + s[1]], i - 1, j - 1,
                 l + 1, g, n + (w[j - 1] == v[i - 1])),
        HGAP_MARK:
            lambda s, i, j, l, g, n:
                ([v[i - 1] + s[0], '-' + s[1]], i - 1, j, l, g + 1, n)
    }

    while [i, j] != [0, 0]:
        length += 1
        state = (seq, i, j, sublen, gapcount, ids)
        seq, i, j, sublen, gapcount, ids = methods[pnt[i][j]](*state)

    return {
        'subject': seq[0][::-1],
        'query': seq[1][::-1],
        'score': optimal[0],
        'gaps': gapcount,
        'length': length,
        'sublength': sublen,
        'identities': ids
    }
