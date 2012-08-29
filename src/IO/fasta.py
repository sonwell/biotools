'''
Functions for manipulating FASTA files.
'''

from biotools.sequence import Sequence, chop


def read(fh):
    '''
    Read sequences in FASTA format; identifiers (names and definition lines) 
    are on lines that begin with carets ('>') and sequence is on lines that 
    intervene between the carets. This function is a generator that yields
    `Sequence` objects.
    '''
    name, defline, seq = '', '', ''
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            if name or seq:
                yield Sequence(name, seq, defline=defline)
            seq = ''
            name = line[1:].split()[0]
            defline = line[1 + len(name):].strip()
            continue
        seq += line
    if name or seq:
        yield Sequence(name, seq, defline=defline)
    raise StopIteration()


def write(fh, s):
    '''
    Write sequences in FASTA format, i.e.,
    
    ```
    >name defline
    sequence ...
    ```
    
    Sequences are wrapped to 70 characters by default.
    '''
    fh.write('>%s %s\n' % (s.name, s.defline) +
             '\n'.join(chop(s.seq, 70)) + '\n')


def probe(fh):
    '''
    Probe a file to determine whether or not it is a FASTA file. That is,
    the first non-empty line should begin with a caret ('>'). If no caret is
    found on the first line, then we conclude that it is not a FASTA file
    and return False, otherwise, we return a dictionary with information
    relevant to the FASTA file type.
    '''
    for line in fh:
        st = line.strip()
        if st:
            fh.close()
            if st[0] == '>':
                return {'type': 'fasta'}
            return False
    fh.close()
    return {'type': 'fasta'}
