'''
Functions for manipulating PHYLIP files.
'''

from biotools.sequence import Sequence, chop


def read(fh):
    '''
    Read sequences in PHYLIP format; This function is a generator that yields
    `Sequence` objects.
    '''
    grid, seqs, names = None, [], []
    
    try:
        line = fh.next().strip()
        while not line:
            line = fh.next.strip()
    except StopIteration:
        fh.close()
        raise StopIteration()
    grid = [int(x) for x in line.strip()]
    while True:
        lines = []
        try:
            for i in xrange(grid[0]):
                line = fh.next().strip()
                while not line:
                    line = fh.next().strip()
                lines.append(line)
            if not names:
                names = [l.split()[0] for l in lines]
                try:
                    seqs = [''.join(l.split()[1:]) for l in lines]
                except IndexError:
                    seqs = [''] * grid[0]
            else:
                temp = [''.join(l.split()) for l in lines]
                seqs = [a + b for a, b in zip(seqs, temp)]
        except StopIteration:
            break
    for name, seq in zip(names, seqs):
        yield Sequence(name, seq)
    fh.close()
    raise StopIteration()


def probe(fh):
    '''
    '''
    for line in fh:
        bits = line.split()
        if len(bits) == 2 and int(bits[0]) > 0 and int(bits[1]) > 0:
            return {'type': 'phylip'}
        else:
            return False
    return {'type': 'phylip'}
