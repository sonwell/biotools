'''
Methods for manipulating clustalw alignment files.
'''

from biotools.sequence import Sequence


def read(fh):
    def clean_alignment(x):
        i = 0
        for c in x:
            if c != "-":
                break
            i += 1
        j = 0
        for c in reversed(x):
            if c != "-":
                break
            j += 1
        return (' ' * i + x[(i or None):(-j or None)] + ' ' * j, i, len(x) - j)

    seqs = {}
    for line in fh:
        if line.startswith(' '):
            continue
        st = line.strip()
        if st:
            bits = st.split()
            if len(bits) != 2:
                continue
            if bits[0] not in seqs:
                seqs[bits[0]] = ''
            seqs[bits[0]] += bits[1]

    for k in seqs:
        seq, start, end = clean_alignment(seqs[k])
        yield Sequence(k, seq, start=start, end=end)
    raise StopIteration()


def probe(fh):
    for line in fh:
        st = line.strip()
        if st:
            if st.startswith('CLUSTAL'):
                return {'type': 'clustalw'}
            return False
    return {'type': 'clustalw'}


def rhook(fh):
    try:
        fh.next()
    except StopIteration:
        pass
