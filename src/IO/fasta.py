from biotools.sequence import Sequence


def read(fh):
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
    fh.close()
    raise StopIteration()


def write(fh, s):
    fh.write('>%s %s\n' % (s.name, s.defline) +
             '\n'.join(chop(s.seq, 70)) + '\n')


def probe(fh):
    for line in fh:
        st = line.strip()
        if st:
            fh.close()
            if st[0] == '>':
                return {'type': 'fasta'}
            return False
    fh.close()
    return {'type': 'fasta'}
