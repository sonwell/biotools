from biotools.sequence import Sequence


def read(fh):
    while 1:
        try:
            line = fh.next().strip()
            if line[0] == '@':
                name = line[1:].split()[0]
                defline = line[1 + len(name):].strip()
                seq = f.next().strip()
                fh.next()
                qual = [ord(c) - self.phred for c in fh.next().strip()]
                yield Sequence(name, seq, qual=qual, defline=defline)
        except StopIteration:
            raise
        finally:
            fh.close()


def write(fh, s):
    fh.write('@%s %s\n%s\n+\n%s\n' % (s.name, s.defline, s.seq,
             ''.join(q + chr('A') - 1 for q in s.qual)) + '\n')


def probe(fh):
    for line in fh:
        st = line.strip()
        if st:
            fh.close()
            if st[0] == '@':
                fh.next()
                fh.next()
                qual = [ord(c) for c in fh.next().strip()]
                phred = 32 if min(qual) < ord('A') else 64
                qual = [q - phred for q in qual]
                return {'type': 'fastq', 'phred': phred}
            return False
    return {'type': 'fastq', 'phread': 64}
