from biotools.annotation import Annotation


def read(fh):
    for line in fh:
        if line[0] != '#':
            yield Annotation(*line.split('\t'))
    raise StopIteration()


def write(fh, a):
    # TODO: improve this...
    fh.write(str(a) + '\n')


def probe(fh):
    for line in fh:
        line = line.strip()
        if line:
            bits = line.split()
            if bits[0] == '##gff-version':
                return {'type': 'gff', 'version': float(bits[1])}
            return False
    return {'type': 'gff', 'version': 3}


def whook(fh):
    fh.write('##gff-version 3\n')
