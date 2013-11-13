#!/usr/bin/env python

'''
A module to manage BLAST databases and interface with the BLAST+ standalone
program available from NCBI.
'''

import biotools.IO as io
import subprocess
from os import sep, getenv, listdir
import shutil


def run(db, sfile, mega_blast=False, **kwargs):
    '''
    Takes a database and a query and runs the appropriate type of BLAST on
    them. The database can be an existing BLAST database or a fasta/fastq
    file. If it is a sequence file, this function will look in the places
    where BLAST would look for an existing database created from that file and
    use that instead. If there is no such database, this function will make
    one for you and then use the newly created database with BLAST.

    Optional named arguments can currently only be `evalue`, `num_threads`,
    `gapopen`, or `gapextend`. The correspond to the BLAST options of the same
    name.
    '''

    cmds = {
        'prot': {
            'prot': 'blastp',
            'nucl': 'tblastn'
        },
        'nucl': {
            'nucl': 'blastn',
            'prot': 'blastx'
        }
    }

    seq = io.open(sfile, 'r').next()
    qtype = seq.type

    rcloc = ''
    for loc in (".:~:" + (getenv("NCBI") or "")).split(':'):
        if loc and loc[-1] == sep:
            loc += sep
        try:
            for line in (l.strip() for l in open(loc + '.ncbirc', 'r')):
                pos = line.find('=')
                if pos >= 0 and line[:pos].strip() == "BLASTDB":
                    rcloc = line[pos + 1:].strip()
        except IOError:
            pass

    dbtype = None
    bdbenv = getenv("BLASTDB")
    dblocations = (":." + ((':' + bdbenv) if bdbenv else '') + 
                   ((':' + rcloc) if rcloc else '')).split(':')
    for loc in dblocations:
        if loc and loc[-1] != sep:
            loc += sep
        try:
            open(loc + db + '.pin', 'r')
            dbtype = 'prot'
            break
        except IOError:
            try:
                open(loc + db + '.nin', 'r')
                dbtype = 'nucl'
                break
            except IOError:
                pass

    if not dbtype:
        odb = db
        pos = db.rfind(".")
        for seq in io.open(db, 'r'):
            dbtype = seq.type
            break
        if not dbtype:
            raise IOError("Database not found: " + odb)

        ndb = None
        sp = db.rfind(sep)
        if sp > -1:
            dbdir, db = db[:sp], db[sp + 1:pos]
        else:
            dbdir, db = '.', db[:pos]

        for file in listdir(dbdir):
            dpos = file.rfind('.')
            if dpos >= 0 and file[dpos + 1:] == dbtype[0] + 'in':
                fh = open(dbdir + sep + file, 'r')
                c = ord(fh.read(12)[-1])
                fname = fh.read(c)
                if fname[0] in ("'", '"'):
                    fname = fname[1:-1]
                if fname.endswith(odb):
                    ndb = dbdir + sep + file[:dpos]
                    break
        if not ndb:
            ndb = '_'.join(db.split())
            try:
                ignore = open('/dev/null', 'w')
            except IOError:
                ignore = open('nul', 'w')

            try:  # possible race condition
                open(ndb, 'r').close()
            except IOError:
                subprocess.call(["makeblastdb", "-in", '"%s"' % odb,
                                 "-out", ndb, "-dbtype", dbtype],
                                 stdout=ignore)
                try:
                    for suff in ['in', 'hr', 'sq']:
                        name = ndb + '.' + dbtype[0] + suff
                        shutil.move(name, dbdir + sep + name)
                except shutil.Error:
                    pass
            db = dbdir + sep + ndb
        else:
            db = ndb
    else:
        raise IOError("Database not found: " + db)
    allowed = set(["evalue", "gapopen", "gapextend", "num_threads"]) & \
        set(kwargs.keys())
    cmd = cmds[qtype][dbtype]
    pn = ["-db", "-query"]
    if mega_blast:
        cmd = "megablast"
        pn = ["-d", "-i"]
        allowed = ["e", "a"]

    proc = subprocess.Popen([cmd, pn[0], db, pn[1], sfile] +
                            [arg for pair in
                             [["-" + k, str(kwargs[k])] for k in allowed]
                             for arg in pair],
                            bufsize=1, stdout=subprocess.PIPE)
    return Result(iter(proc.stdout.readline, ''))


class Result(object):

    '''
    A class which take the raw output from BLAST and generates dictionaries
    from the data from BLAST. This data includes the alignment, percent
    identity, gaps, e-value, score, length of subject, length of query, and
    start and stop positions for both sequences. This class should be used in
    a for loop like so:

    ```python
        for res in Result(file_or_data):
            pass
    ```

    The class instance has a single other property, `headers`, which are the
    lines in BLAST results before the BLAST hits (e.g., citation info, etc.).
    '''

    def __init__(self, file):
        self.file = file
        self.headers = []

    def __iter__(self):
        try:
            ipt = open(self.file, 'r')
        except (IOError, TypeError):
            try:
                ipt = self.file.split('\n')
            except:
                ipt = self.file
        mode = 0
        headers = []
        curr = None
        length = 0

        def sh(sn, qn, l):
            qdl = ''
            space = qn.find(' ')
            if space > -1:
                qn, qdl = qn[:space], qn[space + 1:].lstrip()
            return {
                'subject': {
                    'name':  sn.lstrip(),
                    'defline': '',
                    'start': None,
                    'end':   None,
                    'sequence': ''
                },
                'query': {
                    'name':  qn,
                    'defline': qdl,
                    'start': None,
                    'end':   None,
                    'sequence': ''
                },
                'length':  l
            }

        def ra(sh):
            for res in ('subject', 'query'):
                sh[res]['start'] = int(sh[res]['start'])
                sh[res]['end'] = int(sh[res]['end'])
                sh[res]['length'] = abs(sh[res]['end'] - sh[res]['start'] + 1)
            return sh

        def sh_fmt(l):
            for pairs in (a.strip() for a in l.split(',')):
                l, r = tuple(a.strip() for a in (pairs.split('=')[:2]
                             if '=' in pairs else pairs.split(':')[:2]))
                subheaders[l.lower().split('(')[0]] = r

        for line in ipt:
            line = line.rstrip('\n').lstrip()
            if not line:
                if mode == 4:
                    mode = 5
                continue

            if mode == 0:
                if line[:6] == 'Query=':
                    mode = 1
                    qname = line[6:].lstrip()
                    self.headers = headers
                else:
                    headers.append(line)

            elif mode == 1:
                if line[0] == '>':
                    mode = 3
                    subheaders = sh(line[1:], qname, length)
                elif line[:7] == 'Length=':
                    length = int(''.join(line[7:].strip().split(',')))
                    mode = 2
                elif line[0] == '(' and line.endswith('letters)'):
                    length = int(''.join(line[1:-8].strip().split(',')))
                    mode = 2
                elif line[:6] == 'Query=':
                    qname = line[6:].lstrip()
                else:
                    qname += line

            elif mode == 2:
                if line[0] == '>':
                    mode = 3
                    subheaders = sh(line[1:], qname, length)
                elif line[:6] == 'Query=':
                    qname = line[6:].lstrip()
                    mode = 1

            elif mode == 3:
                if line[:5] == 'Score':
                    snm = subheaders['subject']['name']
                    defline = ''
                    space = snm.find(' ')
                    if space > -1:
                        snm, defline = snm[:space], snm[space + 1:]
                    subheaders['subject']['name'] = snm
                    subheaders['subject']['defline'] = defline
                    sh_fmt(line)
                    mode = 4
                elif line[:7] == 'Length=':
                    pass
                elif line[0] == '(' and line.endswith('letters)'):
                    pass
                else:
                    subheaders['subject']['name'] += line

            elif mode == 4:
                sh_fmt(line)

            elif mode == 5:
                if line[:6] == 'Query=':
                    mode = 1
                    qname = line[6:].lstrip()
                    yield ra(subheaders)
                    continue
                elif line[0] == '>':
                    yield ra(subheaders)
                    subheaders = sh(line[1:], qname, length)
                    mode = 3
                    continue
                elif line[:5] == 'Score':
                    yield ra(subheaders)
                    subheaders = sh(subheaders['subject']['name'], qname,
                                    length)
                    sh_fmt(line)
                    mode = 4
                    continue
                elif line[:5] == 'Sbjct':
                    curr = 'subject'
                elif line[:5] == 'Query':
                    curr = 'query'
                else:
                    continue

                _, start, seq, end = line.split()
                subheaders[curr]['start'] = subheaders[curr]['start'] or start
                subheaders[curr]['end'] = end
                subheaders[curr]['sequence'] += seq

        try:
            yield ra(subheaders)
        except UnboundLocalError:
            pass
        raise StopIteration()

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        output = open(sys.argv[1]).read()
        for result in Result(output):
            print(result)
