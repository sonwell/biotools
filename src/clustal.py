#!/usr/bin/env python
import biotools.IO as io
import subprocess
from os import remove


def run(infile, outfile, **kwargs):
    n = 0
    for seq in io.open(infile, 'r'):
        n += 1
        if n > 1:
            seqtype = seq.type
            break

    if n > 1:
        cmd = "clustalw"

        try:
            ignore = open('/dev/null', 'w')
        except IOError:
            ignore = open('nul', 'w')

        if seqtype == 'nucl':
            defaults = {
                'OUTORDER': 'ALIGNED',
                'GAPOPEN': '10',
                'GAPEXT': '0.1',
                'DNAMATRIX': 'IUB'
            }
            others = ["-%s=%s" % (arg, kwargs.get(arg, defaults[arg]))
                      for arg in set(name.upper() for name in kwargs) &
                      set(defaults.keys())]
            subprocess.call([cmd, "-INFILE=" + infile, "-ALIGN", "-TYPE=DNA",
                            "-OUTFILE=" + outfile] + others, stdout=ignore)
        else:
            defaults = {
                'OUTORDER': 'ALIGNED',
                'GAPOPEN': '10',
                'GAPEXT': '0.1',
                'MATRIX': 'BLOSUM'
            }
            others = ["-%s=%s" % (arg, kwargs.get(arg, defaults[arg]))
                      for arg in set(name.upper() for name in kwargs) &
                      set(defaults.keys())]
            subprocess.call([cmd, "-INFILE=" + infile, "-ALIGN",
                            "-TYPE=PROTEIN", "-OUTFILE=" + outfile] + others,
                            stdout=ignore)
        pos = infile.rfind('.')
        if pos > -1:
            prefix = infile[:pos]
        else:
            prefix = infile
        remove(prefix + '.dnd')
