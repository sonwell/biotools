#!/usr/bin/python
_ref = {'A':'T','T':'A','C':'G','G':'C'}

import sys

def complement(s):
  '''complement( sequence )
Creates the complement of a sequence, which can then be reversed by using seq[::-1], if it needs to be reversed. This function accepts either Sequences or strings.'''

  value = ''.join(_ref.get(c,'N') for c in s.upper())
  try: return s.__class__("complement(%s)" % s.name, value, original=s.original,start=s.start,end=s.end,step=s.step)
  except AttributeError, TypeError: return s.__class__(value)

if __name__ == "__main__":
  print complement(sys.argv[1][::-1])
