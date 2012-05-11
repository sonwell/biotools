'''This is a pretty uninteresting file, the BLOSUM62 matrix I think I ripped from wikipedia, and just
hacked it to pieces to make it do what I want. Lookup looks something like blosum62['A','C'] and would give the value corresponding to the replacement of alanine with cystine.

Basically, you don't need to be using this file, it is just to help in align.py.'''

class blosum62(object):
  def __init__(self):
    self.matrix = [int(x) for x in (\
      ' 9 -1 -1 -3  0 -3 -3 -3 -4 -3 -3 -3 -3 -1 -1 -1 -1 -2 -2 -2 0 -9 ' + \
      '-1  4  1 -1  1  0  1  0  0  0 -1 -1  0 -1 -2 -2 -2 -2 -2 -3 0 -9 ' + \
      '-1  1  4  1 -1  1  0  1  0  0  0 -1  0 -1 -2 -2 -2 -2 -2 -3 0 -9 ' + \
      '-3 -1  1  7 -1 -2 -1 -1 -1 -1 -2 -2 -1 -2 -3 -3 -2 -4 -3 -4 0 -9 ' + \
      ' 0  1 -1 -1  4  0 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -2 -3 0 -9 ' + \
      '-3  0  1 -2  0  6 -2 -1 -2 -2 -2 -2 -2 -3 -4 -4  0 -3 -3 -2 0 -9 ' + \
      '-3  1  0 -2 -2  0  6  1  0  0 -1  0  0 -2 -3 -3 -3 -3 -2 -4 0 -9 ' + \
      '-3  0  1 -1 -2 -1  1  6  2  0 -1 -2 -1 -3 -3 -4 -3 -3 -3 -4 0 -9 ' + \
      '-4  0  0 -1 -1 -2  0  2  5  2  0  0  1 -2 -3 -3 -3 -3 -2 -3 0 -9 ' + \
      '-3  0  0 -1 -1 -2  0  0  2  5  0  1  1  0 -3 -2 -2 -3 -1 -2 0 -9 ' + \
      '-3 -1  0 -2 -2 -2  1  1  0  0  8  0 -1 -2 -3 -3 -2 -1  2 -2 0 -9 ' + \
      '-3 -1 -1 -2 -1 -2  0 -2  0  1  0  5  2 -1 -3 -2 -3 -3 -2 -3 0 -9 ' + \
      '-3  0  0 -1 -1 -2  0 -1  1  1 -1  2  5 -1 -3 -2 -3 -3 -2 -3 0 -9 ' + \
      '-1 -1 -1 -2 -1 -3 -2 -3 -2  0 -2 -1 -1  5  1  2 -2  0 -1 -1 0 -9 ' + \
      '-1 -2 -2 -3 -1 -4 -3 -3 -3 -3 -3 -3 -3  1  4  2  1  0 -1 -3 0 -9 ' + \
      '-1 -2 -2 -3 -1 -4 -3 -4 -3 -2 -3 -2 -2  2  2  4  3  0 -1 -2 0 -9 ' + \
      '-1 -2 -2 -2  0 -3 -3 -3 -2 -2 -3 -3 -2  1  3  1  4 -1 -1 -3 0 -9 ' + \
      '-2 -2 -2 -4 -2 -3 -3 -3 -3 -3 -1 -3 -3  0  0  0 -1  6  3  1 0 -9 ' + \
      '-2 -2 -2 -3 -2 -3 -2 -3 -2 -1  2 -2 -2 -1 -1 -1 -1  3  7  2 0 -9 ' + \
      '-2 -3 -3 -4 -3 -2 -4 -4 -3 -2 -2 -3 -3 -1 -3 -2 -3  1  2 11 0 -9 ' + \
      ' '.join(['0']*22) + ' ' + \
      ' '.join(['-9']*20 + ['0']*2) \
    ).split()]
    acids = 'C S T P A G N D E Q H R K M I L V F Y W X *'.split()
    self.order = dict([(acids[i],i) for i in range(len(acids))])

  def __getitem__(self,q):
    return self.matrix[self.order[q[0]]*22+self.order[q[1]]]
