def _parseAttrs(attr,token='='):
  '''_parseAttrs(attr, token) /internal/
Creates a dictionary from the atrributes (9th column) of a gff file. By default, token is '=', which is the separator used in gff version 3.

In other words, attr "a=b;c=d;" and token '=' will yield the dictionary {'a':'b','c':'d'}. The other separator (';') cannot be changed.

This function is not to be called on its own.'''

  attributes = {}
  attrs = [a.strip() for a in attr.strip().split(';')]
  for attribute in attrs:
    pos = attribute.find(token)
    if pos > -1:
      var,val = attribute[:pos],attribute[pos+1:]
      attributes[var] = attributes.get(var,[]) + [val]
	
  for key in attributes:
    attributes[key] = ','.join(attributes[key])
  return attributes

class Annotation(object):
  '''class Annotation
An object to help with reading and writing GFF files.'''
  unknowns = 0
	
  def __init__(self,ref,src,type,start,end,score,strand,phase,attr,name_token='ID',gff_token='='):
    '''Annotation(ref, src, type, start, end, score,
           strand, phase, attr, name_token, gff_token)
Constructs an Annotation object with the necessary values. The parameters are passed in the same order as the columns from a GFF (version 3) file and the name_token and gff_token parameters are the defaults for a gff version 3 file from phytozome. Just write (e.g.)
	Annotation(*line.split('\\t')) #(splitting on tabs),
and the rest of the work will be done for you. Other sources may require changes to name_tokens and gff_token.

Instantiating an Annotation will generate for it an id of the form SEQNAME_TYPE[START:END], where SEQNAME is the name of the sequence (column 1) from the GFF file, and type is like 'gene' or 'CDS'. If no SEQNAME is provided, then 'X' be used in its place, and if no identifier can be found in the attributes, the Annotation will generate an identifier for itself in the form of unknown #.'''

    start,end = int(start),int(end)
    self.strand = strand
    self.type   = type
    self.source = src
    self.seq    = ref
    self.start  = min(start,end)
    self.end    = max(end,start)
    self.attr   = _parseAttrs(attr,gff_token)
    self.phase  = phase
    self.score  = score
    self.ntoken = name_token
    self.id     = ((self.seq or 'X') + '_' + self.type + "[%d:%d]" % (self.start,self.end))
    try:
      self.name   = self.attr[name_token]
    except:
      Annotation.unknowns += 1
      self.name = "unknown %d" % Annotation.unknowns
    self.parent = None
    self.children = []

  '''Some things that you can do to Annotation objects:
* len(annotation)  => gets the length of the annotation (end-start+1)
* dictionary[annotation] => store annotations as keys of a dictionary or as elements in a set
* annA == annB     => compare two Annotations, they are the same if they have the same id.
* print annotation => prints the annotation as a line of a GFF version 3 file.
'''
		
  def __len__(self):
		return max(self.start,self.end)-min(self.end,self.start)+1
    
  def __hash__(self):
    return self.id.__hash__()
    
  def __eq__(self,other):
    try:
      return self.id == other.id
    except:
      return False

  def __str__(self):
    return '\t'.join((self.seq,self.source, \
      self.type,str(self.start),str(self.end),self.score, \
      self.strand,str(self.phase), \
      ';'.join(k+'='+self.attr[k] for k in self.attr)))

