_char = ""
_pos  = -1

def parse(ipt):
	global _char, _pos

	_char = ""
	_pos  = -1
	leng  = 0
	val   = ""
	ret   = {}

	def advance(test=None):
		global _pos, _char
		curr = _char
		if test and test != curr:
			raise ValueError("Expected %s, saw %s." % (test, curr))

		_pos += 1
		if _pos < leng:
			_char = val[_pos]
		else:
			_char = ''
	
	def whitespace():
		while _char and ord(_char) <= ord(' '):
			advance()

	def variable():
		whitespace()
		start = _pos
		while \
				('0' <= _char <= '9') or \
				('A' <= _char <= 'Z') or \
				('a' <= _char <= 'z'):
			advance()
		stop = _pos
		whitespace()
		advance('=')

		return (val[start:stop], value())

	def value():
		whitespace()
		if _char == '{': return dictionary()
		if _char == '[': return array()
		if _char in ("'", '"'): return string()
		if '0' <= _char <= '9' or \
			_char == '-': return number()
		if _char == 'T':
			advance('r')
			advance('u')
			advance('e')
			return True
		if _char == 'F':
			advance('a')
			advance('l')
			advance('s')
			advance('e')
			return False
		raise ValueError("Unexpected value starting with %s." % _char)

	def dictionary():
		ret = {}
		advance('{')
		whitespace()
		while _char != '}':
			s = value()
			advance(':')
			v = value()
			ret[s] = v
			if _char == ',':
				advance(',')
			else:
				break
		advance('}')
		whitespace()
		return ret

	def array():
		ret = []
		advance('[')
		whitespace()
		while _char != ']':
			ret.append(value())
			if _char == ',':
				advance(',')
			else:
				break
		advance(']')
		whitespace()
		return ret

	def string():
		ret = ""
		q = _char
		advance(q)
		while _char != q:
			if _char == '\\':
				advance('\\')
				if _char == 'n': ret += '\n'
				elif _char == 't': ret += '\t'
				else: ret += _char
				advance()
				continue
			ret += _char
			advance()
		advance(q)
		whitespace()
		return ret

	def number():
		ret = ""
		if _char == '-':
			ret += '-'
			advance('-')
		while '0' <= _char <= '9':
			ret += _char
			advance()
		if _char == '.':
			ret += '.'
			advance('.')
			while '0' <= _char <= '9':
				ret += _char
				advance()
		if _char in ('e', 'E'):
			ret += _char
			advance()
			if _char == '-':
				ret += '-'
				advance('-')
			while '0' <= _char <= '9':
				ret += _char
				advance()
		whitespace()
		return float(ret)

	for line in ipt:
		val  += line
		leng += len(line)

	advance()
	while _pos < leng:
		k, v = variable()
		ret[k] = v

	return ret
