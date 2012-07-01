'''
This is a pretty simple JSON-like parser. Specifically, it can load Python-like
object, list, and other literals, i.e., the sort of stuff you'd get it you
dumped the the string representation of some data into a file.

The real difference is that you must specify a variable name, e.g.:

```python
my_stuff = { ... }
```

These variable names don't need to be on a newline or anything like that, you
should be able to omit any and all whitespace. The result of a successful 
parse is a dictionary:

```python
{'my_stuff': { ... }}
```

This function really only works for `None`, `True`, `False`, numbers, strings, 
dictionaries, and lists.
'''

def parse(ipt):
    par = {'char': '', 'pos': -1}
    leng = 0
    val = ""
    ret = {}

    def advance(test=None):
        curr = par['char']
        if test and test != curr:
            raise ValueError("Expected %s, saw %s." % (test, curr))

        par['pos'] += 1
        if par['pos'] < leng:
            par['char'] = val[par['pos']]
        else:
            par['char'] = ''

    def whitespace():
        while par['char'] and ord(par['char']) <= ord(' '):
            advance()

    def variable():
        whitespace()
        start = par['pos']
        while \
                ('0' <= par['char'] <= '9') or \
                ('A' <= par['char'] <= 'Z') or \
                ('a' <= par['char'] <= 'z'):
            advance()
        stop = par['pos']
        whitespace()
        advance('=')

        return (val[start:stop], value())

    def value():
        whitespace()
        if par['char'] == '{':
            return dictionary()
        if par['char'] == '[':
            return array()
        if par['char'] in ("'", '"'):
            return string()
        if '0' <= par['char'] <= '9' or par['char'] == '-':
            return number()
        at = par['pos']
        if par['char'] == 'T':
            advance('r')
            advance('u')
            advance('e')
            return True
        if par['char'] == 'F':
            advance('a')
            advance('l')
            advance('s')
            advance('e')
            return False
        if par['char'] == 'N':
            advance('o')
            advance('n')
            advance('e')
            return None
        prefix = val[at:par['pos']]
        raise ValueError("Unexpected value starting with %s." % prefix)

    def dictionary():
        ret = {}
        advance('{')
        whitespace()
        while par['char'] != '}':
            s = value()
            advance(':')
            v = value()
            ret[s] = v
            if par['char'] == ',':
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
        while par['char'] != ']':
            ret.append(value())
            if par['char'] == ',':
                advance(',')
            else:
                break
        advance(']')
        whitespace()
        return ret

    def string():
        ret = ""
        q = par['char']
        advance(q)
        while par['char'] != q:
            if par['char'] == '\\':
                advance('\\')
                if par['char'] == 'n':
                    ret += '\n'
                elif par['char'] == 't':
                    ret += '\t'
                else:
                    ret += par['char']
                advance()
                continue
            ret += par['char']
            advance()
        advance(q)
        whitespace()
        return ret

    def number():
        ret = ""
        if par['char'] == '-':
            ret += '-'
            advance('-')
        while '0' <= par['char'] <= '9':
            ret += par['char']
            advance()
        if par['char'] == '.':
            ret += '.'
            advance('.')
            while '0' <= par['char'] <= '9':
                ret += par['char']
                advance()
        if par['char'] in ('e', 'E'):
            ret += par['char']
            advance()
            if par['char'] == '-':
                ret += '-'
                advance('-')
            while '0' <= par['char'] <= '9':
                ret += par['char']
                advance()
        whitespace()
        return float(ret)

    for line in ipt:
        val += line
        leng += len(line)

    advance()
    while par['pos'] < leng:
        k, v = variable()
        ret[k] = v

    return ret
