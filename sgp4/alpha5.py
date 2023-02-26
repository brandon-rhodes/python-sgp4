"""Alpha 5 encoding of satellite numbers to fit in 5 characters."""

def to_alpha5(n):
    if n < 100000:
        return '%05d' % n
    if n > 339999:
        raise ValueError('satellite number cannot exceed 339999,'
                         " whose Alpha 5 encoding is 'Z9999'")
    i, n = divmod(n, 10000)
    i += ord('A') - 10
    if i >= ord('I'): i += 1
    if i >= ord('O'): i += 1
    return '%c%04d' % (i, n)

def from_alpha5(s):
    if not s[0].isalpha():
        return int(s)
    c, s = s[0], s[1:]
    n = ord(c) - ord('A') + 10
    n -= c > 'I'
    n -= c > 'O'
    return n * 10000 + int(s)
