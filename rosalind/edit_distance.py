import sys

def levenshtein(first, second):
    """
    Calculate the levenshtein distance between two sequences. http://en.wikipedia.org/wiki/Levenshtein_distance
    Arguments: Sequence first, Sequence second
    Returns: int
    """
    s = first.sequence
    t = second.sequence
    return _levenshtein_matrix(s, t)[len(s)][len(t)]


def edit_distance_alignment(first, second):
    s = first.sequence
    t = second.sequence
    matrix = _levenshtein_matrix(s, t)
    s1, t1 = _levenshtein_backtrack(s, t, matrix)
    return matrix[-1][-1], s1, t1


def optimal_alignment_count(first, second, modulus=None):
    sys.setrecursionlimit(10000)
    s = first.sequence
    t = second.sequence
    matrix = _levenshtein_matrix(s, t)
    memo = _Memoizer()
    memo.set_params(matrix, s, t, modulus)
    return memo[(len(s), len(t))]


def _levenshtein_matrix(s, t):
    """
    Calculate a matrix showing the levenshtein distance between all prefixes of two strings, with the last element
    being the levenshtein distance between the whole two strings
    Arguments: str s, str t
    Returns: int[][]
    """
    m = len(s)
    n = len(t)
    matrix = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    for i in xrange(m+1):
        matrix[i][0] = i
    for i in xrange(n+1):
        matrix[0][i] = i
    for j in xrange(1, n + 1):
        for i in xrange(1, m + 1):
            if s[i - 1] == t[j - 1]:
                matrix[i][j] = matrix[i - 1][j - 1]
            else:
                matrix[i][j] = min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1, matrix[i - 1][j - 1] + 1)
    return matrix


def _levenshtein_backtrack(s, t, matrix):
    """
    Given two strings and the levenshtein distance matrix between the two, backtrack through and create an alignment.
    Arguments: str s, str t, int[][] matrix
    Returns: str, str
    """
    s1 = ""
    t1 = ""
    i = len(s)
    j = len(t)
    while i > 0 or j > 0:
        if i == 0:
            up = len(s)
        else:
            up = matrix[i - 1][j]
        if j == 0:
            left = len(t)
        else:
            left = matrix[i][j - 1]
        if j == 0 or i == 0:
            diag = max(len(s), len(t))
        else:
            diag = matrix[i - 1][j - 1]
        if diag <= left and diag <= up:
            i -= 1
            j -= 1
            s1 = s[i] + s1
            t1 = t[j] + t1
        elif left <= up:
            j -= 1
            s1 = "-" + s1
            t1 = t[j] + t1
        else:
            i -= 1
            s1 = s[i] + s1
            t1 = "-" + t1
    return s1, t1


class _Memoizer(dict):
    def set_params(self, matrix, s, t, modulus=None):
        self.matrix = matrix
        self.s = s
        self.t = t
        self.modulus = modulus

    def __missing__(self, key):
        i, j = key
        if self.matrix[i][j] == 0:
            self[key] = 1
            return 1
        if i == 0:
            self[key] = self[(i, j - 1)]
        elif j == 0:
            self[key] = self[(i - 1, j)]
        else:
            if self.s[i - 1] == self.t[j - 1]:
                diag = self.matrix[i - 1][j - 1] - 1
            else:
                diag = self.matrix[i - 1][j - 1]
            m = min(diag, self.matrix[i - 1][j], self.matrix[i][j - 1])
            self[key] = 0
            if diag == m:
                self[key] += self[(i - 1, j - 1)]
            if self.matrix[i - 1][j] == m:
                self[key] += self[(i - 1, j)]
            if self.matrix[i][j - 1] == m:
                self[key] += self[(i, j - 1)]
        if self.modulus is not None:
            self[key] %= self.modulus
        return self[key]


class _DefaultMatrix(dict):
    def __missing__(self, key):
        if key[0] == key[1]:
            return 0
        else:
            return 1
