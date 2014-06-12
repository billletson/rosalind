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
    """
    Find an optimal alignment minimizing edit distance between two sequences, return the distance and a pair of
    strings representing an optimal alignment.
    Arguments: Sequence first, Sequence second
    Returns: int, str, str
    """
    s = first.sequence
    t = second.sequence
    matrix = _levenshtein_matrix(s, t)
    s1, t1 = _levenshtein_backtrack(s, t, matrix)
    return matrix[-1][-1], s1, t1


def optimal_alignment_count(first, second, modulus=None):
    """
    Finds the number of valid optimal alignments for a pair of sequences. As this can get large, can return modulo
    some number.
    Arguments: Sequence first, Sequence second, modulus int or None
    Returns: int
    """
    sys.setrecursionlimit(10000)
    s = first.sequence
    t = second.sequence
    matrix = _levenshtein_matrix(s, t)
    return _count_alignments(s, t, matrix, modulus)


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


def _count_alignments(s, t, matrix, modulus=None):
    current = [1] * (len(t) + 1)
    for i in xrange(1, len(s) + 1):
        prev = current
        current = [1] * (len(t) + 1)
        for j in xrange(1, len(t) + 1):
            if s[i - 1] == t[j - 1]:
                diag = matrix[i - 1][j - 1] - 1
            else:
                diag = matrix[i - 1][j - 1]
            m = min(diag, matrix[i - 1][j], matrix[i][j - 1])
            alignments = 0
            if diag == m:
                alignments += prev[j - 1]
            if matrix[i - 1][j] == m:
                alignments += prev[j]
            if matrix[i][j - 1] == m:
                alignments += current[j - 1]
            current[j] = alignments
            if modulus is not None:
                current[j] %= modulus
    return current[-1]


class _DefaultMatrix(dict):
    def __missing__(self, key):
        if key[0] == key[1]:
            return 0
        else:
            return 1
