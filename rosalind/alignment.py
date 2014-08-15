from .io import load_scoring_matrix


def levenshtein(first, second, scoring_matrix=None, gap=1):
    """
    Calculate the levenshtein distance between two sequences. http://en.wikipedia.org/wiki/Levenshtein_distance
    Wraps alignment score, with a negation
    Arguments: Sequence first, Sequence second, str scoring_matrix, int gap
    Returns: int
    """
    return -1 * alignment_score(first, second, scoring_matrix, gap)


def alignment_score(first, second, scoring_matrix=None, gap=1):
    """
    Calculate the alignment score distance between two sequences. http://en.wikipedia.org/wiki/Levenshtein_distance
    Can take a custom scoring matrix and gap penalty
    Arguments: Sequence first, Sequence second, str scoring_matrix, int gap
    Returns: int
    """
    if scoring_matrix is None:
        scoring = _DefaultMatrix()
    else:
        scoring = load_scoring_matrix(scoring_matrix)
    s = first.sequence
    t = second.sequence
    return _alignment_matrix(s, t, scoring, gap)[len(s)][len(t)]


def edit_distance_alignment(first, second, scoring_matrix=None, gap=1):
    """
    Find an optimal alignment minimizing edit distance between two sequences, return the distance and a pair of
    strings representing an optimal alignment.
    Arguments: Sequence first, Sequence second, str scoring_matrix, int gap
    Returns: int, str, str
    """
    if scoring_matrix is None:
        scoring = _DefaultMatrix()
    else:
        scoring = load_scoring_matrix(scoring_matrix)
    s = first.sequence
    t = second.sequence
    matrix = _alignment_matrix(s, t, scoring, gap)
    s1, t1 = _alignment_backtrack(s, t, matrix)
    return -1 * matrix[-1][-1], s1, t1


def optimal_alignment_count(first, second, modulus=None, scoring_matrix=None, gap=1):
    """
    Finds the number of valid optimal alignments for a pair of sequences. As this can get large, can return modulo
    some number.
    Arguments: Sequence first, Sequence second, int modulus, str scoring_matrix, int gap
    Returns: int
    """
    if scoring_matrix is None:
        scoring = _DefaultMatrix()
    else:
        scoring = load_scoring_matrix(scoring_matrix)
    s = first.sequence
    t = second.sequence
    matrix = _alignment_matrix(s, t, scoring, gap)
    return _count_alignments(s, t, matrix, modulus)


def best_local_alignment(first, second, scoring_matrix=None, gap=1):
    """
    Finds the local alignment of a pair of sequences with the shortest edit distance.
    Arguments: Sequence first, Sequence second, str scoring_matrix, int gap
    Returns: ?
    """
    if scoring_matrix is None:
        scoring = _DefaultMatrix()
    else:
        scoring = load_scoring_matrix(scoring_matrix)
    s = first.sequence
    t = second.sequence
    best_local = (float('inf'), 0, 0, 0, 1)
    # there is no way this is the most efficient search
    for i in xrange(len(s)):
        for j in xrange(len(t)):
            if s[i] == t[j]:
                matrix = _alignment_matrix(s[i:], t[j:], scoring, gap)
                for m in xrange(len(matrix)):
                    for n in xrange(len(matrix[0])):
                        candidate = matrix[m][n]
                        if candidate < best_local[0]:
                            best_local = (candidate, i, j, i+m, j+n)
    return best_local[0], s[best_local[1]:best_local[3]], t[best_local[2]:best_local[4]]


def _alignment_matrix(s, t, scoring, gap):
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
        matrix[i][0] = i * -gap
    for i in xrange(n+1):
        matrix[0][i] = i * -gap
    for j in xrange(1, n + 1):
        for i in xrange(1, m + 1):
            matrix[i][j] = max(matrix[i - 1][j] - gap, matrix[i][j - 1] - gap, matrix[i - 1][j - 1] + scoring[(s[i - 1], t[j - 1])])
    return matrix


def _alignment_backtrack(s, t, matrix):
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
        if diag >= left and diag >= up:
            i -= 1
            j -= 1
            s1 = s[i] + s1
            t1 = t[j] + t1
        elif left >= up:
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
                diag = matrix[i - 1][j - 1] + 1
            else:
                diag = matrix[i - 1][j - 1]
            m = max(diag, matrix[i - 1][j], matrix[i][j - 1])
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
            return -1
