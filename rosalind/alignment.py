from .io import load_scoring_matrix


class Directions:
    undefined, diag, up, left = range(4)

    def __init__(self):
        pass


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
    Calculate the alignment score distance between two sequences.
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
    if type(gap) is int:
        matrix = _alignment_matrix(s, t, scoring, gap)
    else:
        matrix, pointer = _affine_alignment_matrix(s, t, scoring, gap)
    return matrix[-1][-1]


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
    matrix = _alignment_matrix(s, t, scoring, gap, True)
    start = _array_max_index(matrix)
    s1, t1 = _alignment_backtrack(s, t, matrix, True, "", start)
    return matrix[start[0]][start[1]], s1, t1


def _alignment_matrix(s, t, scoring, gap, local=False):
    """
    Calculate a matrix showing the levenshtein distance between all prefixes of two strings, with the last element
    being the levenshtein distance between the whole two strings
    Arguments: str s, str t
    Returns: int[][]
    """
    m = len(s)
    n = len(t)
    matrix = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    if not local:
        for i in xrange(m+1):
            matrix[i][0] = i * -gap
        for i in xrange(n+1):
            matrix[0][i] = i * -gap
    for j in xrange(1, n + 1):
        for i in xrange(1, m + 1):
            matrix[i][j] = max(matrix[i - 1][j] - gap, matrix[i][j - 1] - gap, matrix[i - 1][j - 1] + scoring[(s[i - 1], t[j - 1])])
            if local and matrix[i][j] < 0:
                matrix[i][j] = 0
    return matrix


def _affine_alignment_matrix(s, t, scoring, gap, local=False):
    m = len(s)
    n = len(t)
    open_gap, extend_gap = gap
    scores = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    f = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    ii = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    ij = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    pointer = [[Directions.undefined for i in xrange(n + 1)] for j in xrange(m + 1)]

    for i in xrange(1, m+1):
        pointer[i][0] = Directions.up
        if not local:
            scores[i][0] = -(open_gap + (i - 1) * extend_gap)
    for i in xrange(1, n+1):
        pointer[0][i] = Directions.left
        if not local:
            scores[0][i] = -(open_gap + (i - 1) * extend_gap)

    for j in xrange(1, n + 1):
        for i in xrange(1, m + 1):
            f[i][j] = scores[i - 1][j - 1] + scoring[(s[i - 1], t[j - 1])]
            if i == 1:
                ii[i][j] = scores[i - 1][j] - open_gap
            else:
                ii[i][j] = max(ii[i - 1][j] - extend_gap, scores[i - 1][j] - open_gap)
            if j == 1:
                ij[i][j] = scores[i][j - 1] - open_gap
            else:
                ij[i][j] = max(ij[i][j - 1] - extend_gap, scores[i][j - 1] - open_gap)
            scores[i][j] = max(f[i][j], ii[i][j], ij[i][j])
            if scores[i][j] == f[i][j]:
                pointer[i][j] = Directions.diag
            elif scores[i][j] == ii[i][j]:
                pointer[i][j] = Directions.up
            elif scores[i][j] == ij[i][j]:
                pointer[i][j] = Directions.left
            else:
                raise ValueError("Uh Oh")

    return scores, pointer



def _alignment_backtrack(s, t, matrix, local=False, gap_symbol="-", start=None):
    """
    Given two strings and the levenshtein distance matrix between the two, backtrack through and create an alignment.
    Arguments: str s, str t, int[][] matrix
    Returns: str, str
    """
    s1 = ""
    t1 = ""
    if start is None:
        i = len(s)
        j = len(t)
    else:
        i = start[0]
        j = start[1]
    while i > 0 or j > 0:
        if local and matrix[i][j] == 0:
            break
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
            s1 = gap_symbol + s1
            t1 = t[j] + t1
        else:
            i -= 1
            s1 = s[i] + s1
            t1 = gap_symbol + t1
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


def _array_max_index(matrix):
    maxes = []
    indicies = []
    for i in xrange(len(matrix)):
        maxes.append(max(matrix[i]))
        indicies.append(matrix[i].index(maxes[i]))
    i = maxes.index(max(maxes))
    j = indicies[i]
    return i, j
