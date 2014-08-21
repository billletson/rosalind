from .io import load_scoring_matrix


class Directions:
    undefined, diag, up, left = range(4)

    def __init__(self):
        pass


def levenshtein(first, second, scoring_matrix=None, gap=1):
    """
    Calculate the levenshtein distance between two sequences. http://en.wikipedia.org/wiki/Levenshtein_distance
    Wraps alignment score, with a negation
    Arguments: Sequence first, Sequence second, str scoring_matrix, int or (int, int) gap
    Returns: int
    """
    return -1 * alignment_score(first, second, scoring_matrix, gap)


def alignment_score(first, second, scoring_matrix=None, gap=1):
    """
    Calculate the alignment score distance between two sequences.
    Can take a custom scoring matrix and gap penalty
    Arguments: Sequence first, Sequence second, str scoring_matrix, int or (int, int) gap
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
    Arguments: Sequence first, Sequence second, str scoring_matrix, int or (int, int) gap
    Returns: int, str, str
    """
    if scoring_matrix is None:
        scoring = _DefaultMatrix()
    else:
        scoring = load_scoring_matrix(scoring_matrix)
    s = first.sequence
    t = second.sequence
    if type(gap) == int:
        matrix = _alignment_matrix(s, t, scoring, gap)
        s1, t1 = _alignment_backtrack(s, t, matrix, scoring, gap)
    else:
        matrix, pointer = _affine_alignment_matrix(s, t, scoring, gap)
        s1, t1, = _alignment_backtrack_pointer(s, t, matrix, pointer)
    return -1*matrix[-1][-1], s1, t1


def optimal_alignment_count(first, second, modulus=None, scoring_matrix=None, gap=1):
    """
    Finds the number of valid optimal alignments for a pair of sequences. As this can get large, can return modulo
    some number.
    Arguments: Sequence first, Sequence second, int modulus, str scoring_matrix, int or (int, int) gap
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
    Arguments: Sequence first, Sequence second, str scoring_matrix, int or (int, int) gap
    Returns: ?
    """
    if scoring_matrix is None:
        scoring = _DefaultMatrix()
    else:
        scoring = load_scoring_matrix(scoring_matrix)
    s = first.sequence
    t = second.sequence
    if type(gap) is int:
        matrix = _alignment_matrix(s, t, scoring, gap, True)
        start = _array_max_index(matrix)
        s1, t1 = _alignment_backtrack(s, t, matrix, scoring, gap, True, "", start)
    else:
        matrix, pointer = _affine_alignment_matrix(s, t, scoring, gap, True)
        start = _array_max_index(matrix)
        s1, t1 = _alignment_backtrack_pointer(s, t, matrix, pointer, True, "", start)
    return matrix[start[0]][start[1]], s1, t1


def _alignment_matrix(s, t, scoring, gap, local=False):
    """
    Calculate a matrix showing the alignment score between all prefixes of two strings, with the last element
    being the alignment score between the whole two strings. Only works for a linear gap penalty, but uses less space
    than a more general solution.
    Arguments: str s, str t, {(str, str): int} scoring, int gap, bool local
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
    """
    Calculate a matrix showing the alignment score between all prefixes of two strings, with the last element
    being the alignment score between the whole two strings. Works for affine gap penalty (constant gap penalty can
    be expressed as an affine penalty with 0 extension penalty, so this works for constant penalty too.) Uses 5 times
    the space as the linear version, and needs to return a pointer matrix.
    Arguments: str s, str t, {(str, str): int} scoring, int gap, bool local
    Returns: int[][], int[][]
    """
    m = len(s)
    n = len(t)
    open_gap, extend_gap = gap
    scores = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    f = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    ii = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    ij = [[0 for i in xrange(n + 1)] for j in xrange(m + 1)]
    pointer = {'f': [[Directions.undefined for i in xrange(n + 1)] for j in xrange(m + 1)],
               'i': [[Directions.up for i in xrange(n + 1)] for j in xrange(m + 1)],
               'j': [[Directions.left for i in xrange(n + 1)] for j in xrange(m + 1)]}

    for i in xrange(1, m+1):
        pointer['f'][i][0] = Directions.up
        pointer['i'][i][0] = Directions.up
        pointer['j'][i][0] = Directions.up
        if not local:
            scores[i][0] = -(open_gap + (i - 1) * extend_gap)
            ii[i][0] = -(open_gap + (i - 1) * extend_gap)
    for i in xrange(1, n+1):
        pointer['f'][0][i] = Directions.left
        pointer['i'][0][i] = Directions.left
        pointer['j'][0][i] = Directions.left
        if not local:
            scores[0][i] = -(open_gap + (i - 1) * extend_gap)
            ij[0][i] = -(open_gap + (i - 1) * extend_gap)

    for j in xrange(1, n + 1):
        for i in xrange(1, m + 1):
            f[i][j] = scores[i - 1][j - 1] + scoring[(s[i - 1], t[j - 1])]
            if local and f[i][j] < 0:
                f[i][j] = 0
            if scores[i - 1][j - 1] == f[i - 1][j - 1]:
                pointer['f'][i][j] = Directions.diag
            elif scores[i - 1][j - 1] == ii[i - 1][j - 1]:
                pointer['f'][i][j] = Directions.up
            elif scores[i - 1][j - 1] == ij[i - 1][j - 1]:
                pointer['f'][i][j] = Directions.left
            else:
                raise ValueError("Uh Oh")
            if i == 1:
                ii[i][j] = scores[i - 1][j] - open_gap
                if local and ii[i][j] < 0:
                    ii[i][j] = 0
            else:
                ii[i][j] = max(ii[i - 1][j] - extend_gap, scores[i - 1][j] - open_gap)
                if local and ii[i][j] < 0:
                    ii[i][j] = 0
                if ii[i][j] == scores[i - 1][j] - open_gap and scores[i - 1][j] == f[i - 1][j]:
                    pointer['i'][i][j] = Directions.diag
            if j == 1:
                ij[i][j] = scores[i][j] - open_gap
                if local and ij[i][j] < 0:
                    ij[i][j] = 0
            else:
                ij[i][j] = max(ij[i][j - 1] - extend_gap, scores[i][j - 1] - open_gap)
                if local and ij[i][j] < 0:
                    ij[i][j] = 0
                if ij[i][j] == scores[i][j - 1] - open_gap and scores[i][j - 1] == f[i][j - 1]:
                    pointer['j'][i][j] = Directions.diag

            scores[i][j] = max(f[i][j], ii[i][j], ij[i][j])
    if scores[-1][-1] == f[-1][-1]:
        pointer['start'] = 'f'
    elif scores[-1][-1] == ii[-1][-1]:
        pointer['start'] = 'i'
    elif scores[-1][-1] == ij[-1][-1]:
        pointer['start'] = 'j'
    return scores, pointer


def _alignment_backtrack(s, t, matrix, scoring, gap, local=False, gap_symbol="-", start=None):
    """
    Given two strings and the alignment score matrix between the two, backtrack through and create an alignment.
    Arguments: str s, str t, int[][] matrix, {(str, str): int} scoring, int gap, bool local, str gap_symbol, (int, int) start
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
        if j > 0 and i > 0 and matrix[i][j] - matrix[i - 1][j - 1] == scoring[(s[i - 1], t[j - 1])]:
            i -= 1
            j -= 1
            s1 = s[i] + s1
            t1 = t[j] + t1
        elif i > 0 and matrix[i][j] - matrix[i - 1][j] == -gap:
            i -= 1
            s1 = s[i] + s1
            t1 = gap_symbol + t1
        elif j > 0 and matrix[i][j] - matrix[i][j - 1] == -gap:
            j -= 1
            s1 = gap_symbol + s1
            t1 = t[j] + t1
        else:
            raise ValueError("Wrong scoring matrix?")
    return s1, t1


def _alignment_backtrack_pointer(s, t, matrix, pointer, local=False, gap_symbol="-", start=None):
    """
    Given two strings and a set of pointer matricies, backtrack through and create an alignment. The alignment score matrix is
    only really necessary for the local version, but will require a refactor to make optional.
    Arguments: str s, str t, int[][] matrix, int[][] pointer, bool local, str gap_symbol, (int, int) start
    Returns: str, str
    """
    s1 = ""
    t1 = ""
    if start is None:
        i = len(s)
        j = len(t)
        pointer_id = pointer['start']
    else:
        i = start[0]
        j = start[1]
        pointer_id = 'f'
    current_pointer = pointer[pointer_id]
    while i > 0 or j > 0:
        next_dir = current_pointer[i][j]
        if local and matrix[i][j] == 0:
            break
        if pointer_id == 'f':
            i -= 1
            j -= 1
            s1 = s[i] + s1
            t1 = t[j] + t1
        elif pointer_id == 'i':
            i -= 1
            s1 = s[i] + s1
            t1 = gap_symbol + t1
        elif pointer_id == 'j':
            j -= 1
            t1 = t[j] + t1
            s1 = gap_symbol + s1
        else:
            raise ValueError("Invalid pointer matrix")

        if next_dir == Directions.diag:
            pointer_id = 'f'
        elif next_dir == Directions.up:
            pointer_id = 'i'
        elif next_dir == Directions.left:
            pointer_id = 'j'
        else:
            raise ValueError("Invalid pointer matrix")
        current_pointer = pointer[pointer_id]
    return s1, t1


def _count_alignments(s, t, matrix, modulus=None):
    """
    Using an alignment score matrix and the two original strings, count the number of optimal alignments modulo some
    number. This currently does not work for an arbitrary scoring scheme, only the simplest default.
    Arguments: str s, str t, int[][] matrix, int or None modulus
    Returns: int
    """
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
    """
    Dict-like object for the default scoring of -1 for a substitution, 0 for a match.
    """
    def __missing__(self, key):
        if key[0] == key[1]:
            return 0
        else:
            return -1


def _array_max_index(matrix):
    """
    Calculates the index of the maximum of a list of lists of ints.
    Arguments: int[][]
    Returns: int, int
    """
    maxes = []
    indicies = []
    for i in xrange(len(matrix)):
        maxes.append(max(matrix[i]))
        indicies.append(matrix[i].index(maxes[i]))
    i = maxes.index(max(maxes))
    j = indicies[i]
    return i, j
