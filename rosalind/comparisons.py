from collections import defaultdict
from .sequences import *


def hamming(first, second):
    """
    Computes the hamming distance between two sequences
    Arguments: Sequence first, Sequence second
    Returns: int
    """
    return _string_hamming(first.sequence, second.sequence)


def _string_hamming(first, second):
    count = 0
    for x in zip(first, second):
        if x[0] != x[1]:
            count += 1
    return count


def levenshtein(first, second):
    """
    Calculate the levenshtein distance between two sequences. http://en.wikipedia.org/wiki/Levenshtein_distance
    Arguments: Sequence first, Sequence second
    Returns: int
    """
    s = first.sequence
    t = second.sequence
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
    return matrix[m][n]


def transition_transversion(first, second):
    """
    Computes the transition and transversion counts for two DNA sequences
    Arguments: DNA first, DNA second
    Returns: int,int
    """
    transitions = 0
    transversions = 0
    for x, y in zip(first.sequence, second.sequence):
        if x != y:
            if (x in "AG" and y in "AG") or (x in "CT" and y in "CT"):
                transitions += 1
            else:
                transversions += 1
    return transitions, transversions


def consensus(dnas):
    """
    Compute consensus matrix and consensus string for a group of DNA dequences
    Arguments: Sequence[] dnas
    Returns: [{str:int}],str
    REFACTOR
    """
    matrix = []
    for i in xrange(len(dnas[0])):
        column = defaultdict(int)
        for dna in dnas:
            column[dna.sequence[i]] += 1
        matrix.append(column)
    string = "".join(max(c, key=c.get) for c in matrix)
    return matrix, string


def difference_matrix(dnas):
    """
    Construct a difference matrix for a list of sequences. A difference matrix is a matrix of hamming distances/total
    characters for each pair of sequences.
    Arguments: Sequence[]
    Returns: float[[]]
    """
    return [[hamming(x, y) / float(len(x)) for y in dnas] for x in dnas]


def longest_common_substring(dnas):
    """
    Find the longest common substring in a group of Sequences
    Arguments: Sequence[] dnas
    Returns: str
    """
    ordered = sorted([x.sequence for x in dnas])
    short = ordered[0]
    longer = ordered[1:]
    common_substrings = []
    for i in xrange(1, len(short)+1):
        new_substrings = []
        for j in xrange(len(short)-i+1):
            test = short[j:j+i]
            test_found = True
            for l in longer:
                if test not in l:
                    test_found = False
                    break
            if test_found:
                new_substrings.append(test)
        if len(new_substrings) == 0:
            break
        else:
            common_substrings += new_substrings
    return max(common_substrings, key=len)


def longest_monotonic_subsequence(series, decreasing=False):
    """
    Find the longest monotonically increasing or decreasing subsequence of a
    series, recalling that a subsequence is not necessarily contiguous
    Arguments: int[] series, bool decreasing
    Returns: int[]
    """
    subs = [[series[0]]]
    if decreasing:
        series = [-1*x for x in series]
    for n in series:
        for i in xrange(len(subs)-1, -1, -1):
            if n > subs[i][-1]:
                if i == len(subs)-1:
                    subs.append(subs[i]+[n])
                elif n < subs[i+1][-1]:
                    subs[i+1] = subs[i]+[n]
            elif n < subs[i][-1] and (len(subs[i]) < 2 or n > subs[i][-2]):
                subs[i][-1] = n
    if decreasing:
        return [-1*x for x in subs[-1]]
    else:
        return subs[-1]


def absorb_sequence(master, absorbed):
    """
    Given two sequences, combines into one sequence if overlap is longer
    than 1/2 the length of the second sequence. If no sufficiently long
    overlap exists, return None
    Arguments: DNA master, DNA absorbed
    Returns: DNA or None
    """
    l = len(absorbed)
    for i in xrange(l, l/2, -1):
        if master.sequence[:i] == absorbed.sequence[-i:]:
            return DNA(absorbed.name+"/"+master.name, absorbed.sequence[:l-i]+master.sequence)
        if master.sequence[-i:] == absorbed.sequence[:i]:
            return DNA(master.name+"/"+absorbed.name, master.sequence+absorbed.sequence[-1*(l-i):])
    return None


def superstring(dnas):
    """
    Given a list of dna sequences, combine into one superstring based on
    overlaps of at least 1/2 length. Return both superstring, and any leftover
    sequences which do not overlap.
    Arguments: DNA[] dnas
    Returns: DNA, DNA[]
    """
    master = dnas[0]
    dnas = dnas[1:]
    merging = True
    while merging:
        survivors = []
        merging = False
        for dna in dnas:
            new = absorb_sequence(master, dna)
            if new is None:
                survivors.append(dna)
            else:
                master = new
                merging = True
        dnas = survivors[:]
    return master, dnas


def identify_read_errors(dnas):
    """
    Given a list of sequences, find any that only exist once, as either themselves or a reverse complement. Return both
    the incorrect sequence and a corrected version which is in the list twice and is a hamming distance of 1 away.
    Arguments: DNA[] dnas
    Return: (DNA, DNA)[]
    """
    corrections = []
    counts = defaultdict(int)
    for dna in dnas:
        counts[dna.sequence] += 1
        counts[dna.reverse_complement().sequence] += 1
    for dna in dnas:
        if counts[dna.sequence] < 2:
            for correct in counts.keys():
                if counts[correct] > 1:
                    if _string_hamming(dna.sequence, correct) == 1:
                        corrections.append((dna, DNA(dna.name, correct)))
                        break
    return corrections


def minkowski_difference(first, second, sort=True, precision=5):
    """
    Find the Minkowski difference between two multisets. Minkowksi difference is the multiset of all the differences of
    combinations of one element of the first set and one element of the second set.
    Arguments: float[] first, float[] second, bool sort, int precision
    Return: float[]
    """
    mink = [round(x - y, precision) for x in first for y in second]
    if sort:
        mink = sorted(mink, key=mink.count, reverse=True)
    return mink
