from collections import defaultdict
from .sequences import *
import sys


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


def find_reversals(original, target):
    """
    See paper at http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.91.3123&rep=rep1&type=pdf
    """
    if original == target:
        return []
    order = [target.index(x) for x in original]
    return list(reversed(_recurse_find_reversals(order)))


def _recurse_find_reversals(order):
    breakpoints, negative_strip = _find_breakpoints(order)
    if not breakpoints:
        return []
    reversals = []
    for bp in breakpoints:
        if bp == -1:
            swaps = [order.index(0)]
        elif order[bp] == len(order) - 1:
            swaps = [order.index(order[bp] - 1)]
        elif order[bp] == 0:
            swaps = [order.index(order[bp] + 1)]
        else:
            swaps = [order.index(order[bp] - 1), order.index(order[bp] + 1)]
        for swap in swaps:
            candidate_order = order[:bp+1] + list(reversed(order[bp+1:swap+1])) + order[swap+1:]
            candidate_breakpoints, negative_strip = _find_breakpoints(candidate_order)
            reduction = len(breakpoints) - len(candidate_breakpoints)
            if len(candidate_breakpoints) == 0:
                return [(bp + 1, swap)]
            if reduction > 1 or (reduction == 1 and negative_strip):
                reversals.append(_recurse_find_reversals(candidate_order) + [(bp + 1, swap)])
    if reversals:
        return min(reversals, key=len)
    else:
        return [0] * (len(order) + 1)


def _find_breakpoints(order):
    breakpoints = []
    negative_strip = False
    if order[0] != 0:
        breakpoints.append(-1)
    for i in xrange(len(order) - 1):
        if abs(order[i] - order[i + 1]) != 1:
            breakpoints.append(i)
        elif order[i] - order[i + 1] == 1:
            negative_strip = True
    return breakpoints, negative_strip


def longest_common_subsequence(first, second):
    c = _lcs_length(first.sequence, second.sequence)
    sys.setrecursionlimit(1500)
    return DNA("%s/%s LCS" % (first.name, second.name), _lcs_backtrack(c, first.sequence, second.sequence, len(first) - 1, len(second) - 1))


def _lcs_length(first, second):
    c = [[0 for i in xrange(len(second) + 1)] for j in xrange(len(first) + 1)]
    for i in xrange(len(first)):
        for j in xrange(len(second)):
            if first[i] == second[j]:
                if i == 0 or j == 0:
                    c[i][j] = 1
                else:
                    c[i][j] = c[i-1][j-1] + 1
            else:
                if i == 0:
                    up = 0
                else:
                    up = c[i-1][j]
                if j == 0:
                    left = 0
                else:
                    left = c[i][j-1]
                c[i][j] = max(left, up)
    return c


def _lcs_backtrack(c, first, second, i, j):
    if i == -1 or j == -1:
        return ""
    elif first[i] == second[j]:
        return _lcs_backtrack(c, first, second, i-1, j-1) + first[i]
    else:
        if i == 0:
            up = 0
        else:
            up = c[i-1][j]
        if j == 0:
            left = 0
        else:
            left = c[i][j-1]
        if left > up:
            return _lcs_backtrack(c, first, second, i, j-1)
        else:
            return _lcs_backtrack(c, first, second, i-1, j)


def supersequence(first, second):
    lcs = longest_common_subsequence(first, second)
    extra_letters = _find_lcs_interleave_spots(first.sequence, lcs.sequence)
    extra_letters = _find_lcs_interleave_spots(second.sequence, lcs.sequence, extra_letters)
    result = ""
    for i in xrange(len(lcs)):
        result += extra_letters[i]
        result += lcs.sequence[i]
    result += extra_letters[-1]
    return DNA("Supersequence of %s and %s" % (first.name, second.name), result)


def _find_lcs_interleave_spots(main, lcs, extra_letters=None):
    if extra_letters is None:
        extra_letters = []
        for i in xrange(len(lcs) + 1):
            extra_letters.append("")
    idx = 0
    for s in main:
        if idx >= len(lcs):
            extra_letters[idx] += s
        elif s == lcs[idx]:
            idx += 1
        else:
            extra_letters[idx] += s
    return extra_letters


class Trie:
    def __init__(self, words):
        self.max_node = 0
        self.children = [[]]
        self.edge_labels = {}
        for word in words:
            self.insert(word)

    def insert(self, word):
        current_node = 0
        for letter in word:
            for i in self.children[current_node]:
                if letter == self.edge_labels[(current_node, i)]:
                    current_node = i
                    break
            else:
                self.max_node += 1
                self.children[current_node].append(self.max_node)
                self.children.append([])
                self.edge_labels[(current_node, self.max_node)] = letter
                current_node = self.max_node

    def edges_as_strings(self, one_index=False):
        if one_index:
            offset = 1
        else:
            offset = 0
        strings = [" ".join([str(x[0][0] + offset), str(x[0][1] + offset), x[1]]) for x in self.edge_labels.items()]
        return strings