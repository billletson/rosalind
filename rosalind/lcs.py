from .sequences import *
import sys

def longest_common_subsequence(first, second):
    """
    Find the longest common subsequence of two strands.
    Arguments: DNA first, DNA second
    Returns: DNA
    """
    c = _lcs_length(first.sequence, second.sequence)
    return DNA("%s/%s LCS" % (first.name, second.name), _lcs_backtrack(c, first.sequence, second.sequence))


def _lcs_length(first, second):
    """
    Create a matrix with the lengths of the lcs for substrings of the two strings passed in, iteratively building
    from the first characters of each.
    Arguments: str first, str second
    Returns: int[][]
    """
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


def _lcs_backtrack(c, first, second):
    """
    Given a matrix of lcs lengths of substrings, walk back through and find the characters making up the lcs of the
    whole strings.
    Arguments: int[][] c, str first, str second
    Returns: str
    """
    i = len(first) - 1
    j = len(second) - 1
    lcs = ""
    while i >= 0 or j >= 0:
        if first[i] == second[j]:
            lcs = first[i] + lcs
            if i == 0 or j == 0:
                break
            else:
                i -= 1
                j -= 1
        else:
            if i == 0:
                j -= 1
            elif j == 0:
                i -= 1
            elif c[i][j-1] > c[i-1][j]:
                j -= 1
            else:
                i -= 1
    return lcs


def supersequence(first, second):
    """
    Find the shortest common supersequence of two strands.
    Arguments: DNA first, DNA second
    Returns: DNA
    """
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
    """
    Find the letters in the target string that are not also the longest common subsequence. Slot them into buckets
    representing the spaces in between the characters of the lcs. Can pass in an existing set of extra letters.
    Arguments: str main, str lcs, str[] extra_letters
    Return: str[]
    """
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


def interweaveable_matrix(haystack, needles):
    """
    For an array of needles, check if each pair is interweaveable into a haystack.
    Arguments: Sequence or str haystack, Sequence[] or str[] needles
    Returns: int[][]
    """
    matrix = [[0 for x in xrange(len(needles))] for i in xrange(len(needles))]
    for i in xrange(len(needles)):
        for j in xrange(len(needles)):
            if interweaveable(haystack, needles[i], needles[j]):
                matrix[i][j] = 1
    return matrix


def interweaveable(haystack, first, second):
    """
    Check if two needle strands can be interwoven into a haystack stand. Can be passed in as strings or Sequence
    objects, as long as all three are of the same type.
    Arguments: Sequence or str haystack, Sequence or str first, Sequence or str second
    Returns: bool
    """
    try:
        haystack = haystack.sequence
        first = first.sequence
        second = second.sequence
    except AttributeError:
        pass
    total_length = len(first) + len(second)
    for i in xrange(len(haystack) - total_length + 1):
        if _interweaveable_recursion(haystack[i:i + total_length], first, second):
            return True
    else:
        return False


def _interweaveable_recursion(haystack, first, second):
    """
    Given a haystack string and two needle strings (with length haystack equal to the sum of the lengths of the
    needles), determine if the needles can be interwoven into the haystack as disjoint subsequences.
    Arguments: str haystack, str first, str second
    Returns: bool
    """
    if not first:
        if haystack == second:
            return True
        else:
            return False
    if not second:
        if haystack == first:
            return True
        else:
            return False
    h = 0
    i = 0
    j = 0
    l_h = len(haystack)
    l_i = len(first)
    l_j = len(second)
    if l_h != l_i + l_j:
        raise Exception("OH SHIT")
    while i < l_i or j < l_j:
        if i < l_i and j < l_j and haystack[h] == first[i] and haystack[h] == second[j]:
            if _interweaveable_recursion(haystack[h + 1:], first[i + 1:], second[j:]) or \
               _interweaveable_recursion(haystack[h + 1:], first[i:], second[j + 1:]):
                return True
            else:
                return False
        elif i < l_i and haystack[h] == first[i]:
            h += 1
            i += 1
        elif j < l_j and haystack[h] == second[j]:
            h += 1
            j += 1
        else:
            return False
    return True



