import sys
from .sequences import *


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

