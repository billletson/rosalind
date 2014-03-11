import itertools
import math


def punnett_probabilities(k, m, n):
    """
    Takes count(or prob) of homogeneous dominant, heterogeneous,
    and homogeneous recessive organisms and returns
    probabilities of same from random mating, in same order
    Arguments: int/float k, int/float m, int/float n
    Returns: float,float,float
    """
    pop = k+m+n
    # Probabilities of each, considering non-replacement
    pk = float(k)/pop
    pkp1 = float(k)/(pop-1)
    pk1p1 = float(k-1)/(pop-1)
    pm = float(m)/pop
    pmp1 = float(m)/(pop-1)
    pm1p1 = float(m-1)/(pop-1)
    pn = float(n)/pop
    pnp1 = float(n)/(pop-1)
    pn1p1 = float(n-1)/(pop-1)

    hom_dom = pk*pk1p1 + pk*pmp1*0.5 + pm*pkp1*0.5 + pm*pm1p1*0.25
    het = pk*pmp1*0.5 + pk*pnp1 + pm*pkp1*0.5 + pm*pm1p1*0.5 + pm*pnp1*0.5 + pn*pkp1 + pn*pmp1*0.5
    hom_rec = pm*pm1p1*0.25 + pm*pnp1*0.5 + pn*pmp1*0.5 + pn*pn1p1

    return hom_dom, het, hom_rec


def partial_permutation(n, k, mod=0):
    """
    UNKNOWN
    """
    total = 1
    for x in xrange(n, n-k, -1):
        total *= x
        if mod != 0:
            total %= mod
    return total


def signed_permutations(n):
    """
    Enumerates the signed permutations of the integers 1 to n
    The signed permutations are all permutations of the base
    integers with all permutations of sign
    Arguments: int n
    Returns: int[][]
    """
    pos_perms = itertools.permutations(xrange(1, n+1), n)
    signs = itertools.product([1, -1], repeat=n)
    signed_perms_uncombined = itertools.product(pos_perms, signs)
    return [[x*y for x, y in zip(s[0], s[1])] for s in signed_perms_uncombined]


def combinations(n, r):
    import operator as op
    if r == 1:
        return n
    elif n == r or r == 0:
        return 1
    else:
        return reduce(op.mul, xrange(n, n-r, -1))/reduce(op.mul, xrange(1, r+1))


def subset_count(n, min_r=0, max_r=None):
    if max_r is None:
        max_r = n + 1
    return sum(combinations(n, r) for r in xrange(min_r, max_r))


def lexicographic_permutations(letters, n):
    """
    Enumerate the permutations, with replacement, of length n of a list, in lexicographic order
    Arguments: str[] letters, int n
    Returns: str[][]
    """
    perms = list(itertools.product(letters, repeat=n))
    perms.sort(key=lambda x: [letters.index(y) for y in x])
    return ["".join(x) for x in perms]


def multilength_lexicographic_permutations(letters, n):
    perms = []
    for i in xrange(1, n + 1):
        perms += lexicographic_permutations(letters, i)
    return sorted(perms, key=lambda word: [letters.index(c) for c in word])


def expected_restriction_sites(str_len, substring, gc_content):
    same_length_prob = 1
    for l in substring:
        if l in ["G", "C"]:
            same_length_prob *= gc_content / 2
        else:
            same_length_prob *= (1 - gc_content) / 2
    return same_length_prob * (str_len - len(substring) + 1)


def binomial_success_probability(n, k, p, log=False):
    if log:
        return math.log10(combinations(n, k)) + k * math.log10(p) + (n - k) * math.log10(1 - p)
    else:
        return combinations(n, k) * p ** k * p ** (n - k)


def binomial_cdf(n, p, log=False):
    probs = [binomial_success_probability(n, i, p, log) for i in xrange(n + 1)]
    cdf = [probs[0]]
    for i in xrange(1, len(probs)):
        if log:
            cdf.append(logsumexp([probs[i], cdf[-1]]))
        else:
            cdf.append(probs[i] + cdf[-1])
    return cdf


def logsumexp(x):
    a = max(x)
    return a + math.log10(sum([10**(y - a) for y in x]))




