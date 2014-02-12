import re
import urllib2
import math
import itertools
from collections import defaultdict


def load_fasta_file(f, s_type="DNA"):
    """
    Reads a FASTA formatted file (or file-like object) and creates a list of
    sequence objects
    Arguments: file f,str s_type
    Returns: Sequence[] (DNA[],RNA[],or Protein[])
    """
    sequences = []
    tmp_string = ""
    name = ""
    for line in f:
        line = line.strip()
        if line[0] == ">":
            if tmp_string:
                if s_type == "RNA":
                    sequences.append(RNA(name, tmp_string))
                elif s_type == "Protein":
                    sequences.append(Protein(name, tmp_string))
                else:
                    sequences.append(DNA(name, tmp_string))
            name = line[1:]
            tmp_string = ""
        else:
            tmp_string += line
    if tmp_string:
        if s_type == "RNA":
            sequences.append(RNA(name, tmp_string))
        elif s_type == "Protein":
            sequences.append(Protein(name, tmp_string))
        else:
            sequences.append(DNA(name, tmp_string))
    return sequences


def load_fasta_uniprot(uniprot_id, s_type="Protein"):
    """
    Given a uniprot id, creates a sequence object
    Arguments: string uniprot_id, str s_type
    Returns: Sequence (DNA,RNA,or Protein)
    """
    response = urllib2.urlopen("http://www.uniprot.org/uniprot/%s.fasta" % uniprot_id)
    response.next()
    sequence = "".join(line.strip() for line in response)
    if s_type == "RNA":
        return RNA(uniprot_id, sequence)
    elif s_type == "DNA":
        return DNA(uniprot_id, sequence)
    else:
        return Protein(uniprot_id, sequence)


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


def overlap_graph(dnas, o):
    """
    Compute overlap graph given a list of Sequence objects and a
    prefix/suffix length o. Returns list of tuples of directed connections
    Arguments: Sequence[], int
    Returns: [(str,str)]
    """
    return [(x.name, y.name) for x in dnas
                            for y in dnas
                                 if (x.name != y.name or x.sequence != y.sequence)
                                     and x.sequence[-o:] == y.sequence[:o]]


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


def n_connected_subgraphs(nodes, edges):
    """
    Given a set of nodes and undirected edges, assigns each node
    to a subgraph and returns the number of connected subgraphs
    which are not connected to each other
    Arguments: int nodes, (int,int)[] edges
    Returns: int
    """
    assignments = {}
    for n in nodes:
        assignments[n] = n
    for edge in edges:
        assign_to = assignments[edge[0]]
        assign_from = assignments[edge[1]]
        for a in assignments.items():
            if a[1] == assign_from:
                assignments[a[0]] = assign_to
    return len(set(assignments.values()))


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


def unrooted_internal_from_leaves(leaves):
    return leaves - 2


def unrooted_leaves_from_internal(internal):
    return internal + 2


def combinations(n, r):
    import operator as op
    if r == 1:
        return n
    elif n == r or r == 0:
        return 1
    else:
        return reduce(op.mul, xrange(n, n-r, -1))/reduce(op.mul, xrange(1, r+1))


def subset_count(n):
    return sum(combinations(n, r) for r in xrange(0, n+1))


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


def reversal_distance(original, target):
    reversals = 0
    for i in xrange(len(original)):
        print original
        #print target
        if original == target:
            return reversals
        if original[i] != target[i]:
            tmp = original.index(target[i])
            #print tmp, i
            #print original[0:i]
            #print list(reversed(original[i:tmp + 1]))
            #print original[tmp:]
            original = original[0:i] + list(reversed(original[i:tmp + 1])) + original[tmp + 1:]
            reversals += 1
    return reversals


class Sequence:
    codons = {}
    amino_mass = {}

    @classmethod
    def load_codons_file(cls, f):
        codons = {}
        for line in f:
            l = line.strip().split(" ")
            codons[l[0]] = l[1]
        cls.codons = codons

    @classmethod
    def load_mass_file(cls, f):
        mass = {}
        for line in f:
            l = line.strip().split()
            mass[l[0]] = float(l[1])
        cls.amino_mass = mass

    def __init__(self, name, sequence):
        """
        Parent class for DNA strands, RNA strands, and proteins
        """
        self.name = name
        self.sequence = sequence

    def count(self, sub):
        return self.sequence.count(sub)

    def __repr__(self):
        return self.name+"\n"+self.sequence

    def __len__(self):
        return len(self.sequence)

    def __eq__(self, other):
        if self.sequence == other.sequence:
            return True
        else:
            return False

    def find_substring_locations(self, sub, one_based=False):
        pattern = re.compile('(?=(%s))' % sub)
        return [x.start()+int(one_based) for x in pattern.finditer(self.sequence)]

    def find_subsequence_locations(self, sub, one_based=False):
        locations = []
        last_match = 0
        for base in sub:
            for i in xrange(last_match+1, len(self)):
                if base == self.sequence[i]:
                    locations.append(i+int(one_based))
                    last_match = i
                    break
                elif i == len(self):
                    return []
        return locations

    def substring(self, start, length):
        return self.__class__(self.name, self.sequence[start:start+length])

    def failure_array(self):
        failure = [0]
        for i in xrange(1, len(self)):
            if self.sequence[i - failure[-1]:i + 1] == self.sequence[0:failure[-1] + 1]:
                failure.append(failure[-1] + 1)
            elif self.sequence[i - failure[failure[-1] - 1]:i + 1] == self.sequence[0:failure[failure[-1] - 1] + 1]:
                failure.append(failure[failure[-1] - 1] + 1)
            else:
                failure.append(0)
        return failure


class DNA(Sequence):
    def to_rna(self):
        """
        Convert DNA sequence to RNA sequence (T to U)
        Returns RNA object
        """
        return RNA(self.name, str.replace(self.sequence, "T", "U"))

    def __complement(self):
        """
        Complements a stand and returns a string
        """
        pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        return "".join(pairs[x] for x in self.sequence)

    def complement(self):
        """
        Complements a strand and returns a DNA object
        """
        return DNA(self.name, self.__complement())

    def reverse_complement(self):
        """
        Reverse complements a strand and returns a DNA object
        """
        return DNA(self.name, self.__complement()[::-1])

    def reading_frames(self):
        """
        Return three sequences corresponding to the
        three possible reading frames
        """
        return [DNA(self.name, self.sequence[:]),
                DNA(self.name, self.sequence[1:-2]),
                DNA(self.name, self.sequence[2:-1])]

    def reverse_palindromes(self, shortest, longest, one_based=False):
        """
        Returns all reverse palindromes of length l
        min<=l<=max in a list of tuples in format
        (dna object, start loc, length)
        """
        if longest > len(self):
            longest = len(self)
        palindromes = []
        for i in xrange(shortest, longest):
            for j in xrange(len(self)-i+1):
                substring = self.substring(j, i)
                if substring == substring.reverse_complement():
                    palindromes.append((substring, j + int(one_based), i))
        return palindromes

    def log_probability(self, gc=None):
        if gc is None:
            gc = (self.count("G")+self.count("C"))/float(2*len(self))
        else:
            gc /= 2
        prob = 0
        for x in self.sequence:
            if x == "C" or x == "G":
                prob += math.log10(gc)
            else:
                prob += math.log10(0.5-gc)
        return prob

    def probability_with_repeated_attempts(self, attempts, gc=None):
        return 1 - (1.0-10.0**self.log_probability(gc))**attempts

    def k_mer_composition(self, k):
        counts = defaultdict(int)
        for i in xrange(len(self) - k + 1):
            counts[self.sequence[i:i + k]] += 1
        return [counts["".join(x)] for x in itertools.product('ACGT', repeat=k)]


class RNA(Sequence):
    pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}

    def to_dna(self):
        """
        Convert RNA sequence to DNA sequence (U to T)
        Returns DNA object
        """
        return DNA(self.name, str.replace(self.sequence, "U", "T"))

    def to_proteins(self):
        """
        Convert RNA sequnce to a list of protein sequences
        Returns list of protein objects
        """
        proteins = []
        start_flag = False
        tmp_sequence = ""
        for i in xrange(len(self)/3):
            segment = self.sequence[i*3:i*3+3]
            amino = self.codons[segment]
            if start_flag:
                if amino == "Stop":
                    proteins.append(Protein(self.name, tmp_sequence))
                    tmp_sequence = ""
                    start_flag = False
                else:
                    tmp_sequence += amino
            else:
                if amino == "M":
                    tmp_sequence += amino
                    start_flag = True
        subproteins = []
        for p in proteins:
            subproteins += p.subproteins()
        return proteins + subproteins

    def remove_intron(self, intron):
        return RNA(self.name, "".join(self.sequence.split(intron)))

    def perfect_matchings(self):
        cnt_a = self.count("A")
        cnt_u = self.count("U")
        cnt_g = self.count("G")
        cnt_c = self.count("C")

        if cnt_a != cnt_u or cnt_g != cnt_c:
            return 0
        else:
            return math.factorial(cnt_a)*math.factorial(cnt_g)

    def maximum_matchings(self):
        import operator as op
        cnt_a = self.count("A")
        cnt_u = self.count("U")
        cnt_g = self.count("G")
        cnt_c = self.count("C")
        return reduce(op.mul, xrange(max(cnt_a, cnt_u), abs(cnt_a-cnt_u), -1)) * \
               reduce(op.mul, xrange(max(cnt_g, cnt_c), abs(cnt_g-cnt_c), -1))

    def perfect_noncrossing_matchings(self):
        class Crossings(dict):
            def set_sequence(self, sequence):
                self.sequence = sequence

            def __missing__(self, key):
                start = key[0]
                end = key[1]
                self[key] = 0
                if start > end:
                    self[key] = 1
                    return 1
                if (end - start) % 2 == 0:
                    return 0
                for x in xrange(start + 1, end+1):
                    if self.sequence[start] == "A" and self.sequence[x] == "U" or \
                       self.sequence[start] == "U" and self.sequence[x] == "A" or \
                       self.sequence[start] == "C" and self.sequence[x] == "G" or \
                       self.sequence[start] == "G" and self.sequence[x] == "C":
                        self[key] += (self[(start + 1, x - 1)] * self[(x + 1, end)])
                self[key]
                return self[key]

        cross = Crossings()
        cross.set_sequence(self.sequence)
        cross[(0, len(self.sequence) - 1)]
        return cross[(0, len(self.sequence) - 1)]


class Protein(Sequence):
    def infer_rna(self, mod=0):
        """
        For a given protein, find the number of possible
        RNA stands that could produce it, mod some number
        """
        aminos = defaultdict(int)
        for c in self.codons.values():
            aminos[c] += 1
        total = 1
        for a in self.sequence:
            total *= aminos[a]
            if mod > 0:
                total %= mod
        total *= aminos['Stop']
        if mod > 0:
            total %= mod
        return total

    def mass(self):
        """
        Return the mass of a protein, using
        a reference of amino acid masses
        """
        mass = 0
        for a in self.sequence:
            mass += self.amino_mass[a]
        return mass

    def subproteins(self):
        """
        Find all of the valid proteins that are a subsequence
        of this protein
        """
        return [Protein(self.name, self.sequence[x:]) for x in xrange(len(self)) if self.sequence[x] == "M"]

try:		
    Sequence.load_codons_file(open("reference_data/codons.txt"))
except:
    pass
try:
    Sequence.load_mass_file(open("reference_data/mass.txt"))
except:
    pass
