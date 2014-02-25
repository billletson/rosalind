import math
import re
from collections import defaultdict
import itertools
import os


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
                return self[key]

        cross = Crossings()
        cross.set_sequence(self.sequence)
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

    @property
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
    Sequence.load_codons_file(open(os.path.dirname(__file__) + "/reference_data/codons.txt"))
except:
    pass
try:
    Sequence.load_mass_file(open(os.path.dirname(__file__) + "/reference_data/mass.txt"))
except:
    pass


def infer_protein_from_prefix_spectrum(weights):
    diff = []
    for i in xrange(len(weights) - 1):
        diff.append(weights[i + 1] - weights[i])
    aminos = []
    for d in diff:
        for a in Sequence.amino_mass.items():
            if abs(d-a[1])<0.01:
                aminos.append(a[0])
                break
    return Protein("", "".join(aminos))
