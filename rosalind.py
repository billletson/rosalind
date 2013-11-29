import re
import urllib2
import math
import itertools
from collections import defaultdict


def load_FASTA_file(f,type="DNA"):
	"""
	Takes a file-like object as input, returns a list of sequence objects
	"""
	sequences = []
	tmp_string = ""
	for line in f:
		line = line.strip()
		if(line[0]==">"):
			if(tmp_string):
				if(type=="RNA"):
					sequences.append(RNA(name,tmp_string))
				elif(type=="Protein"):
					sequences.append(Protein(name,tmp_string))
				else:
					sequences.append(DNA(name,tmp_string))
			name = line[1:]
			tmp_string = ""
		else:
			tmp_string+=line
	if(tmp_string):
		if(type=="RNA"):
			sequences.append(RNA(name,tmp_string))
		elif(type=="Protein"):
			sequences.append(Protein(name,tmp_string))
		else:
			sequences.append(DNA(name,tmp_string))
	return sequences

def load_FASTA_uniprot(id,type="Protein"):
	base_url = "http://www.uniprot.org/uniprot/%s.fasta"
	response = urllib2.urlopen(base_url % (id))
	response.next()
	sequence = ""
	for line in response:
		sequence += line.strip()
	if(type=="RNA"):
		return RNA(id,sequence)
	elif(type=="DNA"):
		return DNA(id,sequence)
	else:
		return Protein(id,sequence)
	
def hamming(first,second):
	"""
	Computes the hamming distance between two sequences
	"""
	count = 0
	for x in zip(first.sequence,second.sequence):
		if(x[0]!=x[1]):
			count+=1
	return count
	
def transition_transversion(first,second):
	"""
	Computes the transition and transversion counts for two DNA sequences
	"""
	
	transitions = 0
	transversions = 0
	for x,y in zip(first.sequence,second.sequence):
		if(x!=y):
			if((x in "AG" and y in "AG") or (x in "CT" and y in "CT")):
				transitions += 1
			else:
				transversions += 1
	return (transitions, transversions)
		
def consensus(dna_list):
	"""
	Compute consensus matrix and consensus string
	"""
	matrix = []
	for i in xrange(len(dna_list[0])):
		column = defaultdict(int)
		for dna in dna_list:
			column[dna.sequence[i]]+=1
		matrix.append(column)
	string = ""
	for c in matrix:
		string += max(c,key=c.get)
	return [matrix,string]
	
def overlap_graph(dna_list,o):
	"""
	Compute overlap graph given a list of dna objects and a 
	prefix/suffix length o
	Returns list of tuples of directed connections
	"""
	pairs = []
	for x in dna_list:
		for y in dna_list:
			if(x.name!=y.name or x.sequence!=y.sequence):	
				if(x.sequence[-o:]==y.sequence[:o]):
					pairs.append((x.name,y.name))
	return pairs
	
def longest_common_substring(dna_list):
	ordered = sorted([x.sequence for x in dna_list])
	short = ordered[0]
	long = ordered[1:]
	common_substrings = []
	for i in xrange(1,len(short)+1):
		new_substrings = []
		for j in xrange(len(short)-i+1):
			test = short[j:j+i]
			test_found = True
			for l in long:
				if(test not in l):
					test_found = False
					break
			if(test_found):
				new_substrings.append(test)
		if(len(new_substrings)==0):
			break
		else:
			common_substrings += new_substrings
	return max(common_substrings,key=len)
	
def punnett_probabilities(k,m,n):
	"""
	Takes count(or prob) of homogeneous dominant, heterogeneous,
	and homogeneous recessive organisms and returns
	probabilities of same from random mating, in same order as
	a tuple
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
	
	return (hom_dom,het,hom_rec)
	
def partial_permutation(n,k,mod=0):
	total = 1
	for x in xrange(n,n-k,-1):
		total *= x
		if(mod!=0):
			total = total % mod
	return total

def lexographic_permutations(letters,n):
	perms = list(itertools.product(letters,repeat=n))
	perms.sort(key=lambda x:[letters.index(y) for y in x])
	return perms
	
class Sequence:
	codons = {}
	amino_mass = {}
	
	@classmethod
	def load_codons_file(cls,f):
		codons = {}
		for line in f:
			l = line.strip().split(" ")
			codons[l[0]] = l[1]
		cls.codons = codons
		
	@classmethod
	def load_mass_file(cls,f):
		mass = {}
		for line in f:
			l=line.strip().split()
			mass[l[0]] = float(l[1])
		cls.amino_mass = mass
	

	def __init__(self,name,sequence):
		"""
		Parent class for DNA strands, RNA strands, and proteins
		"""
		self.name = name
		self.sequence = sequence
		
	def count(self,sub):
		return self.sequence.count(sub)
		
	def __repr__(self):
		return self.name+"\n"+self.sequence
		
	def __len__(self):
		return len(self.sequence)
		
	def __eq__(self,other):
		if(self.sequence==other.sequence):
			return True
		else:
			return False
		
	def find_substring_locations(self,sub,one_based=False):
		pattern = re.compile('(?=(%s))' % sub)
		return [x.start()+int(one_based) for x in pattern.finditer(self.sequence)]
		
	def find_subsequence_locations(self,sub,one_based=False):
		locations = []
		last_match = 0
		for base in sub:
			for i in xrange(last_match+1,len(self)):
				if(base==self.sequence[i]):
					locations.append(i+int(one_based))
					last_match=i
					break
				elif(i==len(self)):
					return []
		return locations
			
		
	def slice(self,start,length):
		return self.__class__(self.name,self.sequence[start:start+length])
				
class DNA(Sequence):
	def to_RNA(self):
		"""
		Convert DNA sequence to RNA sequence (T to U)
		Returns RNA object
		"""
		return RNA(self.name,str.replace(self.sequence,"T","U"))
		
	def __complement(self):
		"""
		Complements a stand and returns a string
		"""
		new_strand = ""
		pairs = {"A":"T","T":"A","G":"C","C":"G"}
		for x in self.sequence:
			new_strand+=pairs[x]
		return new_strand
	
	def complement(self):
		"""
		Complements a strand and returns a DNA object
		"""
		return DNA(self.name,self.__complement())
		
	def reverse_complement(self):
		"""
		Reverse complements a strand and returns a DNA object
		"""
		return DNA(self.name,self.__complement()[::-1])
		
	def reading_frames(self):
		"""
		Return three sequences corresponding to the
		three possible reading frames
		"""
		dnas = []
		dnas.append(DNA(self.name,self.sequence[:]))
		dnas.append(DNA(self.name,self.sequence[1:-2]))
		dnas.append(DNA(self.name,self.sequence[2:-1]))
		return dnas
		
	def reverse_palindromes(self,min,max,one_based=False):
		"""
		Returns all reverse palindromes of length l
		min<=l<=max in a list of tuples in format 
		(dna object, start loc, length)
		"""
		if max>len(self):
			max=len(self)
		palindromes = []
		for i in xrange(min,max+1):
			for j in xrange(len(self)-i+1):
				slice = self.slice(j,i)
				if(slice==slice.reverse_complement()):
					palindromes.append((slice,j+int(one_based),i))
		return palindromes
		
	def log_probability(self, GC=None):
		if(GC==None):
			GC = (self.count("G")+self.count("C"))/float(2*len(self))
		else:
			GC=GC/2
		prob = 0
		for x in self.sequence:
			if(x=="C" or x=="G"):
				prob += math.log10(GC)
			else:
				prob += math.log10(0.5-GC)
		return prob
						
class RNA(Sequence):
	def to_DNA(self):
		"""
		Convert RNA sequence to DNA sequence (U to T)
		Returns DNA object
		"""
		return DNA(self.name,str.replace(self.sequence,"U","T"))
		
	def to_protein(self):
		"""
		Convert RNA sequnce to a list of protein sequences
		Returns lsit of protein objects
		"""
		proteins = []
		start_flag = False
		tmp_sequence = ""
		for i in xrange(len(self)/3):
			segment = self.sequence[i*3:i*3+3]
			amino = self.codons[segment]
			if(start_flag):
				if(amino == "Stop"):
					proteins.append(Protein(self.name,tmp_sequence))
					tmp_sequence = ""
					start_flag = False
				else:
					tmp_sequence += amino
			else:
				if(amino=="M"):
					tmp_sequence += amino
					start_flag = True
		subproteins = []
		for p in proteins:
			subproteins+= p.subproteins()
		return proteins + subproteins
		
	def remove_intron(self,intron):
		return RNA(self.name,"".join(self.sequence.split(intron)))
		
	def perfect_matchings(self):
		cnt_a = self.count("A")
		cnt_u = self.count("U")
		cnt_g = self.count("G")
		cnt_c = self.count("C")
		
		if(cnt_a!=cnt_u or cnt_g!=cnt_c):
			return 0
		else:
			return math.factorial(cnt_a)*math.factorial(cnt_g)
									
class Protein(Sequence):
	def infer_RNA(self,mod=0):
		"""
		For a given protein, find the number of possible
		RNA stands that could produce it, mod some number
		"""
		aminos = defaultdict(int)
		for c in self.codons.values():
			aminos[c]+=1
		total = 1
		for a in self.sequence:
			total *= aminos[a]
			if(mod>0):
				total = total % mod
		total *= aminos['Stop']
		if(mod>0):
			total = total % mod
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
		subproteins = []
		for x in xrange(len(self)):
			if(self.sequence[x]=="M"):
				subproteins.append(Protein(self.name,self.sequence[x:]))
		return subproteins
			
try:		
	Sequence.load_codons_file(open("reference_data/codons.txt"))
except:
	pass
try:
	Sequence.load_mass_file(open("reference_data/mass.txt"))
except:
	pass
	