import rosalind

f = open("sample_inputs/REVP.txt","rb")
dna = rosalind.load_FASTA_file(f)[0]

for d in dna.reverse_palindromes(4,12,True):
	print " ".join([str(d[1]),str(d[2])])