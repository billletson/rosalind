import rosalind

f = open("sample_inputs/RNA.txt","rb")
dna = rosalind.DNA("", f.next().strip())
rna = dna.to_RNA()
print rna.sequence