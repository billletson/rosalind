import rosalind

f = open("sample_inputs/PROT.txt","rb")
rna = rosalind.RNA("",f.next().strip())
proteins = rna.to_protein()

print proteins[0].sequence