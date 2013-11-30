import rosalind

f = open("sample_inputs/HAMM.txt","rb")
dna1 = rosalind.DNA("", f.next().strip())
dna2 = rosalind.DNA("", f.next().strip())
print rosalind.hamming(dna1, dna2)