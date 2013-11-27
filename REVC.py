import rosalind

f = open("sample_inputs/REVC.txt","rb")
dna = rosalind.DNA("",f.next().strip())
comp = dna.reverse_complement()
print comp.sequence
