import rosalind

f = open("sample_inputs/MRNA.txt","rb")
protein = rosalind.Protein("",f.next().strip())

print protein.infer_RNA(100000)