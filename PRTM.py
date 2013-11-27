import rosalind

f = open("sample_inputs/PRTM.txt","rb")
protein = rosalind.Protein("",f.next().strip())

print protein.mass()