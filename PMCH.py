import rosalind

f = open("sample_inputs/PMCH.txt","rb")

rnas = rosalind.load_FASTA_file(f, "RNA")

print rnas[0].perfect_matchings()