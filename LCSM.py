import rosalind

f = open("sample_inputs/LCSM.txt","rb")

dnas = rosalind.load_FASTA_file(f)

print rosalind.longest_common_substring(dnas)
