import rosalind

f = open("sample_inputs/SSEQ.txt","rb")

dnas = rosalind.load_FASTA_file(f)

print " ".join(map(str, dnas[0].find_subsequence_locations(dnas[1].sequence, True)))


