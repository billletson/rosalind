import rosalind

f = open("sample_inputs/GRPH.txt","rb")
dnas = rosalind.load_FASTA_file(f)
for x in rosalind.overlap_graph(dnas,3):
    print " ".join(x)
