import rosalind

f = open("sample_inputs/TRAN.txt","rb")
dnas = rosalind.load_FASTA_file(f)

transitions, transversions = rosalind.transition_transversion(dnas[0], dnas[1])

print float(transitions)/transversions

