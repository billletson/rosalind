import rosalind

dnas = rosalind.load_FASTA_file(open("sample_inputs/CONS.txt","rb"))

matrix,string = rosalind.consensus(dnas)

print string
print "A: " + " ".join([str(x['A']) for x in matrix])
print "C: " + " ".join([str(x['C']) for x in matrix])
print "G: " + " ".join([str(x['G']) for x in matrix])
print "T: " + " ".join([str(x['T']) for x in matrix])