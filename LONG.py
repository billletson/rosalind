import rosalind

g = open("answer.txt","wb")
dnas = rosalind.load_FASTA_file(open("test_inputs/LONG.txt","rb"))
superstring, leftovers = rosalind.superstring(dnas)
g.write(superstring.sequence)
