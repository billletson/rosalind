import rosalind

f = open("sample_inputs/MPRT.txt","rb")
sub = 'N[^P][ST][^P]'
for line in f:
    protein = rosalind.load_FASTA_uniprot(line.strip())
    matches = [str(x) for x in protein.find_substring_locations(sub,True)]
    if len(matches)>0:
        print protein.name
        print " ".join(matches)