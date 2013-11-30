import rosalind

f = open("sample_inputs/SUBS.txt","rb")
dna = rosalind.DNA("", f.next().strip())
print " ".join([str(x) for x in dna.find_substring_locations(f.next().strip(), True)])