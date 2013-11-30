import rosalind

f = open("sample_inputs/PROB.txt", "rb")

dna = rosalind.DNA("",f.next().strip())
gc = map(float,f.next().strip().split())

print " ".join([str(dna.log_probability(x)) for x in gc])