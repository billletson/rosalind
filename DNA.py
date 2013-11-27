import rosalind

"""
Counting DNA Nucleotides (DNA)
"""
f = open("sample_inputs/DNA.txt","rb")
dna = rosalind.DNA("",f.next().strip())
print " ".join([str(dna.count("A")),
                str(dna.count("C")),
				str(dna.count("G")),
				str(dna.count("T"))])