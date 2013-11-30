import rosalind

sequences = rosalind.load_FASTA_file(open("sample_inputs/GC.txt","rb"))

gc_content = {}

for s in sequences:
    gc_content[s.name] = (s.count("G")+s.count("C"))/float(len(s)) * 100

m = max(gc_content,key=gc_content.get)
print m
print gc_content[m]