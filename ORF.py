import rosalind

f = open("sample_inputs/ORF.txt","rb")

original = rosalind.load_FASTA_file(f)[0]
reverse = original.reverse_complement()

dna_frames = original.reading_frames() + reverse.reading_frames()
rna_frames = [x.to_RNA() for x in dna_frames]

proteins = []
for rna in rna_frames:
	proteins += rna.to_protein()
	
uniques = set([p.sequence for p in proteins])
for u in uniques:
	print u

