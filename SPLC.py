import rosalind

f = open("sample_inputs/SPLC.txt","rb")

dnas = rosalind.load_FASTA_file(f)
rnas = [dna.to_RNA() for dna in dnas]

main_rna = rnas[0]
for rna in rnas[1:]:
	main_rna = main_rna.remove_intron(rna.sequence)
	
protein = main_rna.to_protein()[0]

print protein.sequence