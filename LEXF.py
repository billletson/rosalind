import rosalind

f = open("test_inputs/LEXF.txt","rb")

letters = f.next().strip().split()
n = int(f.next().strip())

for x in rosalind.lexographic_permutations(letters,n):
	print "".join(x)