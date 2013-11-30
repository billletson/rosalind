import rosalind

f = open("sample_inputs/PPER.txt", "rb")
n, k = map(int,f.next().strip().split())

print rosalind.partial_permutation(n,k,1000000)

