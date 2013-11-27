"""
Mendel's First Law (IPRB)
"""
import rosalind

f = open('sample_inputs/IPRB.txt', 'r')
k,m,n = map(int,f.next().strip().split())

genotypes = rosalind.punnett_probabilities(k,m,n)

print genotypes[0]+genotypes[1]