"""
Independent Alleles (LIA)
"""
from math import factorial
f = open('sample_inputs/LIA.txt', 'r')
k,s = map(int,f.readline().strip().split())
n = 2**k
p = 0.25
prob = 0
for i in xrange(n, s-1, -1):
    prob += (factorial(n) * p**i * (1-p)**(n-i)) / float(factorial(i)*factorial(n-i))
print prob