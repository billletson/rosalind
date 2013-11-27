"""
Enumerating Gene Orders (PERM)
"""
import itertools
f=open('sample_inputs/PERM.txt', 'rb')
x = int(f.readline().strip())
integers = [a for a in xrange(1,x+1)]
permutations = list(itertools.permutations(integers))
print len(permutations)
for perm in permutations:
    print " ".join([str(a) for a in list(perm)])