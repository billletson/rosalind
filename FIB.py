"""
Rabbits and Recurrence Relations (FIB)
"""
f = open('sample_inputs/FIB.txt', 'rb')
n,k = map(int,f.next().strip().split())
rabbits = [1,1]
for i in range(2,n):
    rabbits.append(rabbits[i-2]*k + rabbits[i-1])
print rabbits[-1]

