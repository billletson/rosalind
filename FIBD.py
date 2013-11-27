"""
Mortal Fibonacci Rabbits (FIBD)
"""
f = open('sample_inputs/FIBD.txt', 'rb')
n,m = f.readline().strip().split(" ")
n = int(n)
m = int(m)
rabbits = [[0,0,0] for x in xrange(n+m)]
rabbits[0] = [1,0,0]
rabbits[m-1] = [0,0,1]
for x in range(1,n):
    rabbits[x][0] = rabbits[x-1][1]
    rabbits[x][1] = rabbits[x-1][1] + rabbits[x-1][0] - rabbits[x-1][2]
    rabbits[x+m-1][2] = rabbits[x][0]
print (rabbits[n-1][0] + rabbits[n-1][1])