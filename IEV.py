"""
Calculating Expected Offspring (IEV)
"""
f = open('sample_inputs/IEV.txt', 'r')
counts = f.readline().strip().split()
intCounts = [int(x) for x in counts]
expected = sum(intCounts[0:3])*2+intCounts[3]*1.5+intCounts[4]
print expected