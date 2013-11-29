import rosalind

f = open("sample_inputs/LGIS.txt")
g = open("answer.txt", "wb")
f.next()
series = map(int,f.next().strip().split())
inc = rosalind.longest_monotonic_subsequence(series)
dec = rosalind.longest_monotonic_subsequence(series,decreasing=True)
g.write(" ".join(map(str,inc)))
g.write("\n")
g.write(" ".join(map(str,dec)))