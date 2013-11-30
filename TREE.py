import rosalind

f = open("sample_inputs/TREE.txt")
g = open("answer.txt", "wb")

n = int(f.next().strip())
edges = [map(int,x.strip().split()) for x in f]

g.write(str(rosalind.n_connected_subgraphs(list(range(1,n+1)), edges)-1))