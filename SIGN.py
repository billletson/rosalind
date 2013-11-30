import rosalind

f = open("sample_inputs/SIGN.txt")
g = open("answer.txt", "wb")

n = int(f.next().strip())

perms = rosalind.signed_permutations(n)

g.write(str(len(perms))+"\n")
for perm in perms:
    g.write(" ".join(map(str,perm))+"\n")
