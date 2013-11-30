import rosalind


def prob_DNA(mode="sample"):
    """
    Counting DNA Nucleotides
    """
    f = open(mode + "_inputs/DNA.txt", "rb")
    g = open("answer.txt", "wb")
    dna = rosalind.DNA("", f.next().strip())
    g.write(" ".join([str(dna.count("A")),
        str(dna.count("C")),
        str(dna.count("G")),
        str(dna.count("T"))]))
    f.close()
    g.close()


def prob_RNA(mode="sample"):
    """
    Transcribing DNA into RNA
    """
    f = open(mode + "_inputs/RNA.txt","rb")
    g = open("answer.txt", "wb")
    dna = rosalind.DNA("", f.next().strip())
    rna = dna.to_RNA()
    g.write(rna.sequence)
    f.close()
    g.close()


def prob_REVC(mode="sample"):
    """
    Complementing a Strand of DNA
    """
    f = open(mode + "_inputs/REVC.txt","rb")
    g = open("answer.txt", "wb")
    dna = rosalind.DNA("", f.next().strip())
    comp = dna.reverse_complement()
    g.write(comp.sequence)
    f.close()
    g.close()


def prob_FIB(mode="sample"):
    """
    Rabbits and Recurrence Relations
    """
    f = open(mode + "_inputs/FIB.txt", "rb")
    g = open("answer.txt", "wb")
    n, k = map(int, f.next().strip().split())
    rabbits = [1, 1]
    for i in range(2, n):
        rabbits.append(rabbits[i-2]*k + rabbits[i-1])
    g.write(str(rabbits[-1]))
    f.close()
    g.close()


def prob_GC(mode="sample"):
    """
    Computing GC Content
    """
    sequences = rosalind.load_FASTA_file(open(mode + "_inputs/GC.txt", "rb"))
    g = open("answer.txt", "wb")
    gc_content = {}
    for s in sequences:
        gc_content[s.name] = (s.count("G")+s.count("C"))/float(len(s)) * 100
    m = max(gc_content,key=gc_content.get)
    g.write(m + "\r\n")
    g.write("%.6f" % gc_content[m])
    g.close()


def prob_HAMM(mode="sample"):
    """
    Counting Point Mutations
    """
    f = open(mode + "_inputs/HAMM.txt","rb")
    g = open("answer.txt", "wb")
    dna1 = rosalind.DNA("", f.next().strip())
    dna2 = rosalind.DNA("", f.next().strip())
    g.write(str(rosalind.hamming(dna1, dna2)))
    f.close()
    g.close()


def prob_IPRB(mode="sample"):
    """
    Mendel's First Law
    """
    f = open(mode + "_inputs/IPRB.txt", "rb")
    g = open("answer.txt", "wb")
    k, m, n = map(int, f.next().strip().split())
    genotypes = rosalind.punnett_probabilities(k, m, n)
    g.write("%.5f" % (genotypes[0] + genotypes[1]))
    f.close()
    g.close()


def prob_PROT(mode="sample"):
    """
    Translating RNA into Protein
    """
    f = open(mode + "_inputs/PROT.txt", "rb")
    g = open("answer.txt", "wb")
    rna = rosalind.RNA("", f.next().strip())
    proteins = rna.to_protein()
    g.write(proteins[0].sequence)
    f.close()
    g.close()


def prob_SUBS(mode="sample"):
    """
    Finding a Motif in DNA
    """
    f = open(mode + "_inputs/SUBS.txt","rb")
    g = open("answer.txt", "wb")
    dna = rosalind.DNA("", f.next().strip())
    g.write(" ".join([str(x) for x in dna.find_substring_locations(f.next().strip(), True)]))
    f.close()
    g.close()


def prob_CONS(mode="sample"):
    """
    Consensus and Profile
    """
    dnas = rosalind.load_FASTA_file(open(mode + "_inputs/CONS.txt", "rb"))
    g = open("answer.txt", "wb")
    matrix, string = rosalind.consensus(dnas)
    g.write(string)
    g.write("\r\n")
    g.write("A: " + " ".join([str(x['A']) for x in matrix]))
    g.write("\r\n")
    g.write("C: " + " ".join([str(x['C']) for x in matrix]))
    g.write("\r\n")
    g.write("G: " + " ".join([str(x['G']) for x in matrix]))
    g.write("\r\n")
    g.write("T: " + " ".join([str(x['T']) for x in matrix]))
    g.close()


def prob_FIBD(mode="sample"):
    """
    Mortal Fibonacci Rabbits (FIBD)
    """
    f = open(mode + "_inputs/FIBD.txt", "rb")
    g = open("answer.txt", "wb")
    n, m = f.readline().strip().split()
    n = int(n)
    m = int(m)
    rabbits = [[0, 0, 0] for x in xrange(n+m)]
    rabbits[0] = [1, 0, 0]
    rabbits[m-1] = [0, 0, 1]
    for x in range(1, n):
        rabbits[x][0] = rabbits[x-1][1]
        rabbits[x][1] = rabbits[x-1][1] + rabbits[x-1][0] - rabbits[x-1][2]
        rabbits[x+m-1][2] = rabbits[x][0]
    g.write(str(rabbits[n-1][0] + rabbits[n-1][1]))
    f.close()
    g.close()

def prob_GRPH(mode="sample"):
    """
    Overlap Graphs
    """
    f = open(mode + "_inputs/GRPH.txt", "rb")
    g = open("answer.txt", "wb")
    dnas = rosalind.load_FASTA_file(f)
    g.write("\r\n".join([" ".join(x) for x in rosalind.overlap_graph(dnas,3)]))
    f.close()
    g.close()


def prob_IEV(mode="sample"):
    """
    Calculating Expected Offspring
    """
    f = open(mode + "_inputs/IEV.txt", "rb")
    g = open("answer.txt", "wb")
    counts = map(int, f.readline().strip().split())
    expected = sum(counts[0:3])*2+counts[3]*1.5+counts[4]
    g.write(str(expected))
    f.close()
    g.close()


def prob_LCSM(mode="sample"):
    """
    Finding a Shared Motif
    """
    f = open(mode + "_inputs/LCSM.txt", "rb")
    g = open("answer.txt", "wb")
    dnas = rosalind.load_FASTA_file(f)
    g.write(rosalind.longest_common_substring(dnas))
    f.close()
    g.close()


def prob_LIA(mode="sample"):
    """
    Independent Alleles
    """
    f = open(mode + "_inputs/LIA.txt", "rb")
    g = open("answer.txt", "wb")
    from math import factorial
    k, s = map(int, f.readline().strip().split())
    n = 2**k
    p = 0.25
    prob = 0
    for i in xrange(n, s-1, -1):
        prob += (factorial(n) * p**i * (1-p)**(n-i)) / float(factorial(i)*factorial(n-i))
    g.write("%.3f" % prob)
    f.close()
    g.close()


def prob_MPRT(mode="sample"):
    """
    Finding a Protein Motif
    """
    f = open(mode + "_inputs/MPRT.txt", "rb")
    g = open("answer.txt", "wb")
    sub = 'N[^P][ST][^P]'
    write_lines = []
    for line in f:
        protein = rosalind.load_FASTA_uniprot(line.strip())
        matches = [str(x) for x in protein.find_substring_locations(sub,True)]
        if len(matches)>0:
            write_lines.append(protein.name)
            write_lines.append(" ".join(matches))
    g.write("\r\n".join(write_lines))
    f.close()
    g.close()


def prob_MRNA(mode="sample"):
    """
    Inferring mRNA from Protein
    """
    f = open(mode + "_inputs/MRNA.txt", "rb")
    g = open("answer.txt", "wb")
    protein = rosalind.Protein("", f.next().strip())
    g.write(str(protein.infer_RNA(1000000)))
    f.close()
    g.close()


def prob_ORF(mode="sample"):
    """
    Open Reading Frames
    """
    f = open(mode + "_inputs/ORF.txt", "rb")
    g = open("answer.txt", "wb")
    original = rosalind.load_FASTA_file(f)[0]
    reverse = original.reverse_complement()
    dna_frames = original.reading_frames() + reverse.reading_frames()
    rna_frames = [x.to_RNA() for x in dna_frames]
    proteins = []
    for rna in rna_frames:
        proteins += rna.to_protein()
    uniques = set([p.sequence for p in proteins])
    g.write("\r\n".join(uniques))
    f.close()
    g.close()


def prob_PERM(mode="sample"):
    """
    Enumerating Gene Orders
    """
    import itertools
    f = open(mode + "_inputs/PERM.txt", "rb")
    g = open("answer.txt", "wb")
    x = int(f.readline().strip())
    integers = [a for a in xrange(1,x+1)]
    permutations = list(itertools.permutations(integers))
    g.write(str(len(permutations))+"\r\n")
    g.write("\r\n".join([" ".join([str(a) for a in list(perm)]) for perm in permutations]))
    f.close()
    g.close()


def prob_PRTM(mode="sample"):
    """
    Calculating Protein Mass
    """
    f = open(mode + "_inputs/PRTM.txt", "rb")
    g = open("answer.txt", "wb")
    protein = rosalind.Protein("", f.next().strip())
    g.write("%.3f" % protein.mass())
    f.close()
    g.close()


def prob_REVP(mode="sample"):
    """
    Locating Restriction Sites
    """
    f = open(mode + "_inputs/REVP.txt", "rb")
    g = open("answer.txt", "wb")
    dna = rosalind.load_FASTA_file(f)[0]
    g.write("\r\n".join([" ".join([str(d[1]), str(d[2])]) for d in dna.reverse_palindromes(4, 12, True)]))
    f.close()
    g.close()


def prob_SPLC(mode="sample"):
    """
    RNA Splicing
    """
    f = open(mode + "_inputs/SPLC.txt", "rb")
    g = open("answer.txt", "wb")
    dnas = rosalind.load_FASTA_file(f)
    rnas = [dna.to_RNA() for dna in dnas]
    main_rna = rnas[0]
    for rna in rnas[1:]:
        main_rna = main_rna.remove_intron(rna.sequence)
    protein = main_rna.to_protein()[0]
    g.write(protein.sequence)
    f.close()
    g.close()


def prob_LEXF(mode="sample"):
    """
    Enumerating k-mers Lexicographically
    """
    f = open(mode + "_inputs/LEXF.txt", "rb")
    g = open("answer.txt", "wb")
    letters = f.next().strip().split()
    n = int(f.next().strip())
    g.write("\r\n".join(["".join(x) for x in rosalind.lexicographic_permutations(letters,n)]))
    f.close()
    g.close()


def prob_LGIS(mode="sample"):
    """
    Longest Increasing Subsequence
    """
    f = open(mode + "_inputs/LGIS.txt", "rb")
    g = open("answer.txt", "wb")
    f.next()
    series = map(int,f.next().strip().split())
    inc = rosalind.longest_monotonic_subsequence(series)
    dec = rosalind.longest_monotonic_subsequence(series,decreasing=True)
    g.write(" ".join(map(str,inc)))
    g.write("\r\n")
    g.write(" ".join(map(str,dec)))
    f.close()
    g.close()


def prob_LONG(mode="sample"):
    """
    Genome Assembly as Shortest Superstring
    """
    g = open("answer.txt", "wb")
    dnas = rosalind.load_FASTA_file(open(mode + "_inputs/LONG.txt","rb"))
    superstring, leftovers = rosalind.superstring(dnas)
    g.write(superstring.sequence)
    g.close()


def prob_PMCH(mode="sample"):
    """
    Perfect Matchings and RNA Secondary Structures
    """
    f = open(mode + "_inputs/PMCH.txt", "rb")
    g = open("answer.txt", "wb")
    rnas = rosalind.load_FASTA_file(f, "RNA")
    g.write(str(rnas[0].perfect_matchings()))
    f.close()
    g.close()


def prob_PPER(mode="sample"):
    """
    Partial Permutations
    """
    f = open(mode + "_inputs/PPER.txt", "rb")
    g = open("answer.txt", "wb")
    n, k = map(int, f.next().strip().split())
    g.write(str(rosalind.partial_permutation(n, k, 1000000)))
    f.close()
    g.close()


def prob_PROB(mode="sample"):
    """
    Introduction to Random Strings
    """
    f = open(mode + "_inputs/PROB.txt", "rb")
    g = open("answer.txt", "wb")
    dna = rosalind.DNA("", f.next().strip())
    gc = map(float, f.next().strip().split())
    g.write(" ".join(["%.3f" % dna.log_probability(x) for x in gc]))
    f.close()
    g.close()


def prob_SIGN(mode="sample"):
    """
    Enumerating Oriented Gene Orderings
    """
    f = open(mode + "_inputs/SIGN.txt", "rb")
    g = open("answer.txt", "wb")
    n = int(f.next().strip())
    perms = rosalind.signed_permutations(n)
    g.write(str(len(perms))+"\r\n")
    g.write("\r\n".join([" ".join(map(str,perm)) for perm in perms]))
    f.close()
    g.close()


def prob_SSEQ(mode="sample"):
    """
    Finding a Spliced Motif
    """
    f = open(mode + "_inputs/SSEQ.txt", "rb")
    g = open("answer.txt", "wb")
    dnas = rosalind.load_FASTA_file(f)
    g.write(" ".join(map(str, dnas[0].find_subsequence_locations(dnas[1].sequence, True))))
    f.close()
    g.close()


def prob_TRAN(mode="sample"):
    """
    Transitions and Transversions
    """
    f = open(mode + "_inputs/TRAN.txt", "rb")
    g = open("answer.txt", "wb")
    dnas = rosalind.load_FASTA_file(f)
    transitions, transversions = rosalind.transition_transversion(dnas[0], dnas[1])
    g.write(str(float(transitions)/transversions))
    f.close()
    g.close()


def prob_TREE(mode="sample"):
    """
    Completing a Tree
    """
    f = open(mode + "_inputs/TREE.txt", "rb")
    g = open("answer.txt", "wb")
    n = int(f.next().strip())
    edges = [map(int,x.strip().split()) for x in f]
    g.write(str(rosalind.n_connected_subgraphs(list(range(1,n+1)), edges)-1))
    f.close()
    g.close()
