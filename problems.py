import rosalind


def prob_dna(mode="sample"):
    """
    Counting DNA Nucleotides
    """
    with open(mode + "_inputs/DNA.txt", "rU") as f:
        dna = rosalind.DNA("", f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(" ".join([str(dna.count("A")),
                          str(dna.count("C")),
                          str(dna.count("G")),
                          str(dna.count("T"))]))


def prob_rna(mode="sample"):
    """
    Transcribing DNA into RNA
    """
    with open(mode + "_inputs/RNA.txt", "rU") as f:
        dna = rosalind.DNA("", f.next().strip())
    rna = dna.to_rna()
    with open("answer.txt", "wb") as g:
        g.write(rna.sequence)


def prob_revc(mode="sample"):
    """
    Complementing a Strand of DNA
    """
    with open(mode + "_inputs/REVC.txt", "rU") as f:
        dna = rosalind.DNA("", f.next().strip())
    comp = dna.reverse_complement()
    with open("answer.txt", "wb") as g:
        g.write(comp.sequence)


def prob_fib(mode="sample"):
    """
    Rabbits and Recurrence Relations
    """
    with open(mode + "_inputs/FIB.txt", "rU") as f:
        n, k = map(int, f.next().strip().split())
    rabbits = [1, 1]
    for i in range(2, n):
        rabbits.append(rabbits[i - 2] * k + rabbits[i - 1])
    with open("answer.txt", "wb") as g:
        g.write(str(rabbits[-1]))


def prob_gc(mode="sample"):
    """
    Computing GC Content
    """
    with open(mode + "_inputs/GC.txt", "rU") as f:
        sequences = rosalind.load_fasta_file(f)
    gc_content = {}
    for s in sequences:
        gc_content[s.name] = (s.count("G") + s.count("C")) / float(len(s)) * 100
    m = max(gc_content, key=gc_content.get)
    with open("answer.txt", "wb") as g:
        g.write(m + "\r\n")
        g.write("%.6f" % gc_content[m])


def prob_hamm(mode="sample"):
    """
    Counting Point Mutations
    """
    with open(mode + "_inputs/HAMM.txt", "rU") as f:
        dna1 = rosalind.DNA("", f.next().strip())
        dna2 = rosalind.DNA("", f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.hamming(dna1, dna2)))


def prob_iprb(mode="sample"):
    """
    Mendel's First Law
    """
    with open(mode + "_inputs/IPRB.txt", "rU") as f:
        k, m, n = map(int, f.next().strip().split())
    genotypes = rosalind.punnett_probabilities(k, m, n)
    with open("answer.txt", "wb") as g:
        g.write("%.5f" % (genotypes[0] + genotypes[1]))


def prob_prot(mode="sample"):
    """
    Translating RNA into Protein
    """
    with open(mode + "_inputs/PROT.txt", "rU") as f:
        rna = rosalind.RNA("", f.next().strip())
    proteins = rna.to_proteins()
    with open("answer.txt", "wb") as g:
        g.write(proteins[0].sequence)


def prob_subs(mode="sample"):
    """
    Finding a Motif in DNA
    """
    with open(mode + "_inputs/SUBS.txt", "rU") as f:
        dna = rosalind.DNA("", f.next().strip())
        with open("answer.txt", "wb") as g:
            g.write(" ".join([str(x) for x in dna.find_substring_locations(f.next().strip(), True)]))


def prob_cons(mode="sample"):
    """
    Consensus and Profile
    """
    with open(mode + "_inputs/CONS.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    matrix, string = rosalind.consensus(dnas)
    with open("answer.txt", "wb") as g:
        g.write(string)
        g.write("\r\n")
        g.write("A: " + " ".join([str(x['A']) for x in matrix]))
        g.write("\r\n")
        g.write("C: " + " ".join([str(x['C']) for x in matrix]))
        g.write("\r\n")
        g.write("G: " + " ".join([str(x['G']) for x in matrix]))
        g.write("\r\n")
        g.write("T: " + " ".join([str(x['T']) for x in matrix]))


def prob_fibd(mode="sample"):
    """
    Mortal Fibonacci Rabbits (FIBD)
    """
    from collections import defaultdict
    with open(mode + "_inputs/FIBD.txt", "rU") as f:
        n, m = map(int, f.readline().strip().split())
    rabbits = defaultdict(lambda: [0, 0, 0])
    rabbits[0] = [1, 0, 0]
    rabbits[m - 1] = [0, 0, 1]
    for x in range(1, n):
        rabbits[x][0] = rabbits[x - 1][1]
        rabbits[x][1] = rabbits[x - 1][1] + rabbits[x - 1][0] - rabbits[x - 1][2]
        rabbits[x + m - 1][2] = rabbits[x][0]
    with open("answer.txt", "wb") as g:
        g.write(str(rabbits[n - 1][0] + rabbits[n - 1][1]))


def prob_grph(mode="sample"):
    """
    Overlap Graphs
    """
    with open(mode + "_inputs/GRPH.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join([" ".join(x) for x in rosalind.overlap_graph(dnas, 3)]))


def prob_iev(mode="sample"):
    """
    Calculating Expected Offspring
    """
    with open(mode + "_inputs/IEV.txt", "rU") as f:
        counts = map(int, f.readline().strip().split())
    expected = sum(counts[0:3]) * 2 + counts[3] * 1.5 + counts[4]
    with open("answer.txt", "wb") as g:
        g.write(str(expected))


def prob_lcsm(mode="sample"):
    """
    Finding a Shared Motif
    """
    with open(mode + "_inputs/LCSM.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write(rosalind.longest_common_substring(dnas))


def prob_lia(mode="sample"):
    """
    Independent Alleles
    """
    from math import factorial
    with open(mode + "_inputs/LIA.txt", "rU") as f:
        k, s = map(int, f.readline().strip().split())
    n = 2 ** k
    p = 0.25
    prob = 0
    for i in xrange(n, s - 1, -1):
        prob += (factorial(n) * p ** i * (1 - p) ** (n - i)) / float(factorial(i) * factorial(n - i))
    with open("answer.txt", "wb") as g:
        g.write("%.3f" % prob)


def prob_mprt(mode="sample"):
    """
    Finding a Protein Motif
    """
    sub = 'N[^P][ST][^P]'
    write_lines = []
    with open(mode + "_inputs/MPRT.txt", "rU") as f:
        for line in f:
            protein = rosalind.load_fasta_uniprot(line.strip())
            matches = [str(x) for x in protein.find_substring_locations(sub, True)]
            if len(matches) > 0:
                write_lines.append(protein.name)
                write_lines.append(" ".join(matches))
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(write_lines))


def prob_mrna(mode="sample"):
    """
    Inferring mRNA from Protein
    """
    with open(mode + "_inputs/MRNA.txt", "rU") as f:
        protein = rosalind.Protein("", f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(str(protein.infer_rna(1000000)))


def prob_orf(mode="sample"):
    """
    Open Reading Frames
    """
    with open(mode + "_inputs/ORF.txt", "rU") as f:
        original = rosalind.load_fasta_file(f)[0]
    reverse = original.reverse_complement()
    dna_frames = original.reading_frames() + reverse.reading_frames()
    rna_frames = [x.to_rna() for x in dna_frames]
    proteins = []
    for rna in rna_frames:
        proteins += rna.to_proteins()
    uniques = set([p.sequence for p in proteins])
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(uniques))


def prob_perm(mode="sample"):
    """
    Enumerating Gene Orders
    """
    import itertools
    with open(mode + "_inputs/PERM.txt", "rU") as f:
        x = int(f.readline().strip())
    permutations = list(itertools.permutations(xrange(1, x + 1)))
    with open("answer.txt", "wb") as g:
        g.write(str(len(permutations)) + "\r\n")
        g.write("\r\n".join([" ".join([str(a) for a in list(perm)]) for perm in permutations]))


def prob_prtm(mode="sample"):
    """
    Calculating Protein Mass
    """
    with open(mode + "_inputs/PRTM.txt", "rU") as f:
        protein = rosalind.Protein("", f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write("%.3f" % protein.mass)


def prob_revp(mode="sample"):
    """
    Locating Restriction Sites
    """
    with open(mode + "_inputs/REVP.txt", "rU") as f:
        dna = rosalind.load_fasta_file(f)[0]
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join([" ".join([str(d[1]), str(d[2])]) for d in dna.reverse_palindromes(4, 12, True)]))


def prob_splc(mode="sample"):
    """
    RNA Splicing
    """
    with open(mode + "_inputs/SPLC.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    rnas = [dna.to_rna() for dna in dnas]
    main_rna = rnas[0]
    for rna in rnas[1:]:
        main_rna = main_rna.remove_intron(rna.sequence)
    protein = main_rna.to_proteins()[0]
    with open("answer.txt", "wb") as g:
        g.write(protein.sequence)


def prob_lexf(mode="sample"):
    """
    Enumerating k-mers Lexicographically
    """
    with open(mode + "_inputs/LEXF.txt", "rU") as f:
        letters = f.next().strip().split()
        n = int(f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(rosalind.lexicographic_permutations(letters, n)))


def prob_lgis(mode="sample"):
    """
    Longest Increasing Subsequence
    """
    with open(mode + "_inputs/LGIS.txt", "rU") as f:
        f.next()
        series = map(int, f.next().strip().split())
    inc = rosalind.longest_monotonic_subsequence(series)
    dec = rosalind.longest_monotonic_subsequence(series, decreasing=True)
    with open("answer.txt", "wb") as g:
        g.write(" ".join(map(str, inc)))
        g.write("\r\n")
        g.write(" ".join(map(str, dec)))


def prob_long(mode="sample"):
    """
    Genome Assembly as Shortest Superstring
    """
    with open(mode + "_inputs/LONG.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    superstring, leftovers = rosalind.superstring(dnas)
    with open("answer.txt", "wb") as g:
        g.write(superstring.sequence)


def prob_pmch(mode="sample"):
    """
    Perfect Matchings and RNA Secondary Structures
    """
    with open(mode + "_inputs/PMCH.txt", "rU") as f:
        rnas = rosalind.load_fasta_file(f, "RNA")
    with open("answer.txt", "wb") as g:
        g.write(str(rnas[0].perfect_matchings()))


def prob_pper(mode="sample"):
    """
    Partial Permutations
    """
    with open(mode + "_inputs/PPER.txt", "rU") as f:
        n, k = map(int, f.next().strip().split())
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.partial_permutation(n, k, 1000000)))


def prob_prob(mode="sample"):
    """
    Introduction to Random Strings
    """
    with open(mode + "_inputs/PROB.txt", "rU") as f:
        dna = rosalind.DNA("", f.next().strip())
        gc = map(float, f.next().strip().split())
    with open("answer.txt", "wb") as g:
        g.write(" ".join(["%.3f" % dna.log_probability(x) for x in gc]))


def prob_sign(mode="sample"):
    """
    Enumerating Oriented Gene Orderings
    """
    with open(mode + "_inputs/SIGN.txt", "rU") as f:
        n = int(f.next().strip())
    perms = rosalind.signed_permutations(n)
    with open("answer.txt", "wb") as g:
        g.write(str(len(perms)) + "\r\n")
        g.write("\r\n".join([" ".join(map(str, perm)) for perm in perms]))


def prob_sseq(mode="sample"):
    """
    Finding a Spliced Motif
    """
    with open(mode + "_inputs/SSEQ.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write(" ".join(map(str, dnas[0].find_subsequence_locations(dnas[1].sequence, True))))


def prob_tran(mode="sample"):
    """
    Transitions and Transversions
    """
    with open(mode + "_inputs/TRAN.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    transitions, transversions = rosalind.transition_transversion(dnas[0], dnas[1])
    with open("answer.txt", "wb") as g:
        g.write(str(float(transitions) / transversions))


def prob_tree(mode="sample"):
    """
    Completing a Tree
    """
    with open(mode + "_inputs/TREE.txt", "rU") as f:
        n = int(f.next().strip())
        edges = [map(int, x.strip().split()) for x in f]
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.n_connected_subgraphs(list(range(1, n + 1)), edges) - 1))


def prob_inod(mode="sample"):
    """
    Counting Phylogenetic Ancestors
    """
    with open(mode + "_inputs/INOD.txt", "rU") as f:
        leaves = int(f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.unrooted_internal_from_leaves(leaves)))


def prob_mmch(mode="sample"):
    """
    Maximum Matchings and RNA Secondary Structures
    """
    with open(mode + "_inputs/MMCH.txt", "rU") as f:
        rna = rosalind.load_fasta_file(f, "RNA")[0]
    with open("answer.txt", "wb") as g:
        g.write(str(rna.maximum_matchings()))


def prob_sset(mode="sample"):
    """
    Counting Subsets
    """
    with open(mode + "_inputs/SSET.txt", "rU") as f:
        n = int(f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.subset_count(n) % 1000000))


def prob_cat(mode="sample"):
    """
    Catalan Numbers and RNA Secondary Structures
    """
    with open(mode + "_inputs/CAT.txt", "rU") as f:
        rna = rosalind.load_fasta_file(f, "RNA")[0]
    with open("answer.txt", "wb") as g:
        g.write(str(rna.noncrossing_matchings() % 1000000))


def prob_kmp(mode="sample"):
    """
    Speeding Up Motif Finding
    """
    with open(mode + "_inputs/KMP.txt", "rU") as f:
        dna = rosalind.load_fasta_file(f)[0]
    with open("answer.txt", "wb") as g:
        g.write(" ".join(map(str, dna.failure_array())))


def prob_pdst(mode="sample"):
    """
    Creating a Distance Matrix
    """
    with open(mode + "_inputs/PDST.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join([" ".join(map("{:.5f}".format, x)) for x in rosalind.difference_matrix(dnas)]))


def prob_rstr(mode="sample"):
    """
    Matching Random Motifs
    """
    with open(mode + "_inputs/RSTR.txt", "rU") as f:
        attempts, gc = map(float, f.next().strip().split())
        dna = rosalind.DNA("", f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write("%.3f" % dna.probability_with_repeated_attempts(attempts, gc))


def prob_corr(mode="sample"):
    """
    Error Correction in Reads
    """
    with open(mode + "_inputs/CORR.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(["%s->%s" % (x[0].sequence, x[1].sequence) for x in rosalind.identify_read_errors(dnas)]))


def prob_kmer(mode="sample"):
    """
    k-Mer Composition
    """
    with open(mode + "_inputs/KMER.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write(" ".join([str(x) for x in dnas[0].k_mer_composition(4)]))


def prob_lexv(mode="sample"):
    """
    Ordering Strings of Varying Length Lexicographically
    """
    with open(mode + "_inputs/LEXV.txt", "rU") as f:
        alphabet = f.next().strip().split()
        length = int(f.next().strip().split()[0])
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(rosalind.multilength_lexicographic_permutations(alphabet, length)))


def prob_rear(mode="sample"):
    """
    Reversal Distance
    """
    with open(mode + "_inputs/REAR.txt", "rU") as f:
        pairs = []
        tmp = []
        for line in f:
            if line.strip() == "":
                pairs.append(tmp)
                tmp = []
            else:
                tmp.append(line.strip().split())
        pairs.append(tmp)
    with open("answer.txt", "wb") as g:
        g.write(" ".join([str(len(rosalind.find_reversals(x[0], x[1]))) for x in pairs]))


def prob_lcsq(mode="sample"):
    """
    Finding a Shared Spliced Motif
    """
    with open(mode + "_inputs/LCSQ.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write(rosalind.longest_common_subsequence(dnas[0], dnas[1]).sequence)


def prob_sort(mode="sample"):
    """
    Sorting by Reversals
    """
    with open(mode + "_inputs/SORT.txt", "rU") as f:
        original = map(int, f.next().strip().split())
        target = map(int, f.next().strip().split())
    with open("answer.txt", "wb") as g:
        answer = rosalind.find_reversals(original, target)
        answer = [(x[0] + 1, x[1] + 1) for x in answer]
        g.write("\r\n".join([str(len(answer))] + [" ".join(map(str, x)) for x in answer]))


def prob_edit(mode="sample"):
    """
    Edit Distance
    """
    with open(mode + "_inputs/EDIT.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.levenshtein(dnas[0], dnas[1])))


def prob_eval(mode="sample"):
    """
    Expected Number of Restriction Sites
    """
    with open(mode + "_inputs/EVAL.txt", "rU") as f:
        str_len = int(f.next().strip())
        substring = f.next().strip()
        gc_contents = map(float, f.next().strip().split())
    with open("answer.txt", "wb") as g:
        g.write(" ".join(["%0.3f" % rosalind.expected_restriction_sites(str_len, substring, gc) for gc in gc_contents]))


def prob_aspc(mode="sample"):
    """
    Introduction to Alternative Splicing
    """
    with open(mode + "_inputs/ASPC.txt", "rU") as f:
        n, m = map(int, f.next().strip().split())
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.subset_count(n, m) % 1000000))


def prob_seto(mode="sample"):
    """
    Introduction to Set Operations
    """
    with open(mode + "_inputs/SETO.txt", "rU") as f:
        base = set(map(str, xrange(1, int(f.next().strip()) + 1)))
        a = set(f.next().strip().strip("{").strip("}").split(", "))
        b = set(f.next().strip().strip("{").strip("}").split(", "))
    with open("answer.txt", "wb") as g:
        g.write("{" + ", ".join(a.union(b)) + "}")
        g.write("\r\n")
        g.write("{" + ", ".join(a.intersection(b)) + "}")
        g.write("\r\n")
        g.write("{" + ", ".join(a.difference(b)) + "}")
        g.write("\r\n")
        g.write("{" + ", ".join(b.difference(a)) + "}")
        g.write("\r\n")
        g.write("{" + ", ".join(base.difference(a)) + "}")
        g.write("\r\n")
        g.write("{" + ", ".join(base.difference(b)) + "}")


def prob_spec(mode="sample"):
    """
    Inferring Protein from Spectrum
    """
    with open(mode + "_inputs/SPEC.txt", "rU") as f:
        weights = []
        for line in f:
            weights.append(float(line.strip()))
    with open("answer.txt", "wb") as g:
        g.write(rosalind.infer_protein_from_prefix_spectrum(weights).sequence)


def prob_scsp(mode="sample"):
    """
    Interleaving Two Motifs
    """
    with open(mode + "_inputs/SCSP.txt", "rU") as f:
        first = rosalind.DNA("", f.next().strip())
        second = rosalind.DNA("", f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(rosalind.supersequence(first, second).sequence)


def prob_motz(mode="sample"):
    """
    Catalan Numbers and RNA Secondary Structures
    """
    with open(mode + "_inputs/MOTZ.txt", "rU") as f:
        rna = rosalind.load_fasta_file(f, "RNA")[0]
    with open("answer.txt", "wb") as g:
        g.write(str(rna.noncrossing_matchings(False) % 1000000))


def prob_trie(mode="sample"):
    """
    Introduction to Pattern Matching
    """
    with open(mode + "_inputs/TRIE.txt", "rU") as f:
        words = [line.strip() for line in f]
    trie = rosalind.Trie(words)
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(sorted(trie.edges_as_strings(True))))


def prob_nwck(mode="sample"):
    """
    Distances in Trees
    """
    with open(mode + "_inputs/NWCK.txt", "rU") as f:
        sets = []
        newick = ""
        for line in f:
            line = line.strip()
            if line and line[-1] == ";":
                newick = line
            elif newick:
                first, second = line.split()
                sets.append((newick, first, second))
                newick = ""
    with open("answer.txt", "wb") as g:
        g.write(" ".join([str(rosalind.NewickTree(x[0]).find_distance(x[1], x[2])) for x in sets]))


def prob_nkew(mode="sample"):
    """
    Newick Format with Edge Weights
    """
    with open(mode + "_inputs/NKEW.txt", "rU") as f:
        sets = []
        newick = ""
        for line in f:
            line = line.strip()
            if line and line[-1] == ";":
                newick = line
            elif newick:
                first, second = line.split()
                sets.append((newick, first, second))
                newick = ""
    with open("answer.txt", "wb") as g:
        g.write(" ".join([str(int(rosalind.NewickTree(x[0]).find_distance(x[1], x[2]))) for x in sets]))


def prob_ctbl(mode="sample"):
    """
    Creating a Character Table
    """
    with open(mode + "_inputs/CTBL.txt", "rU") as f:
        tree = rosalind.NewickTree(f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(tree.nontrivial_characters()))


def prob_indc(mode="sample"):
    """
    Independent Segregation of Chromosomes
    """
    with open(mode + "_inputs/INDC.txt", "rU") as f:
        n = int(f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(" ".join(["%0.3f" % x for x in reversed(rosalind.binomial_cdf(2*n, 0.5, True)[0:-1])]))


def prob_dbru(mode="sample"):
    """
    Constructing a De Bruijn Graph
    """
    with open(mode + "_inputs/DBRU.txt", "rU") as f:
        dnas = [rosalind.DNA("", line.strip()) for line in f]
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(sorted(["(%s, %s)" % x for x in rosalind.debruijn(dnas)])))


def prob_conv(mode="sample"):
    """
    Comparing Spectra with the Spectral Convolution
    """
    with open(mode + "_inputs/CONV.txt", "rU") as f:
        spectra = [map(float, line.strip().split()) for line in f]
    with open("answer.txt", "wb") as g:
        mink = rosalind.minkowski_difference(spectra[0], spectra[1])
        g.write("\r\n".join([str(mink.count(mink[0])), str(mink[0])]))


def prob_rnas(mode="sample"):
    """
    Wobble Bonding and RNA Secondary Structures
    """
    with open(mode + "_inputs/RNAS.txt", "rU") as f:
        rna = rosalind.RNA("", f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(str(rna.noncrossing_matchings(False, True, 4)))


def prob_full(mode="sample"):
    """
    Inferring Peptide from Full Spectrum
    """
    with open(mode + "_inputs/FULL.txt", "rU") as f:
        total = float(f.next().strip())
        ions = [float(x.strip()) for x in f]
    with open("answer.txt", "wb") as g:
        g.write("".join(rosalind.infer_peptide_from_full_spectrum(total, ions).sequence))


def prob_lrep(mode="sample"):
    """
    Finding the Longest Multiple Repeat
    """
    with open(mode + "_inputs/LREP.txt", "rU") as f:
        tree = rosalind.SuffixTree(f.next().strip())
        k = int(f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write(tree.longest_occuring_k_times(k))


def prob_suff(mode="sample"):
    """
    Encoding Suffix Trees
    """
    with open(mode + "_inputs/SUFF.txt", "rU") as f:
        tree = rosalind.SuffixTree(f.next().strip())
    with open("answer.txt", "wb") as g:
        g.write("\r\n".join(tree.list_of_edges()))


def prob_itwv(mode="sample"):
    """
    Finding Disjoint Motifs in a Gene
    """
    with open(mode + "_inputs/ITWV.txt", "rU") as f:
        haystack = f.next().strip()
        needles = [line.strip() for line in f]
    with open("answer.txt", "wb") as g:
        matrix = rosalind.interweaveable_matrix(haystack, needles)
        g.write("\r\n".join([" ".join(map(str, x)) for x in matrix]))


def prob_edta(mode="sample"):
    """
    Edit Distance Alignment
    """
    with open(mode + "_inputs/EDTA.txt", "rU") as f:
        proteins = rosalind.load_fasta_file(f, s_type="Protein")
    with open("answer.txt", "wb") as g:
        score, s, t = rosalind.global_alignment(proteins[0], proteins[1])
        g.write("\r\n".join([str(-score), s, t]))


def prob_ctea(mode="sample"):
    """
    Counting Optimal Alignments
    """
    with open(mode + "_inputs/CTEA.txt", "rU") as f:
        proteins = rosalind.load_fasta_file(f, s_type="Protein")
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.optimal_alignment_count(proteins[0], proteins[1], 134217727)))


def prob_glob(mode="sample"):
    """
    Global Alignment with Scoring Matrix
    """
    with open(mode + "_inputs/GLOB.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write(str(rosalind.alignment_score(dnas[0], dnas[1], "BLOSUM62", 5)))


def prob_mgap(mode="sample"):
    """
    Maximizing the Gap Symbols of an Optimal Alignment
    """
    with open(mode + "_inputs/MGAP.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        g.write(str(len(dnas[0]) + len(dnas[1]) - len(rosalind.longest_common_subsequence(dnas[0], dnas[1])) * 2))


def prob_loca(mode="sample"):
    """
    Local Alignment with Scoring Matrix
    """
    with open(mode + "_inputs/LOCA.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        score, s, t = rosalind.local_alignment(dnas[0], dnas[1], "PAM250", 5)
        g.write("\r\n".join([str(score), s, t]))


def prob_gcon(mode="sample"):
    """
    Global Alignment with Constant Gap Penalty
    """
    with open(mode + "_inputs/GCON.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        score = rosalind.alignment_score(dnas[0], dnas[1], "BLOSUM62", (5, 0))
        g.write(str(score))


def prob_gaff(mode="sample"):
    """
    Global Alignment with Scoring Matrix and Affine Gap Penalty
    """
    with open(mode + "_inputs/GAFF.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        score, s, t = rosalind.global_alignment(dnas[0], dnas[1], "BLOSUM62", (11, 1))
        g.write("\r\n".join([str(score), s, t]))


def prob_laff(mode="sample"):
    """
    Local Alignment with Affine Gap Penalty
    """
    with open(mode + "_inputs/LAFF.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        score, s, t = rosalind.local_alignment(dnas[0], dnas[1], "BLOSUM62", (11, 1))
        g.write("\r\n".join([str(score), s, t]))


def prob_smgb(mode="sample"):
    """
    Semiglobal Alignment
    """
    with open(mode + "_inputs/SMGB.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        score, s, t = rosalind.semi_global_alignment(dnas[0], dnas[1], (1, -1), 1)
        g.write("\r\n".join([str(score), s, t]))


def prob_oap(mode="sample"):
    """
    Overlapping Reads with Errors
    """
    with open(mode + "_inputs/OAP.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        score, s, t = rosalind.overlap_alignment(dnas[0], dnas[1], (1, -2), 2)
        g.write("\r\n".join([str(score), s, t]))


def prob_sims(mode="sample"):
    """
    Finding a Motif with Modifications
    """
    with open(mode + "_inputs/SIMS.txt", "rU") as f:
        dnas = rosalind.load_fasta_file(f)
    with open("answer.txt", "wb") as g:
        score, s, t = rosalind.semi_global_alignment(dnas[0], dnas[1], (1, -1), 1, True)
        g.write("\r\n".join([str(score), s, t]))


def prob_ksim(mode="sample"):
    """
    Finding All Similar Motifs
    """
    with open(mode + "_inputs/KSIM.txt", "rU") as f:
        k = -int(f.next().strip())
        first = rosalind.DNA("", f.next().strip())
        second = rosalind.DNA("", f.next().strip())
    with open("answer.txt", "wb") as g:
        pairs = rosalind.all_semi_global_alignments(first, second, (0, -1), 1, k)
        g.write("\r\n".join([" ".join([str(pair[0]+1), str(pair[1])]) for pair in pairs]))