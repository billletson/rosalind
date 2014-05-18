import rosalind
import random
import timeit

alphabet = "abcdefghij"
strings = []

for i in [10, 10, 10, 10]:
    strings.append("".join([alphabet[random.randint(0, len(alphabet) - 1)] for _ in xrange(i)]))

times = [timeit.timeit("rosalind.SuffixTree('%s')" % x, "import rosalind", number=1) for x in strings]

print times
#for i in xrange(1, len(times)):
#    print times[i]/times[0]
