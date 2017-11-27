import os, pysam, sys, collections

dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'files/sample2.bam')

samfile = pysam.AlignmentFile(filename, "rb")

for gene in samfile.references:
    print gene
    break

for length in samfile.lengths:
    print length
    break

# temp = samfile.fetch(gene, 0, 21516)

mylist = [''] * 0 + ['a'] + [''] * 0
print mylist

dicta = collections.OrderedDict({
    'a': 1,
    'b': 3
})

s = "Hello"

print s[0]
