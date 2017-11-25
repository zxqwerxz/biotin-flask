import os, pysam, sys

dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'files/sample.sam')

pysam.sort("-o", "output.bam", filename)
pysam.index("output.bam")

samfile = pysam.AlignmentFile("output.bam", "rb")
