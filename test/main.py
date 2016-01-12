import os, pysam

dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'files/sample.sam')

pysam.sort(filename, "output")
pysam.index("output.bam")

samfile = pysam.AlignmentFile("output.bam", "rb")

for col in samfile.pileup(region="IL4:123232123-1234455"):
    print col.n

"""
for col in samfile.pileup(reference="IL4", start=5879, end=5882):
    print ("\ncoverage at base %s = %s" % (col.pos, col.n))
    for pileupread in col.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
            print ('\tbase in read %s = %s' %
                (pileupread.alignment.query_name,
                 pileupread.alignment.query_sequence[pileupread.query_position]))
samfile.close()
"""

# 1: Upload SAM/BAM - assume sam
# 2: Determine if SAM or BAM and procede with creating temporary files
    # bam = TempBam.create(filename)
    # bam.close() -- cleanup operation
# Two main workflows:
    # 1: Nogap (clonal software)
    # 2: Gapped (qc/visualization)
    # 3: ClustalW (multiple alignment)