import os, pysam

dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'files/sample.sam')

#pysam.sort(filename, "output")
#pysam.index("output.bam")

samfile = pysam.AlignmentFile("output.bam", "rb")

#for read in samfile.pileup(reference="IL4", start=5879, end=5882):
#    print str(read)

# 1: Upload SAM/BAM - assume sam
# 2: Determine if SAM or BAM and procede with creating temporary files
    # bam = TempBam.create(filename)
    # bam.close() -- cleanup operation
# Two main workflows:
    # 1: Nogap (clonal software)
    # 2: Gapped (qc/visualization)
    # 3: ClustalW (multiple alignment)