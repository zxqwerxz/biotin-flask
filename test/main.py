import os, pysam, sys, collections

dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'files/sample3.bam')

# samfile = pysam.AlignmentFile(filename, "rb")

""""
for gene in samfile.references:
    print gene
    break

for length in samfile.lengths:
    print length
    break
"""


class BamFile(pysam.AlignmentFile):
    """Wrapper around the AlignmentFile class."""

    def __init__(self, filepath, *args, **kwargs):
        """Construct a BAM file."""
        super(BamFile, self).__init__(filename=filepath, *args, **kwargs)
        self.filepath = filepath


samfile = BamFile(filename, "rb")
for gene in samfile.references:
    print gene
    break

for length in samfile.lengths:
    print length
    break
