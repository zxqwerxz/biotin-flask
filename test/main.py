import os, pysam, sys
sys.path.insert(0, "/Users/jeffrey/src/flask/biotin-flask/")
from biotin_flask.models.pysam import get_indel, ref2querypos

dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'files/sample.sam')

pysam.sort(filename, "output")
pysam.index("output.bam")

samfile = pysam.AlignmentFile("output.bam", "rb")

# THIS IS GOOD CODE
read = None
foc = 3
seq = ""
first = True
for pileup_column in samfile.pileup(region="IL4:5879-5880"):
    count = 0
    for r in pileup_column.pileups:
        if first:
            if count == foc:
                read = r
                first = False
            count = count + 1
        if read is not None and read.alignment.query_name == r.alignment.query_name:
            if not r.is_del:
                print str(r.query_position) + " " + str(ref2querypos(r.alignment, pileup_column.reference_pos)) + " " + str(pileup_column.reference_pos)
print read.alignment.cigarstring
print ref2querypos(read.alignment, 5959.01)

"""
read_dict['RB40W:00428:02640'] = []


ref_spacer = []
for each column:
    max_indel = 0
    for each read:
        find max indel
    ref_spacer.append((ref_pos, max_indel))

# Now generate printing view (numeric header)
fasta_ids = []
for pos in ref_spacer:
    fasta_ids.append(pos[0])
    if pos[1] > 0:
        for i in xrange(pos[1]):
            fasta_ids.append("-")

# Now generate reads
use cigar? -- I like this idea

"""