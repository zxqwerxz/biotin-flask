import os, pysam
from unittest import TestCase
from biotin_flask.models.pysam import get_indel

class Test_Get_Indel(TestCase):
    def setUp(self):
        dir = os.path.dirname(__file__)
        self.filename = os.path.join(dir, 'files/sample.bam')
        self.bamfile = pysam.AlignmentFile(self.filename, 'rb')

    def test_get_indel_matches_pileupread_indel(self):
        read = None
        foc = 3
        first = True
        for col in self.bamfile.pileup(region="IL4:5879-5880"):
            count = 0
            for r in col.pileups:
                if first:
                    if count == foc:
                        read = r
                        first = False
                    count = count + 1
                if read is not None and read.alignment.query_name == r.alignment.query_name:
                    if not r.is_del:
                        self.assertEquals(r.indel, get_indel(r.alignment, r.query_position))