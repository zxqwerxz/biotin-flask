import os, pysam
from unittest import TestCase
from biotin_flask.models.pysam import get_indel, ref2querypos

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

class Test_Ref2QueryPos(TestCase):
    def setUp(self):
        dir = os.path.dirname(__file__)
        self.filename = os.path.join(dir, 'files/sample.bam')
        self.bamfile = pysam.AlignmentFile(self.filename, 'rb')
        self.read = None
        self.foc = 3

    def test_get_indel_matches_pileupread_pos(self):
        first = True
        for col in self.bamfile.pileup(region="IL4:5879-5880"):
            count = 0
            for r in col.pileups:
                if first:
                    if count == self.foc:
                        self.read = r
                        first = False
                    count = count + 1
                if self.read is not None and self.read.alignment.query_name == r.alignment.query_name:
                    if not r.is_del:
                        self.assertEquals(r.query_position, ref2querypos(r.alignment, col.reference_pos))

    def test_fetch_decimal(self):
        first = True
        for col in self.bamfile.pileup(region="IL4:5879-5880"):
            count = 0
            for r in col.pileups:
                if first:
                    if count == self.foc:
                        self.read = r
                        first = False
                    count = count + 1
        self.assertEquals(82, ref2querypos(self.read.alignment, 5959.01))
