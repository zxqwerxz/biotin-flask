import os, pysam
from unittest import TestCase

class ParseSamfile(TestCase):
    def setUp(self):
        dir = os.path.dirname(__file__)
        self.filename = os.path.join(dir, 'files/sample.sam')
        self.samfile = pysam.AlignmentFile(self.filename)

    def test_fileopen(self):
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        self.assertEquals(len(lines), 1125)

    def test_countlines(self):
        count = self.samfile.count(until_eof=True)
        self.assertEquals(count, 1123)

        # In order to loop through the file more than once, use the head method, because the pointer on the file
        # already moved and the samfile iterator/list will be empty otherwise.
        i = 0
        for read in self.samfile.head(count): i = i + 1
        self.assertEquals(i, 1123)

    def test_keysamtools(self):
        # These methods must work on the operating system
        #pysam.sort(filename, "output")
        #pysam.index("output.bam")
        #samfile = pysam.AlignmentFile("output.bam", "rb")
        pass

    def tearDown(self):
        self.samfile.close()
