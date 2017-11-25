"""Miscellaneous common shared utilities."""

import io
import os
import random
import string
import tempfile
import time
import threading
import zipfile

import pysam
from HTMLParser import HTMLParser
from flask import flash

from biotin_flask import app

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou', 'Eric Zhou']
__version__ = '0.0.5'
__status__ = 'Production'


class FastaUpload:
    """Class for handling any FASTA-like file upload."""

    @staticmethod
    def extract_seq(file_h, seqname):
        """Method for extracting the genomic sequence of a certain seqname/gene in the fasta file"""
        found_it = False
        for line in file_h.readlines():
            line = line.rstrip()
            if line == '>' + seqname:
                found_it = True
                continue
            if found_it:
                return line
        return None


class SamUpload:
    """Class for handling any SAM-like file upload. Automatically converts to an indexed BAMfile."""

    def __init__(self, file_h, filename):
        self.filelist = [] # Need to clean this up afterwards
        self.bamfiles = [] # Group of AlignmentFiles
        filefront, extension = os.path.splitext(filename)

        if extension == '.zip':
            # Save the zipfile
            file_h.save(os.path.join(tempfile.gettempdir(), filename))
            self.filelist.append(os.path.join(tempfile.gettempdir(), filename))
            print 'saved file: ' + filename
            print 'now extracting: ' + filename

            # Extract the zipfile
            with zipfile.ZipFile(os.path.join(tempfile.gettempdir(), filename), "r") as z:
                z.extractall(tempfile.gettempdir())
                for file in z.namelist():
                    self.filelist.append(os.path.join(tempfile.gettempdir(), file))
                    filefront, extension = os.path.splitext(file)
                    if extension == '.sam' and filefront[0] != "_":
                        self.__loadsamfile(file, filefront)

        elif extension == '.sam':
            # Save the zipfile
            file_h.save(os.path.join(tempfile.gettempdir(), filename))
            self.filelist.append(os.path.join(tempfile.gettempdir(), filename))
            print 'saved file: ' + filename
            file_path = os.path.join(tempfile.gettempdir(), filename)
            filefront, extension = os.path.splitext(file_path)
            self.__loadsamfile(file_path, filefront)

    def __loadsamfile(self, filename, filefront):
        # Sort the extracted SAM file into a compressed BAM file
        pysam.sort(
            os.path.join(tempfile.gettempdir(), filename),
            os.path.join(tempfile.gettempdir(), filefront)
        )
        self.filelist.append(os.path.join(tempfile.gettempdir(), filefront + '.bam'))

        # Index the generated BAM file
        pysam.index(os.path.join(tempfile.gettempdir(), filefront + '.bam'))
        self.filelist.append(os.path.join(tempfile.gettempdir(), filefront + '.bam.bai'))

        # Load the indexed BAM file as an Alignment File
        try:
            bamfile = pysam.AlignmentFile(os.path.join(tempfile.gettempdir(), filefront + '.bam'))
            self.bamfiles.append(bamfile)
            print 'loaded file: ' + bamfile.filename
        except:
            flash('Invalid SAM file: ' + filename, 'error')

    def __del__(self):
        for file in self.bamfiles:
            file.close()
        for file in self.filelist:
            os.remove(file)
            print 'deleted file: ' + file


class WriteZip:
    """Class for handling any group of files that needs to be written to .zip as a response."""

    def __init__(self, filename):
        self.filelist = [] # Need to clean this up afterwards
        self.filename = filename

    def add_file(self, filename):
        self.filelist.append(filename)

    def send_zipfile(self):
        memory_file = io.BytesIO()
        with zipfile.ZipFile(memory_file, 'w') as zf:
            files = self.filelist
            for file in files:
                zf.write(os.path.join(tempfile.gettempdir(), file), file)
        memory_file.seek(0)
        return memory_file

    def __del__(self):
        for file in self.filelist:
            os.remove(file)
            print 'deleted file: ' + file


class ParsePsy(HTMLParser):
    """Class for turning the html file in a list of strings"""

    def __init__(self):
        # Pass down parent class HTMLParser
        HTMLParser.__init__(self)
        # Initialize data to be extracted from HTML file
        self.data = []

    def handle_data(self,data):
        if data.strip():
            self.data.append(data.replace('\r\n ', ''))


def random_id(size=6, chars=string.ascii_uppercase + string.digits):
    """Generate a random string of letters and digits of a given length.

    Parameters:
        size (int): Length of the random string to be generated.
        cahrs (list): List of valid characters to be used in generation.

    Returns:
        A string of random characters of length `size`.

    """
    return ''.join(random.choice(chars) for _ in range(size))
