import os, tempfile, zipfile, pysam
from io import BytesIO
from flask import flash

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
        memory_file = BytesIO()
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