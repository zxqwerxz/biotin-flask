import os, tempfile, zipfile, pysam
from flask import flash

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

            # Extract the zipfile
            with zipfile.ZipFile(os.path.join(tempfile.gettempdir(), filename), "r") as z:
                z.extractall(tempfile.gettempdir())
                for file in z.namelist():
                    self.filelist.append(os.path.join(tempfile.gettempdir(), file))
                    filefront, extension = os.path.splitext(file)
                    if extension != '.sam':
                        continue
                    self.__loadsamfile(file, filefront)

        if extension == '.sam':
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
        for file in self.filelist:
            os.remove(file)
            print 'deleted file: ' + file
