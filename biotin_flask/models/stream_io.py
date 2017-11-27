# -*- coding: utf-8 -*-
"""Utilities for sending files to the user.

This module primarily deals with files that will not reside on the server
hard drive. Any saved files should be deleted immediately after they are
served to the user.

"""

import csv
import io
import os
import StringIO
import tempfile
import zipfile

from flask import make_response, send_file

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Production'


###############################################################################
# Streams - These files are never saved to the disk
###############################################################################

def csv_response(filename, data, delimiter=','):
    """Return a .csv file as a Flask response.

    Parameters:
        filename (str): The filename to send back to the user.
        data (list): Data to write to csv (see example).
        delimiter (str): The character to use as the csv delimiter.

        data_input_example = [
            ['Hello', 'World', 'This', 'is', 'first', 'row'],
            ['This', 'Is', 'The', 'Second', 'Row'],
            ['Third', 'Row', 'Here']
        ]

    Returns:
        A Flask response (attachment).

    """
    # Write CSV file
    destination = StringIO.StringIO()
    writer = csv.writer(destination, delimiter=delimiter)
    for row in data:
        writer.writerow(row)

    # Make Flask Response
    result = make_response(destination.getvalue())
    result.headers['Content-Disposition'] = 'attachment; filename=' + filename
    return result


def zip_response(filename, filepaths, filenames):
    """Return a .zip file as a Flask response.

    Parameters:
        filename (str): The filename to send back to the user.
        filepaths (list): A string list of file paths to include in the .zip.

    """
    memory_file = io.BytesIO()
    with zipfile.ZipFile(memory_file, 'w') as file_h:
        i = 0
        for filepath in filepaths:
            file_h.write(filepath, filenames[i])
            i = i + 1
    memory_file.seek(0)
    return send_file(
        memory_file,
        attachment_filename=filename,
        as_attachment=True
    )


###############################################################################
# Object Garbarge Removal - These files are cleaned up on object destruction
###############################################################################

def zip_csv_response(filename, data, delimiter=','):
    """Return a .zip file of multiple csv files as a Flask response.

    Parameters:
        filename (str): The filename to send back to the user.
        data (dict): Data to write to csv (see example).
        delimiter (str): The character to use as the csv delimiter.

        data_input_example = {
            'file1.csv': [
                ['Hello', 'World', 'Blah'],
                ['This', 'Is', 'The', 'Second', 'Row'],
                ['Third', 'Row', 'Here']
            ],
            'file2.csv': [
                ['Hello', 'World', 'Blah'],
                ['This', 'Is', 'The', 'Second', 'Row'],
                ['Third', 'Row', 'Here']
            ]
        }

    Returns:
        A Flask response (attachment).

    """
    temp_manager = TemporaryFileManager()

    # Write CSV files
    csvnames = []
    for csvname, csvdata in data.items():
        csvpath = os.path.join(tempfile.gettempdir(), csvname)
        with open(csvpath, 'wb') as file_h:
            w = csv.writer(file_h, delimiter=delimiter)
            for row in csvdata:
                w.writerow(row)
        temp_manager.add_file(csvpath)
        csvnames.append(csvname)

    # Write and send ZIP file
    return zip_response(filename, temp_manager.filelist, csvnames)


###############################################################################
# Temporary File Class
###############################################################################

class TemporaryFileManager:
    """Class that will delete files after the instance falls out of scope.

    Attributes:
        filelist (list): A string list of paths to clean up later.

    """

    def __init__(self, filepath=None):
        """Construct a TemporaryFileManager object."""
        if filepath:
            self.filelist = [filepath]
        else:
            self.filelist = []

    def add_file(self, filepath):
        """Add a new file to the temporary file manager."""
        self.filelist.append(filepath)

    def __del__(self):
        """Destruct the TemporaryFileManager object and delete files."""
        for filepath in self.filelist:
            if os.path.isfile(filepath):
                os.remove(file)
