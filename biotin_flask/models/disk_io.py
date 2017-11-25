# -*- coding: utf-8 -*-
"""Utilities for uploading and interacting with files on disk.

This module primarily deals with temporary files that are uploaded by users to
app.config['UPLOAD_FOLER'] - which is currently set to biotin_flask/tmp/

Since the server hard drive will eventually run out of space if users are
allowed to upload files without restraint, this module implements a threading
mechanism to automatically delete files after a certain age. This age is set
in: app.config['MAX_FILE_AGE'].

"""

import threading
import time
import os

import pysam
from werkzeug import secure_filename

from biotin_flask import app

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'


###############################################################################
# Uploading Files to Disk
###############################################################################

def upload(file_h, filename, folder=None, prefix=None,
           timeout=app.config['MAX_FILE_AGE']):
    """Save a file in app.config['UPLOAD_FOLDER'] and delete after some time.

    Parameter:
        file_h (Flask.Object): The Flask file handle.
        filename (str): The filename (without path) to be saved.
        folder (str): The subfolder in which to be saved. If None, save in the
                root folder.
        timeout (int): Number of seconds before the file will be deleted.
        prefix (str): Append a prefix to the filename.

    Returns:
        The string full path of the file saved.

    """
    if len(filename) == 0:
        raise Exception("Invalid filename in file upload!")
    if prefix:
        filename = prefix + filename
    filename = secure_filename(filename)
    if folder:
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], folder, filename)
    else:
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file_h.save(filepath)
    if timeout:
        _delete_after_time(filepath, timeout)
    return filepath


def upload_sam(file_h, filename, prefix=None):
    """Upload a SAM or BAM file, then sort/index it on disk.

    If the uploaded file is a .sam file, delete it after creating the .bam
    and .bam.bai file to save space.

    Parameter:
        file_h (Flask.Object): The Flask file handle.
        filename (str): The filename (with extension) to be uploaded.
        prefix (str): Append a prefix to the filename.

    Returns:
        The string full path of the bam file saved.

    """
    front, ext = os.path.splitext(filename)
    if not ext.lower().endswith(('.sam', '.bam')):
        raise Exception("Invalid file extension in SAM file upload!")

    # Input file is a SAM file
    if ext.lower() == '.sam':
        samname = front + '.sam'
        sampath = upload(file_h, samname, 'sam', prefix, timeout=None)

        # Sort the extracted SAM file into a compressed BAM file
        bampath = os.path.splitext(sampath)[0] + '.bam'
        pysam.sort('-o', bampath, sampath)
        os.remove(sampath)
        _delete_after_time(bampath)

    # Input file is a BAM file
    else:
        bamname = front + '.bam'
        bampath = upload(file_h, bamname, 'sam', prefix)

    # Index the generated BAM file
    pysam.index(bampath)
    _delete_after_time(bampath + '.bai')


###############################################################################
# Deleting Files from Disk
###############################################################################

def delete_file(folder, filename, prefix=None):
    """Delete the requested filename in app.config['UPLOAD_FOLDER'].

    Parameters:
        folder (str): The subdirectory in which to search. If None, search from
                the root directory.
        filename (ext): The filename to delete.
        prefix (str): A string prefix to match.

    Returns:
        True if any file was deleted, otherwise False.

    """
    if prefix:
        filename = secure_filename(prefix + filename)
    else:
        filename = secure_filename(filename)
    if folder:
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], folder, filename)
    else:
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    if os.path.isfile(filepath):
        os.remove(filepath)
        return True
    return False


def delete_files(folder=None, prefix=None, filefront=None, ext=None):
    """Delete all files in app.config['UPLOAD_FOLDER'] matching criteria.

    Parameters:
        folder (str): The subdirectory to delete files from. If None, delete
                from the root directory.
        prefix (str): A string prefix to match.
        filefront (str): The filename (with no extension) to match.
        ext (tuple): A list of extensions to match. If None, list everything.

    Returns:
        True if any files were deleted, otherwise False.

    """
    # Determine common startswith match string
    if prefix and filefront:
        filestart = secure_filename(prefix + filefront) + '.'
    elif prefix:
        filestart = secure_filename(prefix)
    elif filefront:
        filestart = secure_filename(filefront) + '.'
    else:
        filestart = None

    # List files in a directory
    if folder:
        files = os.listdir(os.path.join(app.config['UPLOAD_FOLDER'], folder))
    else:
        files = os.listdir(app.config['UPLOAD_FOLDER'])
    # Don't delete the .gitignore
    if '.gitignore' in files:
        files.remove('.gitignore')

    # Match files that begin with a certain string
    if filestart:
        result = []
        for f in files:
            if f.startswith(filestart):
                result.append(f)
        files = result

    # Match files that end with an extension
    if ext:
        result = []
        for f in files:
            if f.endswith(ext):
                result.append(f)
        files = result

    # Remove files
    if files:
        for f in files:
            if folder:
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], folder, f)
            else:
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], f)
            os.remove(filepath)
        return True
    return False


###############################################################################
# Finding and Searching Files on Disk
###############################################################################

def list_dir(folder=None, prefix=None, ext=None, remove_prefix=True):
    """List all files in app.config['UPLOAD_FOLDER'] matching criteria.

    Parameters:
        folder (str): The subdirectory to list files from. If None, list from
                the root directory.
        prefix (str): A string prefix to match.
        ext (tuple): A list of extensions to match. If None, list everything.
        remove_prefix (bool): If true, remove the prefix string from results.

    Returns:
        A list of string filenames matching the criteria.

    """
    if folder:
        files = os.listdir(os.path.join(app.config['UPLOAD_FOLDER'], folder))
    else:
        files = os.listdir(app.config['UPLOAD_FOLDER'])
    result = []

    # Query using both prefix and ext
    if prefix and ext:
        for f in files:
            if f.startswith(prefix) and f.endswith(ext):
                if remove_prefix:
                    result.append(f[len(prefix):])
                else:
                    result.append(f)

    # Query using only prefix
    elif prefix:
        for f in files:
            if f.startswith(prefix):
                if remove_prefix:
                    result.append(f[len(prefix):])
                else:
                    result.append(f)

    # Query using only ext
    elif ext:
        for f in files:
            if f.endswith(ext):
                result.append(f)

    # No query provided
    else:
        return files

    return result


def get_firstfile(folder=None, prefix=None, ext=None, remove_prefix=True):
    """Get the first file in app.config['UPLOAD_FOLDER'] matching criteria.

    Parameters:
        folder (str): The subdirectory to list files from. If None, list from
                the root directory.
        prefix (str): A string prefix to match.
        ext (tuple): A list of extensions to match. If None, list everything.
        remove_prefix (bool): If true, remove the prefix string from results.

    Returns:
        A string filename of the first file matching the criteria.

    """
    files = list_dir(folder, prefix, ext, remove_prefix)
    if files:
        return files[0]
    return None


###############################################################################
# Private Methods
###############################################################################

def _sleeper_delete(filepath, timeout):
    """Delete a specified file after a certain amount of time.

    This subroutine is only meant to be called by the_delete_after_time method
    that spawns a new thread.

    Parameter:
        filepath (str): The filename (with path) to be deleted.
        time (int): Number of seconds before the file will be deleted.

    """
    time.sleep(timeout)
    if os.path.isfile(filepath):
        os.remove(filepath)


def _delete_after_time(filepath, timeout=app.config['MAX_FILE_AGE']):
    """Spawns a sleeper thread to delete a file after some time.

    Parameter:
        filepath (str): The filename (with path) to be deleted.
        time (int): Number of seconds before the file will be deleted.

    """
    t = threading.Thread(target=_sleeper_delete, args=(filepath, timeout))
    t.start()
