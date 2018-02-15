# -*- coding: utf-8 -*-
"""Upload and manage a SAM file."""

import os

from flask import escape, flash, session, render_template, request
from werkzeug import secure_filename

from .. import sam
from .. forms.upload_form import UploadForm
from ... models.pysam_ext import BamFile

"""
from biotin_flask import app
from biotin_flask.models.utils import random_id
from biotin_flask.models.disk_io import basename, delete_file, delete_files, get_firstfile, list_dir, upload, upload_sam  # noqa
from biotin_flask.models.fasta import Fasta
from biotin_flask.models.json import json_success, json_notfound, json_notauth
"""

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.1.0'
__status__ = 'Development'


@sam.route('/upload', methods=['GET', 'POST'])
def upload():
    """Render the file manager for handling SAM and FASTA files.

    Methods:
        GET: Render the File Management Page with form.
        POST: Process the form upload and redirect user.

    """
    form = UploadForm()
    if form.validate_on_submit() and form.validate_sam_exts(request):
        # Perform extra validation
        pass
    return render_template('sam/upload.html', form=form)

    """

    # Generate a new session ID if it doesn't already exist
    if 'id' not in session:
        session['id'] = random_id()

    # Get all files associated with this session id
    oldbams = list_dir('sam', session['id'], ('.bam'))
    oldfasta = get_firstfile('sam', session['id'], ('.fa', '.fasta'))

    # Render an empty form for a GET request
    if request.method == 'GET':
        return render_form(session['id'], oldbams, oldfasta)

    # Otherwise validate the form on a POST request
    if request.method == 'POST':
        sams = request.files.getlist('sam')
        fasta = request.files['fasta']

        # Throw error if entire form is empty
        if not fasta and not sams[0]:
            flash('No files were selected for upload.', 'alert-warning')
            return render_form(session['id'], oldbams, oldfasta)

        # Throw error if any file does not have the appropriate extension
        abort = False
        for sam in sams:
            if len(sam.filename) > 0:
                ext = os.path.splitext(secure_filename(sam.filename))[1]
                if ext.lower() not in ['.bam', '.sam']:
                    abort = True
        if fasta:
            ext = os.path.splitext(secure_filename(fasta.filename))[1]
            if ext.lower() not in ['.fa', '.fasta']:
                abort = True
        if abort:
            flash('Invalid file types were uploaded', 'alert-warning')
            return render_form(session['id'], oldbams, oldfasta)

        # Upload SAM files
        if sams:
            error_count = 0
            for sam in sams:
                if len(sam.filename) > 0:
                    try:
                        upload_sam(sam, sam.filename, escape(session['id']))
                    except IOError, message:
                        flash(message, 'alert-warning')
                        error_count = error_count + 1
            if len(sams) == error_count:
                return render_form(session['id'], oldbams, oldfasta)

        # Upload Fasta file; first delete old fasta files if they exist
        if fasta:
            delete_files('sam', session['id'], None, ('.fa', '.fasta'))
            fa = upload(fasta, fasta.filename, 'sam', escape(session['id']))
            try:
                Fasta(fa)
            except:
                flash(message, 'alert-warning')
                if not sams:
                    return render_form(session['id'], oldbams, oldfasta)

        flash('Files were successfully uploaded.', 'alert-success')
        return render_template('sam.html', showHidden=True)
        """


# Quick shortcut method so I don't have to rewrite this code
def render_form(session_id, samfiles=None, fastafile=None):
    """Render the form associated with this route."""
    return render_template(
        'sam/upload.html',
        session_id=escape(session_id),
        # samfiles=basename(samfiles, session_id),
        # fastafile=basename(fastafile, session_id)
    )


@sam.route('/sam/delete', methods=['POST'])
def sam_delete():
    """Delete a specified file and its similar extensions via AJAX."""
    """
    if 'id' not in session:
        return json_notauth()

    # Fetch queried filename from JSON request
    filename = request.json['file']

    # Delete a bam file and associated files
    if filename.endswith(('.bam')):
        front = os.path.splitext(filename)[0]
        exts = ('.bam', '.bam.bai', '.sam')
        if delete_files('sam', session['id'], front, exts):
            return json_success()
        return json_notfound()

    # Delete some other kind of file
    if delete_file('sam', filename, session['id']):
        return json_success()
    return json_notfound()
    """
    pass


@sam.route('/sam/delete_allbam', methods=['POST'])
def sam_delete_allbam():
    """Delete all the non-fasta files that were uploaded in this session."""
    """
    if 'id' not in session:
        return json_notauth()
    exts = ('.bam', '.bam.bai', '.sam')
    if delete_files('sam', session['id'], None, exts):
        return json_success()
    return json_notfound()
    """
    pass


@sam.route('/sam/logout', methods=['GET'])
def sam_logout():
    """Log out of the session and delete all session files."""
    """
    if 'id' in session:
        delete_files('sam', session['id'])
        session.pop('id', None)
    flash('Session was reset.', 'success')
    return render_template('sam.html', showHidden=False)
    """
    pass
