# -*- coding: utf-8 -*-
"""Advanced SAM file viewer."""

import os

from flask import escape, jsonify, flash, session, render_template, request
from werkzeug import secure_filename

from biotin_flask import app
from biotin_flask.models.utils import random_id, list_uploaded_files
from biotin_flask.models.utils import list_unique_files


__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.2'
__status__ = 'Development'


@app.route('/sam/upload', methods=['GET', 'POST'])
def sam_upload():
    """Render the file manager for handling SAM and FASTA files."""
    # Generate a new session ID if it doesn't already exist
    if 'id' not in session:
        session['id'] = random_id()

    # Get all files associated with this session id
    samfiles = list_unique_files(session['id'], ('.sam', '.bam'))
    fastafiles = list_unique_files(session['id'], ('.fasta', '.fa'))
    if fastafiles:
        fastafile = fastafiles[0]
    else:
        fastafile = None

    # Render an empty form for a GET request
    if request.method == 'GET':
        return render_template(
            'sam/index.html',
            session_id=escape(session['id']),
            samfiles=samfiles,
            fastafile=fastafile
        )

    # Otherwise validate the form on a POST request
    if request.method == 'POST':
        sams = request.files.getlist('sam')
        fasta = request.files['fasta']

        # Throw error if required field is missing
        if not fasta and not sams[0]:
            flash('No files were selected for upload.', 'error')
            return render_template(
                'sam/index.html',
                session_id=escape(session['id']),
                samfiles=samfiles,
                fastafile=fastafile
            )

        # Throw error if any file does not have the appropriate extension
        abort = False
        for f in sams:
            if f.filename:
                ext = os.path.splitext(secure_filename(f.filename))[1]
                if ext not in ['.bam', '.sam']:
                    abort = True
        if fasta:
            ext = os.path.splitext(secure_filename(fasta.filename))[1]
            if ext not in ['.fa', '.fasta']:
                abort = True
        if abort:
            flash('Invalid file types were uploaded', 'error')
            return render_template(
                'sam/index.html',
                session_id=escape(session['id']),
                samfiles=samfiles,
                fastafile=fastafile
            )

        # Upload SAM files
        for f in sams:
            fname = escape(session['id']) + secure_filename(f.filename)
            f.save(os.path.join(app.config['UPLOAD_FOLDER'], 'sam', fname))

        # Upload Fasta file
        if fasta:
            fname = escape(session['id']) + secure_filename(fasta.filename)
            f.save(os.path.join(app.config['UPLOAD_FOLDER'], 'sam', fname))

        return render_template('sam.html')


@app.route('/sam/delete', methods=['POST'])
def sam_delete():
    """Delete a specified file and its similar extensions."""
    if 'id' not in session:
        return jsonify({
            'status': 401,
            'error': 'Access denied.'
        }), 401
    found = False
    exts = ['.sam', '.bam', '.fa', '.fasta']
    for ext in exts:
        filename = escape(session['id']) + request.json['file'] + ext
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], 'sam', filename)
        if os.path.isfile(filepath):
            os.remove(filepath)
            found = True
    if found:
        return jsonify({
            'status': 200,
            'message': 'File successfully deleted.'
        })
    return jsonify({
        'status': 404,
        'error': 'File not found.'
    }), 404


@app.route('/sam/delete_all', methods=['POST'])
def sam_delete_all():
    """Delete all the non-fasta files that were uploaded in this session."""
    if 'id' not in session:
        return jsonify({
            'status': 401,
            'error': 'Access denied.'
        }), 401
    files = list_uploaded_files('sam', escape(session['id']))
    found = False
    for f in files:
        if f.endswith(('.sam', '.bam')):
            name = escape(session['id']) + f
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], 'sam', name)
            os.remove(filepath)
            found = True
    if found:
        return jsonify({
            'status': 200,
            'message': 'File(s) successfully deleted.'
        })
    return jsonify({
        'status': 404,
        'error': 'No files found.'
    }), 404
