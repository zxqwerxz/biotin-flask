import os, tempfile, pysam
from flask import render_template, request, flash
from werkzeug import secure_filename

from biotin_flask import app
from biotin_flask.models.utils import SamUpload

@app.route('/sam/pileup', methods=['GET', 'POST'])
def pileup():

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('pileup/form.html')

    # Otherwise validate the form on a POST request
    region = request.form.get('region')
    ids = request.form.get('ids')
    sam_upload = request.files['sam']
    fasta_upload = request.files['fasta']
    SHOWINDEL = SHOWEXT = SHOWCALC = False
    FILTER = FASTA = True
    for option in request.form.getlist('options'):
        if option == 'ShowIndel': SHOWINDEL = True
        if option == 'ShowExt': SHOWEXT = True
        if option == 'ShowCalc': SHOWCALC = True
    if not region or not sam_upload:
        flash('A reqired field is missing', 'error')
        return render_template('pileup/form.html')
    if not ids: FILTER = False
    if not fasta_upload: FASTA = False

    # Now load samfile
    bamfiles = SamUpload(sam_upload, secure_filename(sam_upload.filename))

    return render_template('pileup/form.html')
