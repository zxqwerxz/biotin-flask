import os, tempfile, pysam
from flask import render_template, request, flash
from werkzeug import secure_filename
from collections import OrderedDict

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
    SHOWINDEL = SHOWCALC = False
    FILTER = TRUNCATE = FASTA = True
    for option in request.form.getlist('options'):
        if option == 'ShowIndel': SHOWINDEL = True
        if option == 'ShowExt': TRUNCATE = False
        if option == 'ShowCalc': SHOWCALC = True
    if not region or not sam_upload:
        flash('A reqired field is missing', 'error')
        return render_template('pileup/form.html')
    if not ids: FILTER = False
    if not fasta_upload: FASTA = False

    # Now load samfile
    bam = SamUpload(sam_upload, secure_filename(sam_upload.filename))

    # Later need to incorporate a loop
    samfile = pysam.AlignmentFile(bam.bamfiles[0], "rb")

    # Get reference positions list
    ref_list = []
    for col in samfile.pileup(region=region, truncate=TRUNCATE):
        ref_list.append(col.reference_pos)

    # Use a dictionary of lists to hold the bases and reads
    read_dict = OrderedDict()
    for read in samfile.fetch(region=region):
        read_dict[read.query_name] = []



    return render_template('pileup/form.html')
