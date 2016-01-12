import os, tempfile, pysam
from flask import render_template, request, flash
from werkzeug import secure_filename
from biotools import samtools

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

    """
    filename = secure_filename(file.filename)

    # Validate data
    if len(gene) == 0 or len(ids) == 0:
        flash('Please do not leave any field blank.', 'error')
        return render_template('allele_specific/single_form.html')
    if '.' not in filename or not (filename.rsplit('.', 1)[1] == 'sam' or filename.rsplit('.', 1)[1] == 'bam'):
        flash('Invalid filename.', 'error')
        return render_template('allele_specific/single_form.html')

    # Extract fasta ids
    try:
        ids_list = map(int, ids.splitlines())
        ids_list.sort()
    except:
        flash('Error in the FASTID(s) provided. Place one integer on each line with no extraneous lines.', 'error')
        return render_template('allele_specific/single_form.html')

    # Process sam file
    try:
        file.save(os.path.join(tempfile.gettempdir(), filename))
        samfile = pysam.AlignmentFile(os.path.join(tempfile.gettempdir(), filename))
        flash('Upload successful!')
    except:
        flash('Invalid SAM file.', 'error')
        return render_template('allele_specific/single_form.html')
    readcount = 0
    readlist = []
    for read in samfile:
        if samfile.getrname(read.tid) == gene:
            readcount += 1
            if read.reference_end >= ids_list[0] and read.reference_start <= ids_list[-1]:
                        readlist.append(read)
    if readcount == 0:
        flash('Specified gene (case-sensitive) was not located in the SAM file.', 'error')
        return render_template('allele_specific/single_form.html')

    # Start base by base printing
    rows = []
    for read in readlist:
        row = []
        inrange = False
        for pos in ids_list:
            base = samtools.getBase(read, pos)
            row.append(base)
            if len(base) > 0:
                inrange = True
        if inrange:
            rows.append(row)

    flash(str(len(rows)) + '/' + str(readcount) + ' (' + str(len(rows)*100/readcount) + '%) of reads aligned to ' + gene + ' are covered at selected sites.')
    return render_template('allele_specific/single_results.html', rows=rows, sites=ids_list)
    """
