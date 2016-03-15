import os, tempfile, csv
from flask import render_template, request, flash, send_file
from werkzeug import secure_filename

from biotin_flask import app
from biotin_flask.models.utils import SamUpload, WriteZip
from biotin_flask.models.pysam_ext import get_refbase

@app.route('/sam/pileup', methods=['GET', 'POST'])
def pileup():

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('pileup/form.html')

    # Otherwise validate the form on a POST request
    region = request.form.get('region').encode('ascii', 'ignore')
    ids = request.form.get('ids')
    sam_upload = request.files['sam']
    SHOWINDEL = False
    FILTER = TRUNCATE = True
    for option in request.form.getlist('options'):
        if option == 'ShowIndel': SHOWINDEL = True
        if option == 'ShowExt': TRUNCATE = False
    if not region or not sam_upload:
        flash('A reqired field is missing', 'error')
        return render_template('pileup/form.html')
    if not ids: FILTER = False
    fasta_ids = None
    if FILTER:
        try:
            fasta_ids = map(int, ids.splitlines())
            fasta_ids = [x - 1 for x in fasta_ids]
            fasta_ids.sort()
        except:
            flash('Error in the FASTID(s) provided. Place one integer on each line with no extraneous lines.', 'error')
            return render_template('pileup/form.html')

    # Now load samfile
    bam = SamUpload(sam_upload, secure_filename(sam_upload.filename))
    zip_writer = WriteZip('results_' + sam_upload.filename)

    # Begin loop
    for samfile in bam.bamfiles:

        # Get reference positions list
        ref_list = []
        for col in samfile.pileup(region=region, truncate=TRUNCATE):
            if FILTER:
                if col.reference_pos in fasta_ids:
                    ref_list.append(col.reference_pos)
            else:
                ref_list.append(col.reference_pos)
                if SHOWINDEL:
                    maxindel = 0
                    for read in col.pileups:
                        if read.indel > maxindel:
                            maxindel = read.indel
                    for z in xrange(maxindel):
                        ref_list.append(col.reference_pos + 0.01 * (z + 1))

        # Use a list of lists to hold the bases and reads
        rows = []
        for read in samfile.fetch(region=region):
            base_list = []
            base_list.append(read.query_name)
            for pos in ref_list:
                base_list.append(get_refbase(read, pos))
            rows.append(base_list)
        new_list = []
        new_list.insert(0, "")
        for num in ref_list:
            if float(num).is_integer():
                new_list.append(num + 1)
            else:
                new_list.append("-")

        # If there is only a single file, return an HTML page
        if len(bam.bamfiles) == 1:
            return render_template('pileup/results.html', rows=rows, sites=new_list)
        else:
            # Save a temporary csv
            csvname = os.path.splitext(os.path.split(samfile.filename)[1])[0] + '.csv'
            with open(os.path.join(tempfile.gettempdir(), csvname), 'wb') as csvfile:
                w = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                w.writerow(new_list)
                for row in rows:
                    w.writerow(row)
            zip_writer.add_file(csvname)

    # Send zipfile if applicable
    return send_file(zip_writer.send_zipfile(), attachment_filename=zip_writer.filename, as_attachment=True)

