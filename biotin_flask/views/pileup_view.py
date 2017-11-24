# -*- coding: utf-8 -*-
"""Visualize each read within a SAM file at the given coordinates.

This module is a simple wrapper around the Pysam package that allows for
selected nucleotide positions of a SAM file to visualized on an HTML page.

Either a single SAM file or multiple SAM files within a .zip folder can be
provided. If a single SAM file is provided, the results are printed on the HTML
page. If multiple files are provided, the results are returned as a zipped csv.

Several options are available:
    ShowIndel   Show aligner insertions and deletions on the printed output.
    ShowExt     Show nucleotide positions on reads extending outside of the
                bounds of the specified coordinates.
    PrintCsv    Print the results as a CSV file instead of HTML.
"""

import os
import csv
import tempfile
import StringIO

from flask import render_template, request, flash, send_file, make_response
from werkzeug import secure_filename

from biotin_flask import app
from biotin_flask.models.utils import SamUpload, WriteZip, FastaUpload
from biotin_flask.models.pysam_ext import get_refbase

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Production'


@app.route('/sam/pileup', methods=['GET', 'POST'])
def pileup():
    """Handle GET or POST requests to the pileup URL.

    Form args:
        region:     (str) Coordinate region to iterate over.
        sam:        (File) SAM or ZIP file to process.
        ids:        (str) (Optional) Line-delimited coordinates to display.
        fasta:      (File) (Optional) FASTA file to use for genomic bases.
        options:    (list of str) (Optional) List of options.

    Returns:
        A Flask Response object.

    """
    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('pileup/form.html')

    # Otherwise validate the form on a POST request
    region = request.form.get('region').encode('ascii', 'ignore')
    ids = request.form.get('ids')
    sam_upload = request.files['sam']
    fasta_upload = request.files['fasta']
    SHOWINDEL = PRINTCSV = FASTA = False
    FILTER = TRUNCATE = True
    for option in request.form.getlist('options'):
        if option == 'ShowIndel': SHOWINDEL = True
        if option == 'ShowExt': TRUNCATE = False
        if option == 'PrintCsv': PRINTCSV = True
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
    if ("-" not in region) or (":" not in region):
        flash('Invalid amplicon region provided.', 'error')
        return render_template('pileup/form.html')
    genomic_seq = None
    if fasta_upload:
        FASTA = True
        genomic_seq = FastaUpload.extract_seq(fasta_upload, region.split(":")[0])

    # Now load samfile
    bam = SamUpload(sam_upload, secure_filename(sam_upload.filename))
    zip_writer = WriteZip('results_' + sam_upload.filename)

    # Begin loop
    for samfile in bam.bamfiles:

        # Get reference positions list
        ref_list = []
        try:
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
        except ValueError:
            flash(region + ': Amplicon region gene/chromosome name not found in the SAM file, or invalid coordinate range.', 'error')
            return render_template('pileup/form.html')

        # Use a list of lists to hold the bases and reads
        rows = []
        for read in samfile.fetch(region=region):
            base_list = []
            if read.is_reverse:
                base_list.append("-")
            else:
                base_list.append("+")
            base_list.append(read.query_name)
            for pos in ref_list:
                base_list.append(get_refbase(read, pos))
            rows.append(base_list)

        # Make list of FASTA IDs
        new_list = []
        new_list.insert(0, "FASTA ID:")
        new_list.insert(0, "(+/-)")
        for num in ref_list:
            if float(num).is_integer():
                new_list.append(num + 1)
            else:
                new_list.append("-")

        # Make genomic sequence header if applicable
        gen_list = []
        gen_list.insert(0, region.split(":")[0])
        gen_list.insert(0, "")
        if FASTA:
            for i in new_list[2:]:
                try:
                    num = int(i)
                    gen_list.append(genomic_seq[num - 1])
                except:
                    gen_list.append("-")

        # If there is only a single file, return an HTML page or print a single csv file
        if len(bam.bamfiles) == 1:

            if PRINTCSV:

                # Prepare csv printer
                dest = StringIO.StringIO()
                writer = csv.writer(dest)

                # Print rows
                if FASTA: writer.writerow(gen_list)
                writer.writerow(new_list)
                for row in rows:
                    writer.writerow(row)

                # Make response
                response = make_response(dest.getvalue())
                csvname = os.path.splitext(os.path.split(samfile.filename)[1])[0] + '.csv'
                response.headers["Content-Disposition"] = "attachment; filename=results_" + csvname
                return response

            else:

                # Print HTML file
                if FASTA:
                    return render_template('pileup/results.html', rows=rows, sites=new_list, genomic=gen_list)
                return render_template('pileup/results.html', rows=rows, sites=new_list)

        else:
            # Save a temporary csv
            csvname = os.path.splitext(os.path.split(samfile.filename)[1])[0] + '.csv'
            with open(os.path.join(tempfile.gettempdir(), csvname), 'wb') as csvfile:
                w = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                if FASTA: w.writerow(gen_list)
                w.writerow(new_list)
                for row in rows:
                    w.writerow(row)
            zip_writer.add_file(csvname)

    # Send zipfile if applicable
    return send_file(zip_writer.send_zipfile(), attachment_filename=zip_writer.filename, as_attachment=True)
