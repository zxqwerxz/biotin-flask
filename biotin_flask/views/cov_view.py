"""Matches data from a .cov file to a list of FASTA IDs

This module parses FASTA IDs from NGS data and matches corresponding
percent methylation data to FASTA IDs in a CSV file containing reference FASTA IDs

Either a single csv file or multiple csv files can be provided.
"""

import io, os, csv, StringIO

from flask import render_template, request, flash, send_file, make_response
from werkzeug import secure_filename

from biotin_flask import app
from biotin_flask.models.openpyxl_ext import delete_column, wedge_column_right, csv_to_xlsx, delete_row, surround_border

@app.route('/misc/cov', methods=['GET', 'POST'])
def cov():
    """Handle GET or POST requests to the psq URL.

        Form args:
            f:      (File) list of .csv files uploaded by user
            ref:    (File) .csv file containing list of FASTA IDs uploaded by user

        Returns:
            A csv file
    """

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('cov/form.html')

    # Otherwise validate the form on a POST request and process uploaded files

    # Obtain files from the form
    f = request.files.getlist('input')
    ref = request.files.getlist('reference')

    # Throw error if one of the required files are missing
    if not f[0]:
        flash('A reqired field is missing', 'error')
        return render_template('cov/form.html')
    if not ref[0]:
        flash('A reqired field is missing', 'error')
        return render_template('cov/form.html')

    # Validate NGS Data files
    for file in f:
        filename = secure_filename(file.filename)
        filefront, extension = os.path.splitext(filename)
        if not extension == '.csv':
            flash('Only .csv and .csv files are allowed.', 'error')
            return render_template('cov/form.html')

    # Validate reference file
    if len(ref) > 1:
        flash('Only one reference file should be selected')
        return render_template('cov/form.html')
    filename = secure_filename(ref[0].filename)
    filefront, extension = os.path.splitext(filename)
    if not extension == '.csv':
        flash('Only .csv files are allowed.', 'error')
        return render_template('cov/form.html')

    # Read the reference file
    try:
        stream = io.StringIO(ref[0].stream.read().decode("UTF8"), newline=None)
        ref_reader = csv.reader(stream)
    except:
        flash('Unable to read csv file.', 'error')
        return render_template('cov/form.html')

    fasta_ref = []
    for row in ref_reader:
        for element in row:
            fasta_ref.append(element)

    # Prepare csv printer
    dest = StringIO.StringIO()
    writer = csv.writer(dest)

    # Print header
    writer.writerow(["Sample ID"] + fasta_ref)

    # Loop through samples
    for file in f:
        sample_id = secure_filename(file.filename).split(".")[0]
        try:
            stream = io.StringIO(file.stream.read().decode("UTF8"), newline=None)
            reader = csv.reader(stream)
        except:
            flash('Unable to read csv file.', 'error')
            return render_template('cov/form.html')
        fasta_id = []
        methylation = []
        reads = []
        for row in reader:
            fasta_id.append(row[0] + "." + row[1])
            methylation.append(row[3])
            reads.append(int(row[4])+int(row[5]))
        # Sort the list
        methylation_data = []
        reads_data = []
        for element in fasta_ref:
            try:
                indx = fasta_id.index(element)
                methylation_data.append(methylation[indx])
                reads_data.append(reads[indx])
            except:
                methylation_data.append("n/a")
                reads_data.append("n/a")
        # Write to the output
        writer.writerow([sample_id] + methylation_data)
        writer.writerow(['']+reads_data)

    response = make_response(dest.getvalue())
    response.headers["Content-Disposition"] = "attachment; filename=results.csv"
    return response
