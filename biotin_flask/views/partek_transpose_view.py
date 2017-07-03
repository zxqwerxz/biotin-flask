import io, os, csv, StringIO
from collections import OrderedDict

from flask import render_template, request, flash, make_response
from werkzeug import secure_filename

from biotin_flask import app

def row_to_dict(csv_row):
    rec = OrderedDict()
    try:
        rec['log-odds ratio of SNP against reference'] = csv_row[1]
        rec['genotype call'] = csv_row[4]
        rec['total Non-Reference bases'] = csv_row[5]
        rec['total coverage at locus'] = csv_row[6]
        rec['non-reference average base qualities'] = csv_row[7]
        rec['reference base qualities'] = csv_row[8]
        rec['non-reference average mapping qualities'] = csv_row[9]
        rec['reference average mapping qualities'] = csv_row[10]
        rec['A'] = csv_row[11]
        rec['C'] = csv_row[12]
        rec['G'] = csv_row[13]
        rec['T'] = csv_row[14]
    except:
        flash('Csv file in inproper format.', 'error')
    return rec

@app.route('/misc/partek_transpose', methods=['GET', 'POST'])
def partek_transpose():

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('partek_transpose/form.html')

    # Otherwise validate the form on a POST request and convert uploaded file into a stream
    f = request.files['csv']
    filename = secure_filename(f.filename)
    if not f:
        flash('A reqired field is missing', 'error')
        return render_template('partek_transpose/form.html')
    filefront, extension = os.path.splitext(filename)
    if not extension == '.csv':
        flash('Only .csv files are allowed.', 'error')
        return render_template('partek_transpose/form.html')
    try:
        stream = io.StringIO(f.stream.read().decode("UTF8"), newline=None)
        csv_reader = csv.reader(stream)
    except:
        flash('Unable to read csv file.', 'error')
        return render_template('partek_transpose/form.html')

    # Begin main processing loop
    pos_list = []
    ref_list = []
    samples = {}
    is_header = True
    for row in csv_reader:

        # Skip first row
        if is_header:
            is_header = False
            continue

        # Initialize variables
        record = row_to_dict(row)
        pos = row[0]
        ref_base = row[3]
        sample_id = row[2]

        # Add to dictionary/lists
        if pos not in pos_list:
            pos_list.append(pos)
            ref_list.append(ref_base)
        if sample_id not in samples:
            samples[sample_id] = {}
            samples[sample_id][pos] = record
        elif pos in samples[sample_id]:
            flash('Duplicate record! Problem in csv file.', 'error')
        else:
            samples[sample_id][pos] = record

    # Prepare csv printer
    dest = StringIO.StringIO()
    writer = csv.writer(dest)

    # Print header
    writer.writerow(["", "Position"] + pos_list)
    writer.writerow(["Sample ID", "reference base"] + ref_list)

    # Loop through Samples
    for sample_id, sample_data in sorted(samples.items()):

        # Iterate through each type of measurement (e.g. log-odds ratio)
        for key in sample_data.values()[0].keys():
            row = []
            if key == 'log-odds ratio of SNP against reference':
                row.append(sample_id)
            else:
                row.append("")
            row.append(key)

            # Iterate through each position
            for pos in pos_list:
                if pos in sample_data:
                    row.append(sample_data[pos][key])
                else:
                    row.append("")
            writer.writerow(row)

    # Make response
    response = make_response(dest.getvalue())
    response.headers["Content-Disposition"] = "attachment; filename=results_" + filename
    return response
