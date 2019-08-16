"""Matches EPIC array cg IDs to EpigenDx cg #s
Takes a pre-formatted CSV file as input
And outputs a list of matched EpigenDX cg #s with
EPIC array cg IDs
"""

import os, io, csv, StringIO
from flask import render_template, request, flash, make_response
from werkzeug import secure_filename

from biotin_flask import app

__author__ = 'Eric Zhou'
__copyright__ = 'Copyright (C) 2018, EpigenDx Inc.'
__credits__ = ['Eric Zhou']
__version__ = '0.0.1'
__status__ = 'Production'

@app.route('/misc/epic', methods=['GET', 'POST'])
def epic():

    """Handle GET or POST requests to the psq URL.

    Form args:
        f:      (File) a CSV uploaded by user

    Returns:
        A CSV file.

    """

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('epic/form.html')

    # Otherwise validate the form on a POST request and process uploaded files
    f = request.files.getlist('csv')
    try:
        if request.form['strand'] == 'forward':
            strand = 0
        elif request.form['strand'] == 'reverse':
            strand = 1
    except:
        flash('A required field is missing', 'alert-warning')
        return render_template('epic/form.html')

    if not f[0]:
        flash('A required field is missing', 'alert-warning')
        return render_template('epic/form.html')

    if len(f) > 1:
        flash('Only one file input is accepted', 'alert-warning')
        return render_template('epic/form.html')

    for file in f:
        filename = secure_filename(file.filename)
        filefront, extension = os.path.splitext(filename)
        if not extension == '.csv':
            flash('Only .csv files are allowed.', 'alert-warning')
            return render_template('epic/form.html')

    try:
        stream = io.StringIO(f[0].stream.read().decode('utf-8-sig'), newline=None)
        csv_reader = csv.reader(stream)
    except:
        flash('Unable to read csv file.', 'alert-warning')
        return render_template('epic/form.html')

    # Begin main processing loop
    edx_cg = []
    edx_coord = []
    epic_cg = []
    epic_coord = []
    for row in csv_reader:
        try:
            # Add to lists
            edx_cg.append(row[0])
            edx_coord.append(int(row[1].split(":")[1]))
            if row[2]:
                epic_cg.append(row[2])
                epic_coord.append(int(row[3])+strand)
        except:
            break

    cg = []
    for c, coord in enumerate(epic_coord):
        try:
            index = edx_coord.index(coord)
            cg.append([edx_cg[index], epic_cg[c], edx_coord[index], edx_cg[index] + ' ' + epic_cg[c]])
        except:
            cg.append(['n/a',epic_cg[c],coord,'n/a'])

    for c, cg_num in enumerate(edx_cg):
        try:
            cg[c].append(cg_num)
        except:
            cg.append(['','','','',cg_num])

    # Prepare csv printer
    dest = StringIO.StringIO()
    writer = csv.writer(dest)

    for entry in cg:
        writer.writerow(entry)

    # Make response
    response = make_response(dest.getvalue())
    response.headers["Content-Disposition"] = "attachment; filename=results_" + filename
    return response