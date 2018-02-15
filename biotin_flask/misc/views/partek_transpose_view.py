# -*- coding: utf-8 -*-
"""Partek transpose view."""

from collections import OrderedDict

from flask import render_template, request, flash
from werkzeug import secure_filename

from .. import misc
from .. forms.partek_transpose_forms import PartekTransposeForm
from ... models.stream_io import csv_receive, csv_response

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.1.0'
__status__ = 'Development'


@misc.route('/partek_transpose', methods=['GET', 'POST'])
def partek_transpose():
    """Render the partek transpose form or process the response."""
    form = PartekTransposeForm()

    # If POST request - Validate form and process data
    if request.method == 'POST' and form.validate_on_submit():
        csv_file = form.csv.data
        csv_filename = secure_filename(csv_file.filename)

        # Try to load the csv file
        try:
            csv_reader = csv_receive(csv_file)
        except:
            flash('Unable to read csv file.', 'alert-warning')
            form.csv.errors.append('Error in the (.csv) file.')
            return render_template('misc/partek_transpose.html', form=form)

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
            record = _row_to_dict(row)
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
                flash('Duplicate record! Csv file problem.', 'alert-warning')
            else:
                samples[sample_id][pos] = record

        # Prepare the header for printing
        data = []
        data.append(['', 'Position'] + pos_list)
        data.append(['Sample ID', 'reference base'] + ref_list)

        # Loop through Samples
        for sample_id, sample_data in sorted(samples.items()):

            # Iterate through each type of measurement (e.g. log-odds ratio)
            for key in sample_data.values()[0].keys():
                row = []
                if key == 'log-odds ratio of SNP against reference':
                    row.append(sample_id)
                else:
                    row.append('')
                row.append(key)

                # Iterate through each position
                for pos in pos_list:
                    if pos in sample_data:
                        row.append(sample_data[pos][key])
                    else:
                        row.append('')
                data.append(row)

        return csv_response('results_' + csv_filename, data)

    # GET request, or the POST validation failed
    return render_template('misc/partek_transpose.html', form=form)


def _row_to_dict(csv_row):
    """Convert a single csv.reader() row into a python dictionary."""
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
        flash('Csv file in inproper format.', 'alert-warning')
    return rec
