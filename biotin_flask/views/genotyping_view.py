"""Inputs NGS genotyping data and outputs genotyping results

This module parses excel sheets uploaded by the user and reads variant call frequencies per sample.
It then outputs the genotyping results of each sample in the form of an excel sheet

Either a single excel file or multiple excel files can be
provided.
"""

import io, os, csv, StringIO
from collections import OrderedDict

from flask import render_template, request, flash, send_file, make_response
from werkzeug import secure_filename
from openpyxl import load_workbook, Workbook

from biotin_flask import app
from biotin_flask.models.openpyxl_ext import delete_column


@app.route('/misc/genotyping', methods=['GET', 'POST'])
def genotyping():
    """Handle GET or POST requests to the psq URL.

        Form args:
            f:      (File) list of .xlsx files uploaded by user
            ref:    (File) .xlsx file containing list of SNP IDs uploaded by user

        Returns:
            An excel file

        """

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('genotyping/form.html')

    # Otherwise validate the form on a POST request and process uploaded files

    # Obtain files from the form
    f = request.files.getlist('input')
    ref = request.files.getlist('reference')

    # Throw error if one of the required files are missing
    if not f[0] or not ref[0]:
        flash('A reqired field is missing', 'error')
        return render_template('genotyping/form.html')

    # Throw error if any file does not have .xlsx extension
    for file in f:
        filename = secure_filename(file.filename)
        filefront, extension = os.path.splitext(filename)
        if not extension == '.xlsx':
            flash('Only .xlsx files are allowed.', 'error')
            return render_template('genotyping/form.html')
    for file in ref:
        filename = secure_filename(file.filename)
        filefront, extension = os.path.splitext(filename)
        if not extension == '.xlsx':
            flash('Only .xlsx files are allowed.', 'error')
            return render_template('genotyping/form.html')

    # Save the files to the tmp folder
    ref[0].save(os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(ref[0].filename)))
    for file in f:
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

    # Open the ref file and store SNP IDs in a list
    list_snp = load_workbook(os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(ref[0].filename)))
    snp_id = []
    for row in range(1, list_snp.active.max_row + 1):
        snp_id.append(list_snp.active['A{}'.format(row)].value)

    # Create the output excel file
    output = Workbook()

    # Open the input files and begin main processing loop
    file_num = 1
    for file in f:
        file_wb = load_workbook(os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(file.filename)))
        file_ws = file_wb.active

        # Deactivate the filter
        file_ws.auto_filter.ref = None

        # Validate the spreadsheet
        fields = ('Chrom', 'Position', 'Gene', 'Assay ID', 'SNP ID', 'Ref Allele', 'Var Allele',
                  'Allele Call', 'Var Frequency', 'Quality', 'Type', 'Allele Source', 'Coverage',
                  'Coverage+', 'Coverage-', 'Allele Cov', 'Allele Cov+', 'Allele Cov-', 'Strand Bias',
                  'Sample Name', 'Barcode')
        for column in range(1, 21):
            if file_ws.cell(row=1, column=column).value != fields[column-1]:
                delete_files(f, ref)
                flash('File not formatted correctly')
                return render_template('genotyping/form.html')


        # Calculate the genotyping result
        # Sort the worksheet based off SNP ID
        for row in range(2, file_ws.max_row + 1):
            cell = file_ws['V{}'.format(row)] # genotyping result

            # columns that represent various parameters
            var_freq = file_ws['I{}'.format(row)].value
            ref_allele = str(file_ws['F{}'.format(row)].value)
            var_allele = str(file_ws['G{}'.format(row)].value)
            coverage = file_ws['M{}'.format(row)].value

            # Calculate the genotype
            if var_freq < 10:
                cell.value = ref_allele + '/' + ref_allele
            elif var_freq > 90:
                cell.value = var_allele + '/' + var_allele
            else:
                cell.value = ref_allele + '/' + var_allele
            if coverage <= 5:
                cell.value = 'NA'
            elif coverage <= 30:
                cell.value = 'MANUAL'

            # Search SNP_ID in snp_id list and assign index to each row
            index_cell = file_ws['W{}'.format(row)]
            rs_id = str(file_ws['E{}'.format(row)].value)

            try:
                index_cell.value = snp_id.index(rs_id) + 1
            except:
                index_cell.value = 'na'

        # Create new sheet for sorting
        file_processed = file_wb.create_sheet('processed')

        # Copy first row to the new sheet
        for column in range(1, file_ws.max_column + 1):
            file_processed.cell(row=1, column=column).value = file_ws.cell(row=1, column=column).value
        file_processed.cell(row=1, column=file_processed.max_column+1).value = 'Genotyping'

        # Sort rows into the new processed sheet
        sort_row = 1
        for snp in range(1, len(snp_id) + 1):
            for row in range(2, file_ws.max_row + 1):
                if file_ws['W{}'.format(row)].value == snp:
                    sort_row += 1
                    # Copy row into new sheet
                    for column in range(1, file_ws.max_column + 1):
                        file_processed.cell(row=sort_row, column=column).value = file_ws.cell(row=row, column=column).value

        # Delete any unnecessary columns
        for row in range(2, file_processed.max_row + 1):
            file_processed = delete_column(file_processed, 'W') # This is the index column

        # Create a new worksheet in output
        output_sheet = output.create_sheet(secure_filename(file.filename))
        file_num += 1

        # Copy the processed worksheet into output
        for row in range(1, file_processed.max_row + 1):
            for column in range(1, file_processed.max_column): # not including the last emptied column
                output_sheet.cell(row=row, column=column).value = file_processed.cell(row=row, column=column).value

    # Delete the first sheet of the output file
    output.remove_sheet(output.active)

    # Save the output file
    if (os.path.isfile(os.path.join(app.config['UPLOAD_FOLDER'], 'results_genotyping.xlsx'))):
        os.remove(os.path.join(app.config['UPLOAD_FOLDER'], 'results_genotyping.xlsx'))
    output.save(os.path.join(app.config['UPLOAD_FOLDER'], 'results_genotyping.xlsx'))

    # Delete files from server
    delete_files(f, ref)

    return send_file(os.path.join(app.config['UPLOAD_FOLDER'], 'results_genotyping.xlsx'), attachment_filename='results.xlsx', as_attachment=True)


def delete_files(f, ref):
    """

    :param f:       (File) list of files
    :param ref:     (File) list of files
    :return:        True if successful, False is unsuccessful
    """

    try:
        for file in f:
            filename = secure_filename(file.filename)
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        for file in ref:
            filename = secure_filename(file.filename)
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        return True
    except:
        return False