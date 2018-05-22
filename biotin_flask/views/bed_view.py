"""This module grabs SNP Information from dbSNP. It uses rs numbers
to acquire information about global coordinates,
variant alleles, and the identity of the base upstream of the SNP.
The input is a csv file that is a list of rs id's.
The output is a csv file.
"""

import os, io, csv, StringIO
from flask import render_template, request, flash, make_response
from werkzeug import secure_filename
from Bio import Entrez
from xml.dom import minidom

from biotin_flask import app

__author__ = 'Eric Zhou'
__copyright__ = 'Copyright (C) 2018, EpigenDx Inc.'
__credits__ = ['Eric Zhou']
__version__ = '0.0.1'
__status__ = 'Production'

Entrez.email = 'eyzhou@college.harvard.edu'

@app.route('/misc/bed', methods=['GET', 'POST'])
def bed():
    """Handle GET or POST requests to the psq URL.

    Form args:
        f:      (File) a CSV uploaded by user

    Returns:
        A CSV file.

    """

    # Render an empty form for GET request
    if request.method == 'GET':
        return render_template('bed/form.html')

    # Otherwise validate the form on a POST request and process uploaded files
    f = request.files.getlist('csv')

    if not f[0]:
        flash('A required field is missing', 'alert-warning')
        return render_template('bed/form.html')

    if len(f) > 1:
        flash('Only one file input is accepted', 'alert-warning')
        return render_template('bed/form.html')

    for file in f:
        filename = secure_filename(file.filename)
        filefront, extension = os.path.splitext(filename)
        if not extension == '.csv':
            flash('Only .csv files are allowed.', 'alert-warning')
            return render_template('bed/form.html')

    try:
        stream = io.StringIO(f[0].stream.read().decode('utf-8-sig'), newline=None)
        csv_reader = csv.reader(stream)
    except:
        flash('Unable to read csv file.', 'alert-warning')
        return render_template('bed/form.html')

    # Header row
    track = str(request.form.get('track'))
    description = str(request.form.get('description'))

    # Begin main processing loop
    data = []
    for row in csv_reader:
        rsid = row[0]
        try:
            response = Entrez.efetch(db='SNP', id=rsid[2:], rettype='flt', retmode='xml')
        except:
            flash('NCBI was unable to fetch one of the rs numbers. Either NCBI was overloaded,'
                  'or one of the rs numbers was invalid.', 'alert-warning')

            # write an error row
            row = ['error','','',rsid,'','','','']
            data.append(row)
            continue

        doc = minidom.parseString(response.read())

        root = doc.getElementsByTagName('Rs')[0]

        # Find out what strand it is
        ass = root.getElementsByTagName('Assembly')[0]
        map = ass.getElementsByTagName('MapLoc')[0]
        strand = map.getAttribute('orient')

        # Get the anchor nucleotide
        # What is one basepair upstream of the SNP
        anc = ''
        for element in root.getElementsByTagName('Ss'):
            if element.getAttribute('orient') == strand:
                anc = element.getElementsByTagName('Seq5')[0].firstChild.nodeValue[-1]
                anc = anc.upper()
                break

        # Get the flanking sequences
        seq = root.getElementsByTagName('Sequence')[0]
        seq5 = seq.getElementsByTagName('Seq5')[0].firstChild.nodeValue
        seq3 = seq.getElementsByTagName('Seq3')[0].firstChild.nodeValue

        # Get the chromosome number
        comp = ass.getElementsByTagName('Component')[0]
        chr = comp.attributes['chromosome'].value

        # Get the HGVS ID's
        nc = []
        nc_ver = []
        for hgvs in root.getElementsByTagName('hgvs'):
            id = hgvs.firstChild.nodeValue
            if id[:2] == 'CM':
                # Get the hg38 coordinates here
                # I think the CM numbers (GenBank accession numbers)
                # are more reliable for the coordinates
                descr = str(id.split(':g.')[1])
                coord = int(''.join(filter(str.isdigit, descr)))
                if '_' in descr:
                    coord = int(descr.split('_')[0])
            elif id[:3] == 'NC_':
                nc.append(id)
                nc_ver.append(int(id.split(".")[1].replace(':g','')))
            else:
                break

        try:
            max_ver = max(nc_ver)
            for entry in nc:
                row = []
                ver = int(entry.split('.')[1].replace(':g',''))
                if ver == max_ver:
                    row.extend(['chr' + str(chr),coord, coord + 1,rsid])
                    descr = str(entry.split(':g.')[1])
                    # If it's a SNP
                    try:
                        ref = descr.split('>')[0][-1]
                        obs = descr.split('>')[1]
                        row.extend(['REF=' + ref + ';OBS=' + obs + ';ANCHOR=' + anc, '.', seq5, seq3])
                        data.append(row)
                    except:
                        # If it's an insertion
                        try:
                            ref = ''
                            obs = descr.split('ins')[1]
                            row.extend(['REF=' + ref + ';OBS=' + obs + ';ANCHOR=' + anc, '.', seq5, seq3])
                            data.append(row)
                        except:
                            # If it's not a SNP or insertion
                            # (if it's a deletion or wildtype allele)
                            if '=' not in descr:
                                row.append('.')
                                row.extend(['.',seq5, seq3])
                                data.append(row)
                            else:
                                pass

        # If there is no NC_ number, then get the CM number
        # The code is mostly copy-pasted from above
        # I was too lazy to make another function...
        except:
            for hgvs in root.getElementsByTagName('hgvs'):
                row = []
                id = hgvs.firstChild.nodeValue
                if id[:2] == 'CM':
                    row.extend(['chr' + str(chr), coord, coord + 1, rsid])
                    # If it's a SNP
                    try:
                        ref = descr.split('>')[0][-1]
                        obs = descr.split('>')[1]
                        row.extend(['REF=' + ref + ';OBS=' + obs + ';ANCHOR=' + anc,'.',seq5, seq3])
                        data.append(row)
                    except:
                        # If it's an insertion
                        try:
                            ref = ''
                            obs = descr.split('ins')[1]
                            row.extend(['REF=' + ref + ';OBS=' + obs + ';ANCHOR=' + anc, '.', seq5, seq3])
                            data.append(row)
                        except:
                            # If it's not a SNP or insertion
                            # (if it's a deletion or wildtype allele)
                            if '=' not in descr:
                                row.append('.')
                                row.extend(['.',seq5, seq3])
                                data.append(row)
                            else:
                                pass
                else:
                    break


    # Prepare csv printer
    dest = StringIO.StringIO()
    writer = csv.writer(dest)

    # Write the header
    writer.writerow(['track name=\"' + track + '\"','description=\"' + description + '\"',
                     'type=bedDetail','','db=hg38','reference=GRCh38.p2'])
    for row in data:
        writer.writerow(row)

    # Make response
    response = make_response(dest.getvalue())
    response.headers["Content-Disposition"] = "attachment; filename=results_" + filename
    return response