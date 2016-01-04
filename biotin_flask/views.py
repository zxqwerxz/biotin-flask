import os, tempfile, csv
from flask import render_template, request, flash, send_file
from werkzeug import secure_filename
from biotin_flask import app
from biotools import samtools
from io import BytesIO
import pysam, zipfile

@app.route('/')
def index():
    return render_template('index.html')
    
@app.route('/sam/')
def sam():
    return render_template('sam.html')
    
@app.route('/sam/upload/', methods=['GET', 'POST'])
def sam_upload():
    if request.method == 'POST':
        # SAM file has been uploaded
        # Initialize variables
        file = request.files['file']
        seqname = request.form['seqname']
        start = int(request.form['range'].split('-')[0])
        end = int(request.form['range'].split('-')[1])
        filename = secure_filename(file.filename)
        file.save(os.path.join(tempfile.gettempdir(), filename))
        
        samfile = pysam.AlignmentFile(os.path.join(tempfile.gettempdir(), filename), "r")
        flash('Upload successful!')

    return render_template('sam.html')

@app.route('/sam/clonal/single', methods=['GET', 'POST'])
def clonal_single():
    error = None
    if request.method == 'POST':
        # SAM file has been uploaded
        # Initialize variables
        file = request.files['file']
        ids = request.form['ids']
        gene = request.form['gene']
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
    return render_template('allele_specific/single_form.html')

@app.route('/sam/clonal/multiple', methods=['GET', 'POST'])
def clonal_multiple():
    error = None
    if request.method == 'POST':
        # SAM file has been uploaded
        # Initialize variables
        zip = request.files['file']
        ids = request.form['ids']
        gene = request.form['gene']
        zipname = secure_filename(zip.filename)

        # Validate data
        if len(gene) == 0 or len(ids) == 0:
            flash('Please do not leave any field blank.', 'error')
            return render_template('allele_specific/multiple_form.html')
        if '.' not in zipname or zipname.rsplit('.', 1)[1] != 'zip':
            flash('Invalid filename.' +zipname, 'error')
            return render_template('allele_specific/multiple_form.html')

        # Extract fasta ids
        try:
            ids_list = map(int, ids.splitlines())
            ids_list.sort()
        except:
            flash('Error in the FASTID(s) provided. Place one integer on each line with no extraneous lines.', 'error')
            return render_template('allele_specific/multiple_form.html')

        # Open zipfile
        csvnames = []
        zip.save(os.path.join(tempfile.gettempdir(), zipname))
        with zipfile.ZipFile(os.path.join(tempfile.gettempdir(), zipname), "r") as z:
            z.extractall(tempfile.gettempdir())
            for filename in z.namelist():
                if '.' not in filename or not (filename.rsplit('.', 1)[1] == 'sam' or filename.rsplit('.', 1)[1] == 'bam'):
                    continue
                try:
                    samfile = pysam.AlignmentFile(os.path.join(tempfile.gettempdir(), filename))
                except:
                    flash('Invalid SAM file: ' + filename, 'error')
                    continue
                readlist = []
                for read in samfile:
                    if samfile.getrname(read.tid) == gene:
                        if read.reference_end >= ids_list[0] and read.reference_start <= ids_list[-1]:
                                    readlist.append(read)
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
                # Save a temporary csv
                csvname = filename.rsplit('.', 1)[0] + '.csv'
                with open(os.path.join(tempfile.gettempdir(), csvname), 'wb') as csvfile:
                    w = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    w.writerow(ids_list)
                    for row in rows:
                        w.writerow(row)
                    csvnames.append(csvname)

        # Make response
        memory_file = BytesIO()
        with zipfile.ZipFile(memory_file, 'w') as zf:
            files = csvnames
            for file in files:
                zf.write(os.path.join(tempfile.gettempdir(), file), file)
        memory_file.seek(0)
        return send_file(memory_file, attachment_filename='results.zip', as_attachment=True)

    return render_template('allele_specific/multiple_form.html')