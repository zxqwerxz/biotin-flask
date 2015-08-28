from flask import render_template, request
from biotin_flask import app
import pysam

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
        
        samfile = pysam.AlignmentFile(file, "r")