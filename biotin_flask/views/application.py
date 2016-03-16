from flask import render_template
from biotin_flask import app

# Go to biotin_flask/__init__.py for application configuration

@app.route('/')
def index():
    return render_template('index.html')
    
@app.route('/sam/')
def sam():
    return render_template('sam.html')

@app.route('/misc/')
def misc():
    return render_template('misc.html')
