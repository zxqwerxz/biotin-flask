from flask import render_template
from biotin_flask import app

@app.route('/')
def index():
    return render_template('index.html')