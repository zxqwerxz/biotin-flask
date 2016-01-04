import os
from flask import Flask

UPLOAD_FOLDER = '/Users/jeffrey/tmp/'

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'Big Secret!'

import biotin_flask.views
