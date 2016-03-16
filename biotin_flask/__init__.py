import os
from flask import Flask

UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), 'tmp')

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'Big Secret!'

import biotin_flask.views.application
import biotin_flask.views.pileup_view
import biotin_flask.views.partek_transpose_view
