import os
from flask import Flask
app = Flask(__name__)
app.secret_key = 'Big Secret!'

import biotin_flask.views
