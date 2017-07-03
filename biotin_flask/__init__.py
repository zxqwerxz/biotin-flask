# -*- coding: utf-8 -*-
"""Initialize the biotin-flask application.

This module initializes and configures the flask application. Then, it imports
the views that contain the business logic.
"""

import os
from flask import Flask

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou', 'Eric Zhou']
__version__ = '0.1.1'
__status__ = 'Production'

UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), 'tmp')

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = os.environ.get('SECRET', 'Big Secret!')

import biotin_flask.views.application
import biotin_flask.views.pileup_view
import biotin_flask.views.partek_transpose_view
import biotin_flask.views.psq_view
