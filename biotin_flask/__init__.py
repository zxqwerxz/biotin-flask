# -*- coding: utf-8 -*-
"""Initialize the biotin-flask application.

This module initializes and configures the flask application. Then, it imports
the views that contain the business logic.
"""

import logging
import os
import sys
from flask import Flask

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou', 'Eric Zhou']
__version__ = '0.1.1'
__status__ = 'Production'

UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), 'tmp')

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_FILE_AGE'] = 7200
app.secret_key = os.environ.get('SECRET', 'Big Secret!')
app.logger.addHandler(logging.StreamHandler(sys.stdout))
app.logger.setLevel(logging.ERROR)

# Clean temporary file folders
from biotin_flask.models.disk_io import delete_files
delete_files('sam')

import biotin_flask.views.application
import biotin_flask.views.partek_transpose_view
import biotin_flask.views.psq_view
import biotin_flask.views.genotyping_view
import biotin_flask.views.cov_view
import biotin_flask.views.sam
import biotin_flask.views.epic_view
import biotin_flask.views.bed_view
import biotin_flask.views.cpg_view
import biotin_flask.views.tngbs_view
import biotin_flask.views.fasta_view
import biotin_flask.views.coord_view