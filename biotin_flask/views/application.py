# -*- coding: utf-8 -*-
"""Home page URL routes for the biotin-flask application.

This module contains the routes for the home page of the app. All of the views
contained within this module are for static HTML pages only.

If you are making a view with dynamic or complicated logic, please make a
new python module within this directory and import it in:
biotin_flask/__init__.py

"""

import logging
import sys

from flask import render_template, session

from biotin_flask import app
from biotin_flask.models.disk_io import list_dir

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.2'
__status__ = 'Production'


app.logger.addHandler(logging.StreamHandler(sys.stdout))
app.logger.setLevel(logging.ERROR)


@app.route('/')
def index():
    """Render the home page."""
    return render_template('index.html')


@app.route('/sam/')
def sam():
    """Render the landing page for SAM file analysis."""
    if 'id' in session and len(list_dir('sam', session['id'], ('.bam'))) > 0:
        # At least one SAM file is uploaded, so show extended menu
        return render_template('sam.html', showHidden=True)
    return render_template('sam.html', showHidden=False)


@app.route('/misc/')
def misc():
    """Render the landing page for miscellaneous analysis."""
    return render_template('misc.html')
