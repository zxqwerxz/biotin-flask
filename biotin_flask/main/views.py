# -*- coding: utf-8 -*-
"""Main Blueprint views."""

from flask import render_template

from . import main

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Production'


@main.route('/')
def index():
    """Render the home page."""
    return render_template('main/index.html')
