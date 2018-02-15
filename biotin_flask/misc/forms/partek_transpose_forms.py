# -*- coding: utf-8 -*-
"""Partek transpose forms."""

from flask_wtf import FlaskForm
from flask_wtf.file import FileAllowed, FileField, FileRequired

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'

required = FileRequired('This field is required.')


class PartekTransposeForm(FlaskForm):
    """Flask-WTF form for the Partek Transpose View.

    Attributes:
        csv (File): Required. Filename must end in `.csv`

    """

    csv = FileField('File (.csv)', validators=[
        FileRequired('This field is required.'),
        FileAllowed(['csv'], 'Only (.csv) files are accepted.')
    ])
