# -*- coding: utf-8 -*-
"""Sam upload form."""

import os

from flask_wtf import FlaskForm
from flask_wtf.file import FileAllowed, FileField
from wtforms.validators import Optional
from werkzeug import secure_filename

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'

_sam_field = FileField(
    'SAM File(s) (.sam/.bam)',
    validators=[
        Optional(),
        FileAllowed(['sam', 'bam'], 'Only (.sam/.bam) files are allowed.')
    ],
    description='All files must be aligned to the same reference sequence.'
)

_fasta_field = FileField(
    'FASTA File (.fa/.fasta)',
    validators=[
        Optional(),
        FileAllowed(['fa', 'fasta'], 'Only (.fa/.fasta) files are allowed.')
    ],
    description='The reference sequence the sam/bam files were aligned to.'
)


class UploadForm(FlaskForm):
    """Flask-WTF form for the Upload View.

    Attributes:
        sam (File): Filename must end in `.sam` or `.bam`
        fasta (File): Filename must end in `.fa` or `.fasta`

    Note:
        Flask-WTF cannot validate multiple files. It only validates the first
        file in a multiple file upload. Therefore, it's necessary to write a
        manual validator (is_bad_exts).

    """

    sam = _sam_field
    fasta = _fasta_field

    def validate(self):
        """Override default auto-validation."""
        # Call the base validator
        valid = FlaskForm.validate(self)
        if not valid:
            return False
        # Ensure that both `sam` and `fasta` are not empty.
        if not (self.fasta.data or self.sam.data):
            self.sam.errors.append('Cannot submit an empty form.')
            self.fasta.errors.append('Cannot submit an empty form.')
            return False
        return True

    def validate_sam_exts(self, request):
        """Return true if one of the `sam` filename extensions are bad."""
        sams = request.files.getlist('sam')
        for sam in sams:
            ext = os.path.splitext(secure_filename(sam.filename))[1]
            if ext.lower() not in ('.bam', '.sam'):
                self.sam.errors.append('Only (.sam/.bam) files are allowed.')
                return False
        return True
