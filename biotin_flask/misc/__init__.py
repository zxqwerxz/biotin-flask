"""Setup Miscellaneous Blueprint."""

from flask import Blueprint, render_template

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'

misc = Blueprint('misc', __name__)


@misc.route('/')
def index():
    """Render the landing page for miscellaneous programs."""
    return render_template('misc/index.html')

# Import other views
from .views import partek_transpose_view  # noqa
from .views import psq_view  # noqa
from .views import genotyping_view  # noqa
from .views import cov_view  # noqa
