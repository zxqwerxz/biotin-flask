"""Setup Sam Blueprint."""

from flask import Blueprint, render_template

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'

sam = Blueprint('sam', __name__)


@sam.route('/')
def index():
    """Render the landing page for the sam application."""
    return render_template('sam/index.html')

# Import other views
from .views import upload  # noqa
# from .views import psq_view  # noqa
# from .views import genotyping_view  # noqa
