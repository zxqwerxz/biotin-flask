"""Setup Main Blueprint."""

from flask import Blueprint

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'

main = Blueprint('main', __name__)

from . import views  # noqa
