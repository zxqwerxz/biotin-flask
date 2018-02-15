# -*- coding: utf-8 -*-
"""Initialize an instance of the biotin-flask application.

This module contains the factory method for setting up the application. See for
reference: http://flask.pocoo.org/docs/0.12/patterns/appfactories/

To accomplish the following activities, go to the following locations:

    * "App Config" - See `settings.py` in the root folder.
    * "Start Database" - See `manage.py` in the root folder.
    * "Run Server" - See `runserver.sh` in the root folder.
    * "Blueprints" - See http://flask.pocoo.org/docs/0.12/blueprints/

"""

import logging
import sys

from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from settings import config

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou', 'Eric Zhou']
__version__ = '0.2.0'
__status__ = 'Development'

db = SQLAlchemy()


def create_app(mode):
    """Create flask application."""
    app = Flask(__name__)
    app.config.from_object(config[mode])

    # Register extensions
    db.init_app(app)

    # Register logger
    app.logger.addHandler(logging.StreamHandler(sys.stdout))
    app.logger.setLevel(logging.ERROR)

    # Clean temporary file folders
    # from biotin_flask.models.disk_io import delete_files
    # delete_files('sam')

    # Register app blueprints
    from .main import main as main_blueprint
    app.register_blueprint(main_blueprint)

    from .sam import sam as sam_blueprint
    app.register_blueprint(sam_blueprint, url_prefix='/sam')

    from .misc import misc as misc_blueprint
    app.register_blueprint(misc_blueprint, url_prefix='/misc')

    return app
