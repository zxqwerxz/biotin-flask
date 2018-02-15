# -*- coding: utf-8 -*-
"""Biotin-Flask Application Settings.

This module contains all of the settings for the Biotin-Flask app. The purpose
of this module is to help manage different settings in each environment the app
may run. For instance, an app run on `Development` mode has different database
connections and security considerations than an app run on `Production` mode.

The currently supported modes are as follows:

    * 'default' - Default mode selected if no options provided. (development)
    * 'development' - Use this mode on your local computer.

See for reference: http://flask.pocoo.org/docs/0.12/config/

"""

import os

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'

BASEDIR = os.path.abspath(os.path.dirname(__file__))


class BaseConfig:
    """Base configuration settings applied in all cases."""

    APP_NAME = 'Biotin-Flask'
    SECRET_KEY = os.environ.get('SECRET_KEY', 'SECRET_KEY_NOT_SET')
    UPLOAD_FOLDER = os.path.join(BASEDIR, 'tmp')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    MAX_FILE_AGE = 7200


class DevelopmentConfig(BaseConfig):
    """Configuration settings for a development server."""

    DEBUG = True
    SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(BASEDIR, 'dev.db')


# This dictionary of classes is exported for application configuration
config = {
    'default': DevelopmentConfig,
    'development': DevelopmentConfig
}
