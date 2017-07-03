#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Ubuntu-only utility script to start flask server.

Use this script only for starting the server on a LAMP stack (Ubuntu server).
As of 2017, this is only used for the server inside the local network.
"""

import sys

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Production'

sys.path.insert(0, "/var/www/biotin-flask/")
activate_this = '/var/www/biotin-flask/env/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))

from biotin_flask import app as application
