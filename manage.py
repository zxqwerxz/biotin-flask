#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Various scripts for managing the biotin-flask app.

This module contains methods for customizing the built in flask-cli interface.
Since Flask v0.11, Flask recommends starting Flask development servers using
the cli (http://flask.pocoo.org/docs/0.12/quickstart/). Example:

    export FLASK_APP=manage.py
    flask run

For more complex apps (including Biotin-Flask), it becomes necessary to take
some additional management steps prior to running the development server. The
`Sam Blueprint` utilizes a SQLAlchemy/sqlite3 database in development mdoe. If
no `*.db` exists in the root directory, then it must be created as follows:

    flask initdb

Both of these steps are automated in: `runserver.sh`. Therefore, if you are
ever unsure of how to start the server, it is always safe to run the bash
script, because all the batteries are include there. Example:

    chmod +x runserver.sh
    ./runserver.sh

See also: http://flask.pocoo.org/docs/0.12/cli/

"""

import os

import click

from biotin_flask import create_app, db

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'

app = create_app(os.getenv('FLASK_CONFIG') or 'default')


@app.cli.command()
def initdb():
    """Drop and initialize the local database."""
    click.echo('Now dropping and initializing the database.')
    db.drop_all()
    db.create_all()
