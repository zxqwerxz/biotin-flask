#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""All-purpose utility script to start flask server.

Use this script for starting the biotin-flask server for any of these usages:
* Testing and staging the application on a local development machine.
* Deploying the flask server on Heroku.

Examples:
    python runserver.py
    python runserver.py --deploy

"""

import argparse
import os
from biotin_flask import app

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.1.1'
__status__ = 'Production'


def main(deploy):
    """Execute main method.

    Args:
        start_option (str): Option to start the flask server with.

    Returns:
        None

    """
    if deploy:
        port = int(os.environ.get("PORT", 5000))
        app.run(host='0.0.0.0', port=port)
    else:
        app.run(debug=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--deploy',
                        help='Start flask server in production mode.',
                        action='store_true')
    args = parser.parse_args()
    main(args.deploy)
