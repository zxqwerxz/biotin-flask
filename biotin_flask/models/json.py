# -*- coding: utf-8 -*-
"""Utilities for JSON and ajax api usages."""

from flask import jsonify

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.1'
__status__ = 'Development'


###############################################################################
# Status codes
###############################################################################

def json_success():
    """Return a JSON success message."""
    return jsonify({
        'status': 200,
        'message': 'Success.'
    })


def json_notfound():
    """Return a JSON not found message."""
    return jsonify({
        'status': 404,
        'error': 'Not found.'
    }), 404


def json_notauth():
    """Return a JSON not authorized message."""
    return jsonify({
        'status': 401,
        'error': 'Access denied.'
    }), 401
