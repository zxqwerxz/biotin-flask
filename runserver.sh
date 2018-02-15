#!/bin/bash
# File: runserver.sh
# Summary: Starts the Biotin-Flask server in development mode.
# See also: `manage.py` in root directory.
# Instructions:
#     chmod +x runserver.sh
#     ./runserver.sh

# Start the server
export FLASK_APP=manage.py
flask run
