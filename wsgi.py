import sys
sys.path.insert(0, "/var/www/biotin-flask/")

activate_this = '/var/www/biotin-flask/env/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))

from biotin_flask import app as application
