Biotin-Flask
============

Web Interface for Tools for DNA methylation analysis, developed in Flask.

Installation Instructions
-------------------------

Note that the configuration file is biotin_flask/__init__.py

Please configure this file before deploying. 

### Brand new install ###

To install this, copy the .conf file to /etc/apache2/sites-available

Make sure that the logs folder exists

Make sure the wsgi.py file is executable

Additionally, you must set up a virtualenv

    sudo virtualenv env
    sudo ./env/bin/pip install -r requirements.txt

### After the first time ###

Then, a2dissite whatever is active, and a2ensite biotin-flask.

Sudo service apache2 restart

For Heroku
----------

Use default heroku installation instructions.

The Procfile is for Heroku purposes.