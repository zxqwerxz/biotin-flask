Biotin-Flask
============

Web Interface for Tools for DNA methylation analysis, developed in Flask.

Installation Instructions
-------------------------

Virtualenvs do not work with the biotin server -- this is because the default server and mod_wsgi is UCS2 python, but
the virtualenv installation is UCS4. It just won't work.

Currently I am using 