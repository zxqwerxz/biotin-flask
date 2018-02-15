Biotin-Flask
============

Web Interface for Tools for DNA methylation analysis, developed in Flask.

Development Server Startup Instructions
---------------------------------------

First clone this repository to your local machine.

Make sure your computer has the following installed:
* Python 2.7 (and pip)
* [virtualenv](https://virtualenv.pypa.io/en/stable/)

Second, prepare a virtualenv to hold the biotin-flask dependencies.

```bash
cd biotin-flask
virtualenv env
. env/bin/activate
pip install -r requirements.txt
```

Whenever you wish to run the biotin-flask development server, you need to make
sure the virtualenv is running. (env) should be displayed on the far left of
the terminal prompt.

Then, you can start the development server:

```bash
python runserver.py
```

When you are finished using biotin-flask, you can deactivate the virtualenv:

```bash
deactivate
```

The next time you need to use biotin-flask again, reactivate the virtualenv:

```bash
. env/bin/activate
```

Heroku Instructions
-------------------

The link to the application on Heroku can be [found here](https://biotin-flask.herokuapp.com/).

The heroku application is currently managed by: jyzhou@epigendx.com

The following files in the root directory are used for heroku deployment:
* Procfile
* runserver.py

First log into heroku

```bash
heroku login
```

Make sure that the heroku origin is set, and if not, add it:

```bash
git remote -v
git remote add heroku https://git.heroku.com/biotin-flask.git
heroku git:remote -a biotin-flask
```

Then the application secret key should be set as an environment variable:

```bash
heroku config:set SECRET=SECRET_KEY_GOES_HERE
heroku config
```

Then push the application to heroku master:

```bash
git push heroku master
heroku ps:scale web=1
```


Local Server Install Instructions
---------------------------------

The biotin-flask is also running inside the local network at: 198.168.1.13

The following files in the root directory are used for local deployment:
* wsgi.py
* biotin-flask.conf

### Brand new install ###

First the application secret key should be set as an environment variable.

```bash
export SECRET=SECRET_KEY_GOES_HERE
```

Then, copy the .conf file to /etc/apache2/sites-available

Make sure that the logs folder exists

Make sure the wsgi.py file is executable

Additionally, you must set up a virtualenv

    sudo virtualenv env
    sudo ./env/bin/pip install -r requirements.txt

### After the first time ###

Then, a2dissite whatever is active, and a2ensite biotin-flask.

Sudo service apache2 restart


Dependencies
------------
* pip
  * flask
  * Flask-WTF
  * Flask-SQLAlchemy
  * pysam
  * openpyxl


Credits
-------
* Bootstrap v4 Theme - "Sandstone" by Bootswatch
