"""Miscellaneous common shared utilities."""

import random
import string

from HTMLParser import HTMLParser

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou', 'Eric Zhou']
__version__ = '0.0.6'
__status__ = 'Production'


class ParsePsy(HTMLParser):
    """Class for turning the html file in a list of strings"""

    def __init__(self):
        # Pass down parent class HTMLParser
        HTMLParser.__init__(self)
        # Initialize data to be extracted from HTML file
        self.data = []

    def handle_data(self,data):
        if data.strip():
            self.data.append(data.replace('\r\n ', ''))


def random_id(size=6, chars=string.ascii_uppercase + string.digits):
    """Generate a random string of letters and digits of a given length.

    Parameters:
        size (int): Length of the random string to be generated.
        cahrs (list): List of valid characters to be used in generation.

    Returns:
        A string of random characters of length `size`.

    """
    return ''.join(random.choice(chars) for _ in range(size))
