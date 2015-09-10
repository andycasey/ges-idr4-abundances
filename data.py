#!/usr/bin/python

""" Load the Gaia-ESO Survey iDR4 node abundances. """

from glob import glob
from astropy.io import fits

def get_elements(extension):
    """
    Return the elemental abundance columns available in the extension.
    """

    return { column: comment for _, column, comment in extension.header.cards \
        if "Abundance" in comment }
    


def load_node(filename):

    return fits.open(filename)