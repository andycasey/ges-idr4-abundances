#!/usr/bin/python

"""
Remove duplicate lines in all the individual line abundance files.
"""

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import os
from glob import glob

filenames = glob("*/*.dat")
N = len(filenames)

for i, filename in enumerate(filenames):
    print("{0}/{1}: {2}".format(i, N, filename))

    os.system("uniq {0} > a".format(filename))
    os.system("mv a {0}".format(filename))

