# Cython compile instructions

from setuptools import setup

from setuptools import Extension, find_packages
import os
import glob
from sys import platform
import sys
import sysconfig


setup(
    name="TinyPeel",
    version="0.0.1",
    author="Andrew Whalen",
    author_email="awhalen@roslin.ed.ac.uk",
    description="Multilocus and hybrid peeling",
    long_description="This is a package for performing multi-locus and hybrid peeling in large pedigree populations. For more information see Hybrid peeling for fast and accurate calling, phasing, and imputation with sequence data of any coverage in pedigrees (Whalen et al, 2018; doi: 10.1186/s12711-018-0438-2) ",
    long_description_content_type="text/markdown",
    packages=['tinypeel', 'tinypeel.tinyhouse','tinypeel.Peeling'],
    package_dir={'': 'src'},

    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    entry_points = {
    'console_scripts': [
        'TinyPeel=tinypeel.tinypeel:main',
        ],
    },
    install_requires=[
        'numpy',
        'numba',
    ]
)
