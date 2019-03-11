# Cython compile instructions

from setuptools import setup

from setuptools import Extension, find_packages
import numpy
import os
import glob
from sys import platform
import sys
import sysconfig
py_modules = ['tinyImpute', 'tinyPeel', 'tinyAssign', 'tinyMgsAssign']

src_modules = []
src_modules += glob.glob(os.path.join('src','General', '*.py'))
src_modules += glob.glob(os.path.join('src','Imputation', '*.py'))
src_modules += glob.glob(os.path.join('src','Assign', '*.py'))
src_modules += glob.glob(os.path.join('src','Peeling', '*.py'))

src_modules = [os.path.splitext(file)[0] for file in src_modules]
py_modules += src_modules

setup(
    name="TinyPeel",
    version="0.0.1",
    author="Andrew Whalen",
    author_email="awhalen@roslin.ed.ac.uk",
    description="Multilocus and hybrid peeling",
    long_description="This is a package for performing multi-locus and hybrid peeling in large pedigree populations. For more information see Hybrid peeling for fast and accurate calling, phasing, and imputation with sequence data of any coverage in pedigrees (Whalen et al, 2018; doi: 10.1186/s12711-018-0438-2) ",
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(),
    package_dir={'': 'src'},
    
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    entry_points = {
    'console_scripts': [
        'TinyPeel=tinyPeel:main',
        ],
    },
    install_requires=[
        'numpy',
        'numba',
    ]
)
