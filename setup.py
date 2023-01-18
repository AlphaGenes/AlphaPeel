from setuptools import setup, Extension, find_packages

setup(
    name="AlphaPeel",
    version="1.1.1",
    author="Andrew Whalen",
    author_email="awhalen@roslin.ed.ac.uk",
    description="Multilocus and hybrid peeling",
    long_description="This is a package for performing multi-locus and hybrid peeling in large pedigree populations. For more information see Hybrid peeling for fast and accurate calling, phasing, and imputation with sequence data of any coverage in pedigrees (Whalen et al, 2018; doi: 10.1186/s12711-018-0438-2) ",
    long_description_content_type="text/markdown",
    packages=['tinypeel', 'tinypeel.tinyhouse','tinypeel.Peeling'],
    package_dir={'': 'src'},
    license="MIT license",

    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    entry_points = {
    'console_scripts': [
        'AlphaPeel=tinypeel.tinypeel:main',
        ],
    },
    install_requires=[
        'numpy>=1.19',
        'numba>=0.50.0'
    ]
)
