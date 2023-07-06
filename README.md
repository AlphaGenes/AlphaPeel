# AlphaPeel

AlphaPeel is a software package for calling, phasing, and imputing genotype and sequence data in pedigree populations. This program implements single-locus peeling, multi-locus peeling, and hybrid peeling. A complete description of these methods is given in the suggested citation:

Whalen, A, Ros-Freixedes, R, Wilson, DL, Gorjanc, G, Hickey, JM. (2018). Hybrid peeling for fast and accurate calling, phasing, and imputation with sequence data of any coverage in pedigrees. Genetics Selection Evolution; doi: <a href="https://doi.org/10.1186/s12711-018-0438-2">https://doi.org/10.1186/s12711-018-0438-2</a>

## User guide

See https://alphapeel.readthedocs.io/en/latest

## Conditions of use

AlphaPeel is fully and freely available for all use under the MIT License.

## Requirements

* Python 3
* NumPy
* Numba

## Installation

AlphaPeel is available on [PyPI](https://pypi.org/project/AlphaPeel): 

    pip install AlphaPeel

## Distribution

Wheel distribution and source distribution are both available on [PyPI](https://pypi.org/project/AlphaPeel/#files).

## Build wheels

 If you want to build the wheels yourself, we require an installation of Python 3 and you need to fork and clone the repository to your local directory first (see [user guide](https://alphapeel.readthedocs.io/en/latest/contribute.html#fork-the-repository))

 Run the following to build the wheel distribution and the source distribution of the package.

First, create the virtual environment:

    python3 -m venv AlphaPeel_env

Next, activate the environment:

    source AlphaPeel_env/bin/activate

Install build:

    python3 -m pip install --upgrade build

Build the distribution:

    python3 -m build

Now, the distributions of AlphaPeel should be available in ``dist/``.

Finally, deactivate the environment:

    deactivate
