===============
Getting Started
===============

Install |Software|
-----------------

Install via pip
===============

|Software| is available on `PyPI <https://pypi.org/project/AlphaPeel>`_:

.. code-block:: bash

    pip install AlphaPeel

Install locally
===============

You can also install |Software| locally. Below are the local install instructions:

Clone the repository:

.. code-block:: bash

    git clone --recurse-submodules https://github.com/AlphaGenes/AlphaPeel.git

Move to the directory:

.. code-block:: bash

    cd AlphaPeel

Create virtual environment:

.. code-block:: bash

    python3 -m venv AlphaPeel_env

Activate the environment:

.. code-block:: bash

    source AlphaPeel-env/bin/activate

Upgrade pip:

.. code-block:: bash

    python3 -m pip install --upgrade pip

Upgrade build:

.. code-block:: bash

    python3 -m pip install --upgrade build

Build the distributions:

.. code-block:: bash

    python3 -m build

Install the package by using the built wheel distribution:

.. code-block:: bash

    python3 -m pip install dist/alphapeel*.whl

Move to example folder to run an example:

.. code-block:: bash

    cd example

Run the example:

.. code-block:: bash

    bash runScript.sh

Deactivate the environment:

.. code-block:: bash

    deactivate

Graphical representation
------------------------

Still in progress...

An example
----------

Still in progress...

.. |Software| replace:: AlphaPeel