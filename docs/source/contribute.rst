Contribution and Development Guide
==================================

Welcome to AlphaPeel contribution and development guide.

This guide will give you an overview of the contribution and development workflow
via the AlphaPeel GitHub repository at `<https://github.com/AlphaGenes/AlphaPeel>`_.

Critically, see also a list of issues at `<https://github.com/AlphaGenes/AlphaPeel/issues>`_.

Fork the repository
-------------------

First, you should fork the `repository <https://github.com/AlphaGenes/AlphaPeel>`_.

For more information, see the `GitHub Docs <https://docs.github.com/en/get-started/quickstart/fork-a-repo#forking-a-repository>`_.

Clone your forked repository
----------------------------

Clone your forked repository into a local directory and
initialise submodules at the same time by running the following command in your terminal:

.. code-block:: bash

    git clone --recurse-submodules URL_of_your_forked_repository

Depending on the type of code change,
you should use different branches.
First, check the available branches:

.. code-block:: bash

    cd AlphaPeel
    git branch # check available branches

Large code changes should go to dedicated development branches,
which will be later merged into the ``devel`` branch by maintainers:

.. code-block:: bash

    git checkout devel # start from the development branch
    git branch issue_GitHubIssueNumber # work on issue with GitHub number GitHubIssueNumber
    git checkout issue_GitHubIssueNumber
    # now work on your code changes

It is a good practice that you first open an issue to document what you plan to do,
then follow the above process.

Small code changes can go directly to the ``devel`` branch,
which will eventually be merged into the main branch by maintainers,
but check this with the maintainers when you open/discuss the issue:

.. code-block:: bash

    git checkout devel
    # now work on your code changes

Stable code for wider use will be published is in the ``main`` branch.
While most changes will be happening on the ``devel`` branch,
critical bugfixes, can go to the ``main`` branch,
but check this with the maintainers when you open/discuss the issue:

.. code-block:: bash

    git checkout main
    # now work on your code changes

Make changes in your clone
--------------------------

Make changes to the code and commit them to your local clone repository.
Adding ``AlphaGenes/AlphaPeel#GitHubIssueNumber`` in the message will link the commit with the issue page.

Before you commit the changes,
make sure you test your changes by running the tests and examples.
To this end, you should install ``pytest`` and ``pytest-benchmark``
(see `pytest Documentation <https://docs.pytest.org/en/stable/getting-started.html>`_ and
`pytest-benchmark Documentation <https://pytest-benchmark.readthedocs.io/en/latest/installation.html>`_) and
run ``pytest`` on the distribution built on your modified code to see if the code passes all the tests.

To install ``pytest`` and ``pytest-benchmark``:

.. code-block:: bash

    pip install pytest
    pip install pytest-benchmark


To run ``pytest``:

.. code-block:: bash

    pytest

If the tests run successfully, you are expected to see output similar to the following:

.. code-block::

    ============================= test session starts ==============================
    platform linux -- Python 3.11.14, pytest-9.0.2, pluggy-1.6.0
    benchmark: 5.2.3 (defaults: timer=time.perf_counter disable_gc=False min_rounds=5 min_time=0.000005 max_time=1.0 calibration_precision=10 warmup=False warmup_iterations=100000)
    rootdir: /home/runner/work/AlphaPeel/AlphaPeel
    configfile: pyproject.toml
    plugins: benchmark-5.2.3
    collected 33 items

    tests/accuracy_tests/run_accu_test.py ....................               [ 60%]
    tests/functional_tests/run_func_test.py .............                    [100%]
    ...

The instructions of building your own distribution is available at :ref:`dist-install`.

Instructions on running the examples are at :ref:`run-examples`.

If tests and examples pass, finally install ``pre-commit`` and the ``pre-commit`` hooks for code formatting.

.. code-block:: bash

    pip install pre-commit
    pre-commit install

Later on, the ``pre-commit`` hooks will automatically run on the files you have changed when you commit your changes.

An example output of the ``pre-commit`` hooks is as follows:

.. code-block:: bash

    black....................................................................Passed
    flake8...................................................................Passed

For more information, see `pre-commit Documentation <https://pre-commit.com/#quick-start>`_.

To commit your changes, run the following commands in your terminal:

.. code-block:: bash

    # after saving your code changes
    git status # check which files you have changed
    git diff fileThatYouHaveChanged # review changes
    git add fileThatYouHaveChanged
    git commit -m "Informative short message AlphaGenes/AlphaPeel#GitHubIssueNumber"

In the ``git add`` line above, don't use ``git add .``
because this last command will add all changes files to your commit,
including temporary files that might not belong in the repository.
Are you aware of `.gitignore file <https://git-scm.com/docs/gitignore>`_?

In the `git add` line above, don't use `git add .` because this last command will add all changes files to your commit, including temporary files that might not belong in the repository. Are you aware of <https://git-scm.com/docs/gitignore>_?

Update submodules?
------------------

Sometimes you have to update the submodules in line with
your code changes in the AlphaPeel or the submodules.

First, check the current state of the submodule:

.. code-block:: bash

    git submodule status

Next, check the latest commit in the submodule's remote repository:

.. code-block:: bash

    cd src/tinypeel/tinyhouse
    git log --oneline --max-count=1 origin/main
    cd ../../..

If the commit hashes match, then the submodule reference is up to date.
If you want to use the old submodule version, then a mismatch is ok.
Otherwise, update the reference using:

.. code-block:: bash

    git submodule update --remote
    git commit -m "Updated submodule reference to X.Y.Z AlphaGenes/AlphaPeel#GitHubIssueNumber"
    # provide submodules version (X.Y.Z) or commit hash

Create a pull request
---------------------

`Create a pull request (PR) <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request#creating-the-pull-request>`_ to propose your changes to the repository. Maintainers will review your PR.

Update the version of the package to publish the package
--------------------------------------------------------

.. note::

    This section is only for the maintainers to publish a new package version.

To release a new package version, we must update the ``version`` in ``pyproject.toml``.
For example, if the current version of the package is ``1.1.3`` and
the updated version should be ``1.1.4``, run:

.. code-block:: bash

    vi pyproject.toml

modify the following:

.. code-block:: toml

    ...
    [project]
    version = "1.1.3"
    ...

to

.. code-block:: toml

    ...
    [project]
    version = "1.1.4"
    ...


Remember to also update the version number in the test workflow file ``.github/workflows/tests.yml``, which is used to test the distribution built on the modified code:

.. code-block:: bash

    vi .github/workflows/tests.yml

modify the following:

.. code-block:: yaml

    - name: Install AlphaPeel
        run: pip install dist/alphapeel-1.1.3-py3-none-any.whl

to

.. code-block:: yaml

    - name: Install AlphaPeel
        run: pip install dist/alphapeel-1.1.4-py3-none-any.whl

commit the change:

.. code-block:: bash

    git commit -m "Bumped version to 1.1.4"

create the release with new version number according to `GitHub Docs <https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository#creating-a-release>`_.

The above will trigger workflow actions to publish the package on PyPI and documentation on Read the Docs:

  * `PyPI <https://pypi.org/project/AlphaPeel>`_
  * `Read the Docs <https://alphapeel.readthedocs.io/en/stable/index.html>`_

