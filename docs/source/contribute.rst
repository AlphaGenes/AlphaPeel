Contribution and Development Guide
==================================

Welcome to ``AlphaPeel`` contribution and development guide.

This guide will give you an overview of the contribution and development workflow
via the ``AlphaPeel`` GitHub repository at `<https://github.com/AlphaGenes/AlphaPeel>`_.

Critically, see also a list of issues at `<https://github.com/AlphaGenes/AlphaPeel/issues>`_.


Introducton to collaborative development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before diving into the technical aspects of contribution, 
it's essential to understand our collaborative development workflow. 
This section will outline the rules and best practices for creating issues, branches, and pull requests, 
as well as our overall workflow from defining tasks to reviewing code.

Issue
-----

To efficiently manage our development process, 
we use tags to categorize issues.
Following are two main issue tags that we use to label issues:

1. Discussion/Idea

    * Label: ``discussion``

    * No strict format

    * Can stay messy

2. Actionable Task

    * Label: ``task``
    
    * Must follow a template

**Rule**: Only task issues can have branches/PRs


.. _task-issue-template:

Task issue template
^^^^^^^^^^^^^^^^^^^

This template is available when you create a new issue.

.. code-block:: markdown

    ## Goal
    What are we trying to achieve?

    ## Proposed approach
    (brief, not perfect)

    ## Scope
    What files/modules are affected?

    ## Done when
    Clear condition for completion

Branch
------

We use branches to work on issues.

Rules
^^^^^

* **Branches** are for development, **Version tags** are for releases

    * ``devel`` → integration branch

    * ``main`` → release branch

    * version tags only on ``main``

* **Flow**:

    * Feature → PR into ``devel``

    * Stabilise ``devel``

    * Merge ``devel`` → ``main``

    * Version tag on ``main``: ``v1.2.0``, ``latest``, ``stable``


.. _branch-naming:

Branch naming
^^^^^^^^^^^^^

A suggested branch naming convention for a task issue on your fork is as follows:

.. code-block:: 

    <type>/issue-<id>-short-description

* Examples:
    * ``feat/issue-42-add-xchrom``
    * ``fix/issue-17-type-error``
    * ``refactor/issue-33-clean-io``
    * ``doc/issue-223-add-collab-guide``

* Types:
    * ``feat``: new functionality
    * ``fix``: bug fix
    * ``refactor``: code cleanup
    * ``doc``: documentation updates
    * etc.

.. _PR-practice:

Pull request
------------

A pull request (PR) is a request to merge your code changes from your branch into the main repository's branch (e.g., ``devel``).

The followings are the templates for PRs of ``AlphaPeel`` and ``tinyhouse``, 
which the latter is the submodule referenced by the former. 
The PR template should be followed when you create a PR to either 
the ``AlphaPeel`` repository or the ``tinyhouse`` repository.

We encourage you to open PRs early, 
even when the code is not fully ready, 
to facilitate early feedback and discussion. 
You can use draft PRs to indicate that the PR is a work in progress.

Template for PR of ``AlphaPeel``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: markdown

    ## Related Issue
    Closes AlphaGenes/AlphaPeel#GitHubIssueNumber

    ## What changed
    - Brief summary of changes
    - Key files/modules affected

    ## Why this change
    - What problem does this solve?
    - What functionality does it enable?

    ## Submodule changes
    - Submodule repository: <name>
    - PR: <link> (if applicable)
    - Commit: <hash>
    - Why needed: <what functionality depends on this>

    ## Notes / Risks
    - Anything reviewers should pay attention to


Template for PR of ``tinyhouse``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: markdown

    ## Related Issue
    <submodule-issue-id or main-repository issue link>

    ## What changed
    - Summary of changes

    ## Why this change
    - What functionality this enables

    ## Used by
    - Main repository PR: <link> (if applicable)

    ## Notes
    - Any compatibility considerations

Full workflow
-------------

Below is the full workflow for contribution and development, 
from defining the work to merging a pull request.

1. Create and define the work
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    In ``AlphaPeel`` repository:

    * Create a **task issue**
    * Clearly define:

        * Goal
        * Scope
        * “Done when”

    If ``tinyhouse`` changes are needed:

    * Create a **linked issue in** ``tinyhouse``


2. Create branches (both repos)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    In your fork:

    * ``AlphaPeel`` repository example:

        .. code-block:: bash

            feat/issue-42-new-analysis
        

    * ``tinyhouse`` repository example:

        .. code-block:: bash

            feat/issue-15-support-analysis
    
    For instructions, see :ref:`branch-instructions`.


3. Implement ``tinyhouse`` changes first
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    Steps:

    * Implement changes in ``tinyhouse`` branch
    * Push to your fork
    * Open **PR to** ``devel`` **branch of** ``tinyhouse`` **repository**


4. Use updated ``tinyhouse`` reference in ``AlphaPeel`` repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    In your ``AlphaPeel`` repository branch:

        * Point ``tinyhouse`` submodule to your **branch commit**
        * Continue development and testing the changes in ``AlphaPeel`` repository


5. Open PRs early
^^^^^^^^^^^^^^^^^

    Open:

    * ``tinyhouse`` PR → target: ``AlphaGenes/tinyhouse:devel``
    * ``AlphaPeel`` PR → target: ``AlphaGenes/AlphaPeel:devel``

    Link them both ways:

    * ``AlphaPeel`` PR → links ``tinyhouse`` PR
    * ``tinyhouse`` PR → links ``AlphaPeel`` PR


6. Development + syncing
^^^^^^^^^^^^^^^^^^^^^^^^

    During development, always keep your branch synced with upstream. 
    
        * ``rebase`` to the latest ``devel`` branch of the AlphaGenes repository to keep up with the latest changes and avoid merge conflicts later on.

            * Instructions: :ref:`rebase_instructions`

7. Review phase
^^^^^^^^^^^^^^^

    * First: ``tinyhouse`` PR

        * Review and merge ``tinyhouse`` PR

        * Now you have a **stable commit hash**

    * Then: Update ``AlphaPeel`` repository
    
        * Update ``tinyhouse`` submodule pointer to the **merged commit**

        * Push update to ``AlphaPeel`` repository branch

    * Then: ``AlphaPeel`` repository PR review

8. Merge flow
^^^^^^^^^^^^^

    * Before merge: 

        * ``squash`` the commits to simplify the history if applicable

            * This is optional

            * Be careful when squashing the merged, rebased, or squashed commits, as it can cause issues with the commit history and may cause GitHub to not able to merge the PR

            * Instructions: :ref:`squash-commits-instructions`

        * Final check the new changes by others and update with ``rebase`` if needed

            * Instructions: :ref:`rebase_instructions`


Developer meetings
------------------

We have routine developer meetings to discuss the ongoing work and plan for the next steps.

In the meeting:

* Each person shares:

    * What issue they’re working on

    * Any blockers

    * Any upcoming PRs

* Draft future plans:

    * New tasks to define


Technical Contribution Guide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fork the repository
-------------------

First, you should fork the `AlphaGenes repository <https://github.com/AlphaGenes/AlphaPeel>`_.

For more information, see the `GitHub Docs <https://docs.github.com/en/get-started/quickstart/fork-a-repo#forking-a-repository>`_.

Clone your forked repository
----------------------------

Clone your forked repository into a local directory and
initialise submodules at the same time by running the following command in your terminal:

.. code-block:: bash

    git clone --recurse-submodules URL_of_your_forked_repository

.. _use-correct-version-submodule:

Use correct version of the submodule
------------------------------------

In the above command, the ``--recurse-submodules`` flag will automatically initialise and update the submodule ``tinyhouse`` 
in the repository according to the ``.gitmodules`` file, it is expected to be the commit that is stable for the current version of ``AlphaPeel``,
but it might not be the latest version of the submodule.

Usually the most up-to-date version of the code of the submodule ``tinyhouse`` is in the ``devel`` branch of its repository, 
so if you need to use the lastest version, you can run the following command to checkout the ``devel`` branch of the submodules:  

.. code-block:: bash

    cd AlphaPeel
    cd src/tinypeel/tinyhouse
    git fetch origin
    git checkout devel # checkout the devel branch of the submodule
    cd ../../..

Here we are using the ``devel`` branch of the submodule as an example,
but you can checkout to any branch or commit of the submodule that 
you want to use for your development in the ``AlphaPeel`` repository.
If you are working with a specific branch of commit of your fork of the submodule, 
you can checkout to that branch or commit instead.

.. code-block:: bash

    cd AlphaPeel
    cd src/tinypeel/tinyhouse
    git remote add fork URL_of_your_forked_submodule_repository
    git fetch fork
    git checkout commit_or_branch_of_your_forked_submodule 
    cd ../../..

.. _branch-instructions:

Work with branches
------------------

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
    git branch <your-branch> # create a new branch for your issue
    git checkout <your-branch>
    # now work on your code changes

For branch naming convention, see :ref:`branch-naming`.

It is a good practice that you first open a task issue to document what you plan to do,
then follow the above process. For the task issue template, see :ref:`task-issue-template`.

Small code changes can go directly to the ``devel`` branch,
which will eventually be merged into the main branch by maintainers,
but check this with the maintainers when you open/discuss the issue:

.. code-block:: bash

    git checkout devel
    # now work on your code changes

.. _rebase_instructions:

Rebase to keep up with the latest changes
-----------------------------------------

During development, always keep your branch synced with upstream. 
You can do this by rebasing to the latest ``devel`` branch of the AlphaGenes repository to 
keep up with the latest changes and avoid merge conflicts later on.

To rebase to the AlphaGenes repository, you need to first add the AlphaGenes repository as an upstream remote:

.. code-block:: bash

    git remote add AlphaGenes https://github.com/AlphaGenes/AlphaPeel
    git fetch AlphaGenes

Then, you can rebase your branch to the latest ``devel`` branch of the AlphaGenes repository:

.. code-block:: bash

    git checkout <your-branch>
    git rebase AlphaGenes/devel

.. _changes_instructions:

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

To build and reinstall the package on your modified code:

.. code-block:: bash

    python -m build
    python -m pip uninstall AlphaPeel -y
    python -m pip install dist/*.whl

To run ``pytest``:

.. code-block:: bash

    pytest

If you want to run a specific functional test, such as the ``test_files``, you can run like the following:

.. code-block:: bash

    pytest tests/functional_tests/run_func_test.py::TestClass::test_files

If some functioanl tests fail and you want to see the output of the tests, you can add the ``-s`` flag, note that the accuracy test report cannot be generated with the ``-s`` flag:

.. code-block:: bash

    pytest -s tests/functional_tests/run_func_test.py

If the tests run successfully, you are expected to see the head of the output similar to the following:

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

In the ``git add`` line above, don't use ``git add .`` because this last command will add all changes files to your commit, 
including temporary files that might not belong in the repository. Are you aware of `<https://git-scm.com/docs/gitignore>`_?

.. _documentation-changes:

Documentation changes
---------------------

If you make documenation changes, you can build the documentation locally to check if the changes are correctly reflected in the documentation.

.. code-block:: bash

    cd docs
    make html

Then, you can open the generated HTML files in the ``build/html`` directory.

.. _submodule-changes:

Submodule consistency
---------------------

Before you push your changes to your fork, make sure the submodule pointer is updated to the commit 
that you want to use for your development in the ``AlphaPeel`` repository.

It is recommended to check if your code changes in the ``AlphaPeel`` is consistent with the latest version of the submodule 
``tinyhouse`` in the ``devel`` branch of its repository.

You can find the instruction to update the submodule pointer in the section :ref:`use-correct-version-submodule` above.

After update the submodule pointer, make sure to test your code changes by running the tests on 
the distribution built on your modified code to check if the code changes are consistent with the latest version of the submodule.

Before open a pull request
--------------------------

Before you open a pull request, make sure you have:

    * Tested your code changes by running the tests on the distribution built on your modified code

        * Instructions: :ref:`changes_instructions`

    * Committed your changes with informative commit messages

        * Instructions: :ref:`changes_instructions`

    * Updated the submodule pointer to the commit that you want to use for your development in the ``AlphaPeel`` repository

        * Instructions: :ref:`submodule-changes`

    * If you have made documentation changes, built the documentation locally to check if the changes are correctly reflected in the documentation 

        * Instructions: :ref:`documentation-changes`

    * You can add description of your change in the ``docs/changelog.rst`` file to keep track of the changes you have made, but this is optional and can be done later by maintainers when merging the PR

Also, make sure you have the lastest code changes from the AlphaGenes repository by rebasing to the latest ``devel`` branch of the AlphaGenes repository (see :ref:`rebase_instructions`).

.. _squash-commits-instructions:

Squash commits
^^^^^^^^^^^^^^

An optional method for cleaning up the commit history is to squash the commits, 
but be careful when squashing the merged, rebased, or squashed commits, 
as it can cause issues with the commit history and may cause GitHub to not be able to merge the PR.
To squash the last N commits, you can run the following command:

.. code-block:: bash

    git reset --soft HEAD~N
    git commit -m "Informative short message AlphaGenes/AlphaPeel#GitHubIssueNumber"

If you have already pushed the commits to your fork, you can use the following command to force push the squashed commit:

.. code-block:: bash

    git push --force-with-lease origin <your-branch>


Create a pull request
---------------------

`Create a pull request (PR) <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request#creating-the-pull-request>`_ to propose your changes to the repository. Maintainers will review your PR.

Please follow the PR best practices outline in the section :ref:`PR-practice` above when you create a PR.


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

