Contribution and Development Guide
==================================

Welcome to AlphaPeel contribution and development guide.

This guide will give you an overview of the contribution and development workflow via the AlphaPeel GitHub repository.

Fork the repository
-------------------

First, you should fork the `repository <https://github.com/AlphaGenes/AlphaPeel>`_.

For more information, see the `GitHub Docs <https://docs.github.com/en/get-started/quickstart/fork-a-repo#forking-a-repository>`_.

Clone your forked repository
----------------------------

Clone your forked repository into a local directory and initialise submodules at the same time by running the following command in your terminal:

.. code-block:: bash

    git clone --recurse-submodules URL_of_your_forked_repository

Depending on the type of code change, you should use different branches. First, check the available branches:

.. code-block:: bash

    cd AlphaPeel
    git branch # check available branches

Large code changes should go to dedicated development branches, which will be later merged into the `devel` branch by maintainers:

.. code-block:: bash

    git checkout devel # start from the development branch
    git branch fix_GitHubIssueNumber # fixing issue with GitHub number GitHubIssueNumber 
    git checkout fix_GitHubIssueNumber
    # now work on your code changes

If you are fixing an unknown issue or adding a new feature, open an issue first to document what you plan to do, then follow the above process.

Small code changes can go directly to the `devel` branch, which will eventually be merged into the main branch by maintainers:

.. code-block:: bash

    git checkout devel
    # now work on your code changes

Stable code for wider use that will be published is in the `main` branches. While most changes will be happening on the `devel` branch, critical bugfixes, can go to `main` branch:

.. code-block:: bash

    git checkout main
    # now work on your code changes

Make changes in your clone 
--------------------------

Make changes to the code and commit them to your local clone repository. Adding `#GitHubIssueNumber` in the message will link the commit with the issue page.

.. code-block:: bash

    # after saving your code changes
    git status # check which files you have changed
    git diff fileThatYouHaveChanged # review changes
    git add fileThatYouHaveChanged
    git commit -m "Informative short message #GitHubIssueNumber"

In the `git add` line above, don't use `git add .` because this last command will add all changes files to your commit, including temporary files that might not belong in the repository. Are you aware of <https://git-scm.com/docs/gitignore>_?

Update submodules?
------------------

Sometimes you have to update the submodules in line with your code changes in the AlphaPeel or the submodules.

First, check the current state of the submodule:

.. code-block:: bash

    git submodule status

Next, check the latest commit in the submodule's remote repository:

.. code-block:: bash
    
    cd src/tinypeel/tinyhouse
    git log --oneline --max-count=1 origin/main
    cd ../../..

If the commit hashes match, then the submodule reference is up to date. If you want to use the old submodule version, then a mismatch is ok. Otherwise, update the reference using:

.. code-block:: bash

    git submodule update --remote
    git commit -m "Updated submodule reference to X.Y.Z #GitHubIssueNumber"
    # provide submodules version (X.Y.Z) or commit hash

Create a pull request
---------------------

`Create a pull request (PR) <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request#creating-the-pull-request>`_ to propose your changes to the repository. Maintainers will review your PR.

Update the version of the package to publish the package
--------------------------------------------------------

.. note:: 

    This section is only for the repository maintainers to publish a new package version.

To release a new package version, we must update the ``version`` in ``pyproject.toml``. For example, if the current version of the package is ``1.1.3`` and the updated version should be ``1.1.4``, run:

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

commit the change:

.. code-block:: bash

    git commit -m "Bumped version to 1.1.4"

tag the version:

.. code-block:: bash

    git tag 1.1.4
    # git tag 1.1.4 --force # if you are reusing the tag

and push:

.. code-block:: bash

    git push # push code changes
    git push --tags # push tag changes
    # git push --tags --force # if you are reusing the tag

The above will trigger workflow actions to publish the package on PyPi and documentation on Read the Docs:

  * <https://pypi.org/project/AlphaPeel>_
  * <https://alphapeel.readthedocs.io/en/stable/index.html>_
