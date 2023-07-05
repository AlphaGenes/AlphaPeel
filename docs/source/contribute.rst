Contributing Guide
==================

Welcome to AlphaPeel contributing guide.

In this guide you will get an overview of the contribution workflow via the AlphaPeel GitHub repository.

Fork the repository
-------------------

First you should fork the `repository <https://github.com/AlphaGenes/AlphaPeel>`_.

For more information, see the `GitHub Docs <https://docs.github.com/en/get-started/quickstart/fork-a-repo#forking-a-repository>_`.

Clone your forked repository
----------------------------

Clone your forked repository into a local directory and initialise submodules at the same time by running the following command in your terminal:

.. code-block:: console

    git clone --recurse-submodules https://github.com/AlphaGenes/AlphaPeel.git

Depending on the type of change you should use either main, devel, or other branch

.. code-block:: console

    cd AlphaPeel
    git branch # check branches
    # git checkout main # 
    git checkout devel # most changes should go here and later merged to main by maintainers
    git branch fix_GitHubIssueNumber # fixing issue with GitHub number GitHubIssueNumber 
    git checkout fix_GitHubIssueNumber

Make changes in your clone 
--------------------------

Make changes and commit them to your local clone repository. Adding #GitHubIssueNumber in the message links with the issue page.

.. code-block:: console

    git commit -m "Informative short message #GitHubIssueNumber"

Update submodules?
------------------

Sometimes you have to update the submodules in line with your code changes in the AlphaPeel or in the submodules.

First check the current state of the submodule:

.. code-block:: console

    git submodule status

Next, check the lateset commit in the submodule's remote repository:

.. code-block:: console
    
    cd src/tinypeel/tinyhouse
    git log --oneline --max-count=1 origin/main
    cd ../../..

If the commit hashes match, then the submodule reference is up to date. If you want to use the old submodule version, then missmatch is ok. Otherwise, update the referece using:

.. code-block:: console

    git submodule update --remote
    git commit -m "Updated submodule reference to X.Y.Z #GitHubIssueNumber" # provide version or hash

Create a pull request
-----------------------

`Create a pull request (PR) <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request#creating-the-pull-request>` to propose your changes to the repository. Your PR will be reviewed by maintainers.

Update the version of the package to publish the package
--------------------------------------------------------

.. note:: 

    This section is only for the repository maintainers to publish a new version of the package.

To release a new version of the package, we must update the ``version`` in ``pyproject.toml``. For example, if the current version of the package is ``1.1.3`` and the updated version should be ``1.1.4``, run:

.. code-block:: console

    vi pyproject.toml

modify the following:

.. code-block:: toml

    ...
    [project]
    name = "AlphaPeel"
    version = "1.1.3"
    ...

to 

.. code-block:: toml

    ...
    [project]
    name = "AlphaPeel"
    version = "1.1.4"
    ...

commit the change:

.. code-block:: console

    git commit -m "Bumped version to 1.1.4"

tag the version:

.. code-block:: console

    git tag 1.1.4
    # git tag 1.1.4 --force # if you are reusing the tag

and push:

.. code-block:: console

    git push # do we need this one or just the next one?
    git push --tags
    # git push --tags --force # if you are reusing the tag

The above will trigger workflow actions to publish the package on PyPi and documentation on Read the Docs:

  * <https://pypi.org/project/AlphaPeel>_
  * <https://alphapeel.readthedocs.io/en/stable/index.html>_
