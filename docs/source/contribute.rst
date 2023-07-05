Contributing Guide
==================

Welcome to AlphaPeel contributing guide.

In this guide you will get an overview of the contribution workflow.

Forking the repository
----------------------

Fork the `repository <https://github.com/AlphaGenes/AlphaPeel>`_.

* For more information, see the `GitHub Docs <https://docs.github.com/en/get-started/quickstart/fork-a-repo?tool=webui&platform=mac#forking-a-repository>`_.

Cloning your forked repository
------------------------------

Clone your forked repository into a local directory, initialized the submodule at the same time by copying the following command to the terminal:

.. code-block:: console

    git clone --recurse-submodules https://github.com/AlphaGenes/AlphaPeel.git


Making changes to your fork
---------------------------

Commit the changes you made.

Check submodule update
----------------------

Before open a pull request, check if the submodule is up to date. This can be done by firstly checking the current state of the submodule:

.. code-block:: console

    git submodule status

Next, verify the lateset commit in the submodule's remote repository:

.. code-block:: console
    
    cd src/tinypeel/tinyhouse
    git log --oneline --max-count=1 origin/main

If the commit hashes match, then the submodule reference is up to date. Otherwise, please update the referece using the following command:

.. code-block:: console

    git submodule update --remote

Remember to commit the changes.

Update the version of the package
---------------------------------

.. note:: 

    This section is only for the repository maintainers when a new version of the package is required to release, otherwise you can skip it.

To release a new version of the package, we also need to update the ``version`` in ``pyproject.toml``. For example, if the current version of the package is ``1.1.3`` and the updated version should be ``1.1.4``, modify the following:

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

Remember to commit the changes.

Creating a pull request
-----------------------

`Create a pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request#creating-the-pull-request>`_ to propose your changes to the repository. Your PR will be reviewed by some of the maintainers.

Publish via actions
-------------------

.. note::

    This section is only for the repository administrators. Others should stop here.

Finally, we need to publish the updated package. In AlphaPeel, this is done by trigering the workflow by a tagged commit.

First, tag the commit with the corresponding version number: 

.. code-block:: console

    git tag 1.1.4

Then push the tagged commit:

.. code-block:: console

    git push --tags

This will trigger the workflow and publish a new version of the project.