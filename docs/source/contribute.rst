Contribute
==========

To contribute to the package, please follow the procedures below.

Check submodule update
----------------------

Firstly, check if the submodule is up to date. This can be done by firstly checking the current state of the submodule:

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

We also need to update the ``version`` in ``pyproject.toml``. For example, if the current version of the package is ``1.1.3`` and the updated version should be ``1.1.4``, modify the following:

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

Publish via workflow (for administrators)
-----------------------------------------

.. note::

    This section is only for the repository administrators. Others should open a pull request instead.

Finally, we need to publish the updated package. In AlphaPeel, this is done by trigering the workflow by a tagged commit.

First, tag the commit with the corresponding version number. 

.. code-block:: console
    git tag 1.1.4

Then push the tagged commit.

.. code-block:: console
    git push --tags

This will trigger the workflow and publish a newer version of the project.