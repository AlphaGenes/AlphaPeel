.. _pheno_prob_file_format:

Phenotype probability
=====================

Format
~~~~~~

There are ``1+k`` lines per individual
(an empty line and ``k`` lines with probabilities for ``0:k`` phenotypes).
The first value in each line is the individual ID.
The remaining values are *phenotype probabilities*.
Current phenotype functionality works only with one locus and one trait!

Example with four individuals and two phenotype states for a binary trait:

::

  id1 0.9000
  id1 0.1000

  id2 0.9000
  id2 0.1000

  id3 0.2453
  id3 0.7547

  id4 0.9000
  id4 0.1000

Output details
~~~~~~~~~~~~~~

The file is saved as ``.pheno_prob.txt``. 
The file contains all individuals from all inputs.

.. TODO: remove the empty lines in the output files? (1+k will be k in text!) 
.. https://github.com/AlphaGenes/AlphaPeel/issues/247

