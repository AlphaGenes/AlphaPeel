.. _pheno_file_format:

Phenotype
=========

Format
~~~~~~

The first value in each line is the individual's ID.
The remaining values are the *phenotypes* of the individual for a specific trait,
coded from ``0`` onwards.
For example, a binary trait should be coded as ``0`` and ``1``,
while a trait with three states should be coded as ``0``, ``1``, and ``2``.
Current phenotype functionality works only with one locus and one trait!

Example with four individuals and their phenotypes for a binary trait:

::

  id1 0
  id2 1
  id3 1
  id4 0

Input details
~~~~~~~~~~~~~

This file has one or multiple lines with
*observed phenotypes* for each phenotyped individual.
Multiple lines per individual support repeated phenotyping.
The file does not need to include all individuals present in other files.
Non-phenotyped individuals should not be present in the file; 
that is, missing phenotypes are not to be included in the file.

If provided, values from this file are translated to :ref:`penetrance <penetrance>` 
component of genotype probabilities.
