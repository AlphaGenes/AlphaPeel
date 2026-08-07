.. _seq_file_format:

Sequence allele read counts
===========================

Format
~~~~~~

There are two lines per individual.
The first line gives the individual's ID and *allele read counts* for the :ref:`reference allele <zero_one_two_etc>` at each locus.
The second line gives the individual's ID and *allele read counts* for the :ref:`alternative allele <zero_one_two_etc>` at each locus.

.. TODO: work with Yuni on X seq data input
.. https://github.com/AlphaGenes/AlphaPeel/issues/248

The same format works for the autosomes and the X chromosome.
Only loci on one chromosome should be provided!

Example with four individuals and their allele read counts at four loci:

::

  id1 4 0 0 7 # Reference allele read counts for id1
  id1 0 3 0 0 # Alternative allele read counts for id1
  id2 1 3 4 3
  id2 1 1 6 2
  id3 0 3 0 1
  id3 5 0 2 0
  id4 2 0 6 7
  id4 0 7 7 0

Input details
~~~~~~~~~~~~~

If provided, values from this file are translated to the 
:ref:`penetrance <penetrance>` component of genotype probabilities.
The file does not need to include all individuals present in other files.
Only loci on one chromosome should be provided!
