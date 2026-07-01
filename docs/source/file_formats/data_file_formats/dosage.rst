.. _dosage_file_format:

Dosage
======

Format
~~~~~~

There is one line per individual.
The first value in each line is the individual ID.
The remaining values are *expected allele dosages* at each locus.
These values will be between ``0`` and ``2``, inclusive
(:ref:`see the section on encoding alleles and genotypes <zero_one_two_etc>`).
These expected dosages are obtained by

(1) taking the :ref:`allele dosages <zero_one_two_etc>` (``0``, ``1``, and ``2``) of each individual at each locus,
(2) weighting them with corresponding :ref:`genotype probabilities <geno_prob_file_format>`, and
(3) summing the weighted dosages.

Example with four individuals and four loci:

::

  id1 0.2089 1.7820 1.4725 0.0000
  id2 1.0000 1.0000 1.0000 0.9999
  id3 0.8293 1.1329 1.9999 0.0000
  id4 0.0001 1.9999 1.0000 0.0000

Output details
~~~~~~~~~~~~~~

The file is saved as ``.dosage.txt``.
The file contains *expected dosages for the alternative allele*
for all individuals from all inputs.

