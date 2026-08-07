.. _geno_file_format:

Genotype
========

Format
~~~~~~

There is one line per individual.
The first value in each line is the individual's ID.
The remaining values are *genotypes* encoded as *allele dosages* at each locus
(:ref:`see the section on encoding alleles and genotypes <zero_one_two_etc>`).  

Example with four individuals and their genotypes at four loci:

::

  id1 0 2 9 0
  id2 1 1 1 1
  id3 2 0 2 0
  id4 0 2 1 0

Example with four individuals and their X chromosome genotypes at four loci:

id1 and id3 are males, while id2 and id4 are females:

::

  id1 0 1 9 0
  id2 1 1 1 1
  id3 1 0 1 0
  id4 0 2 1 0

.. _genotype_input:

Input details
~~~~~~~~~~~~~

This file has one line with
*observed genotypes* for each genotyped individual.
If provided, the values from this file initialises 
the :ref:`penetrance <penetrance>` component of genotype probabilities.
The file does not need to include all individuals present in other files.
Only loci on one chromosome should be provided!

Output details
~~~~~~~~~~~~~~

The file is saved as ``.geno_THRESHOLD.txt``.
The file contains all individuals from all inputs.
When probability is too low to make the call,
genotypes are encoded as missing/unknown
(:ref:`see the note on encoding alleles and genotypes <zero_one_two_etc>`).
