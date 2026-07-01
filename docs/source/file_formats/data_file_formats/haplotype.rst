.. _hap_file_format:

Haplotype
=========

Format
~~~~~~

There are two lines per individual.
The first line provides information on the paternal haplotype, and
the second line provides information on the maternal haplotype.
The first value in each line is the individual's ID.
The remaining values are alleles at each locus forming *haplotypes*
(:ref:`see the section on encoding alleles <zero_one_two_etc>`).

Example with four individuals and four loci:

::

  id1 0 1 1 0 # Paternal haplotype
  id1 0 1 1 0 # Maternal haplotype
  id2 1 0 1 0
  id2 0 1 0 1
  id3 0 1 1 0
  id3 1 0 1 0
  id4 0 1 1 0
  id4 0 1 0 0

When working with the X chromosome, for a female individual, the interpretation
is as for an autosomal chromosome above. For a male individual, 
the first line will have all alleles encoded with ``9`` 
(because males don't inherit X chromosome from the father), 
while the second line provides information on the maternal haplotype.  

Example with four individuals and their X chromosome genotypes at four loci
(id1 and id3 are males, while id2 and id4 are females):

::

  id1 9 9 9 9 # Paternal haplotype
  id1 0 1 1 0 # Maternal haplotype
  id2 1 0 1 0
  id2 0 1 0 1
  id3 9 9 9 9
  id3 1 0 1 0
  id4 0 1 1 0
  id4 0 1 0 0

Input details
~~~~~~~~~~~~~

The file contains *observed haplotypes* for each phased genotyped individual; 
possibly from previous analyses.
If provided, values from this file are translated to genotypes,
and the phase information is not preserved 
(see :ref:`genotype input details <genotype_input>`).
Therefore, the order of maternal and paternal haplotypes is not important.
The file does not need to include all individuals present in other files.
Only loci on one chromosome should be provided!

.. relate to issue: https://github.com/AlphaGenes/AlphaPeel/issues/279

The X chromosome functionality for haplotype input is not tested, so do not use it.

Output details
~~~~~~~~~~~~~~

The file is saved as ``.hap_THRESHOLD.txt``.
The file contains *called haplotypes* for all individuals from all inputs.
When probability is too low to make the call,
haplotypes are encoded as missing/unknown
(:ref:`see the note on encoding alleles and genotypes <zero_one_two_etc>`).