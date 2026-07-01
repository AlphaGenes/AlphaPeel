.. _phased_geno_prob_file_format:

Phased genotype probability
===========================

Format
~~~~~~

There are four lines per individual, corresponding to probabilities for
``aa``, ``aA``, ``Aa``, and ``AA`` phased genotypes
(:ref:`see the section on encoding alleles and genotypes for parental order of alleles <zero_one_two_etc>`).
The first value in each line is the individual ID.
The remaining values are *phased genotype probabilities* at each locus.

Example with four individuals and four loci:

::

  id1 0.7912 0.0000 0.0000 1.0000
  id1 0.1044 0.1090 0.2637 0.0000
  id1 0.1044 0.1090 0.2637 0.0000
  id1 0.0000 0.7820 0.4725 0.0000
  id2 0.0000 0.0000 0.0000 0.0001
  id2 0.3764 0.6611 0.0000 0.9628
  id2 0.6236 0.3388 1.0000 0.0371
  id2 0.0000 0.0000 0.0000 0.0000
  id3 0.3784 0.2171 0.0000 1.0000
  id3 0.4140 0.0000 0.0001 0.0000
  id3 0.0000 0.4328 0.0000 0.0000
  id3 0.2076 0.3500 0.9999 0.0000
  id4 0.9999 0.0000 0.0000 1.0000
  id4 0.0000 0.0000 0.2912 0.0000
  id4 0.0000 0.0000 0.7088 0.0000
  id4 0.0000 0.9999 0.0000 0.0000

When working with the X chromosome, for a female individual, the interpretation
is as for an autosomal chromosome above. For a male individual, the first line and 
the last line respectively correspond to ``a`` and ``A`` genotypes, while the middle 
two lines are placeholders with all probabilities set to 0.

Example with four individuals and their X chromosome genotypes at four loci
(id1 and id3 are males, while id2 and id4 are females):

::

  id1 0.7912 0.2179 0.5274 1.0000
  id1 0.0000 0.0000 0.0000 0.0000
  id1 0.0000 0.0000 0.0000 0.0000
  id1 0.2088 0.7820 0.4725 0.0000
  id2 0.0000 0.0000 0.0000 0.0001
  id2 0.3764 0.6611 0.0000 0.9628
  id2 0.6236 0.3388 1.0000 0.0371
  id2 0.0000 0.0000 0.0000 0.0000
  id3 0.3784 0.2171 0.0001 1.0000
  id3 0.0000 0.0000 0.0000 0.0000
  id3 0.0000 0.0000 0.0000 0.0000
  id3 0.6216 0.7829 0.9999 0.0000
  id4 0.9999 0.0000 0.0000 1.0000
  id4 0.0000 0.0000 0.2912 0.0000
  id4 0.0000 0.0000 0.7088 0.0000
  id4 0.0000 0.9999 0.0000 0.0000


Input details
~~~~~~~~~~~~~

The file contains
*phased genotype probabilities* for each genotyped and phased individual;
possibly from previous analyses.
The file does not need to include all individuals present in other files.
These phased genotype probabilities are assumed without an error and
used directly as individual's phased genotype probability state. 
For more information about the internal state, 
check :ref:`details in peeling basics <peeling_basics>`.
Only loci on one chromosome should be provided!

Output details
~~~~~~~~~~~~~~

The file is saved as ``.phased_geno_prob.txt``.
The file contains all individuals from all inputs.
