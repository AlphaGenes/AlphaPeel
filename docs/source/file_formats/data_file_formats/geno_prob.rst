.. _geno_prob_file_format:

Genotype probability
====================

Format
~~~~~~

There are three lines per individual, corresponding to probabilities for
``aa``, ``aA`` or ``Aa``, and ``AA`` genotypes.
The first value in each line is the individual ID.
The remaining values are *unphased genotype probabilities* at each locus.

Example with four individuals and four loci:

::

  id1 0.7912 0.0000 0.0000 1.0000
  id1 0.2088 0.2179 0.5274 0.0000
  id1 0.0000 0.7820 0.4725 0.0000
  id2 0.0000 0.0000 0.0000 0.0001
  id2 1.0000 1.0000 1.0000 0.9999
  id2 0.0000 0.0000 0.0000 0.0000
  id3 0.3784 0.2171 0.0000 1.0000
  id3 0.4140 0.4329 0.0001 0.0000
  id3 0.2076 0.3500 0.9999 0.0000
  id4 0.9999 0.0000 0.0000 1.0000
  id4 0.0001 0.0001 0.9999 0.0000
  id4 0.0000 0.9999 0.0000 0.0000

When working with the X chromosome, for a female individual, the interpretation
is as for an autosomal chromosome above. For a male individual, the first line and 
the last line respectively correspond to ``a`` and ``A`` genotypes, while the middle
line is a placeholder with all probabilities set to 0. 

Example with four individuals and their X chromosome genotypes at four loci  
(id1 and id3 are males, while id2 and id4 are females):  


::

  id1 0.7912 0.2179 0.5274 1.0000
  id1 0.0000 0.0000 0.0000 0.0000
  id1 0.2088 0.7820 0.4725 0.0000
  id2 0.0000 0.0000 0.0000 0.0001
  id2 1.0000 1.0000 1.0000 0.9999
  id2 0.0000 0.0000 0.0000 0.0000
  id3 0.3784 0.2171 0.0001 1.0000
  id3 0.0000 0.0000 0.0000 0.0000
  id3 0.6216 0.7829 0.9999 0.0000
  id4 0.9999 0.0000 0.0000 1.0000
  id4 0.0001 0.0001 0.9999 0.0000
  id4 0.0000 0.9999 0.0000 0.0000

Output details
~~~~~~~~~~~~~~

The file is saved as ``.geno_prob.txt``. 
The file contains all individuals from all inputs.
