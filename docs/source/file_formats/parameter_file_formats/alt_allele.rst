.. _alt_allele_prob_file_format:

Alternative allele probability
==============================

Format
~~~~~~

There is one line per locus, with an optional header line at the top of the file.  
The header contains the metafounder ID(s),
which can be skipped when there is just one metafounder / no metafounders.
The remaining lines contain the metafounder-specific *alternative allele probabilities* for each locus.

Example with no metafounder and four loci:

::

  0.468005
  0.195520
  0.733061
  0.145847


Example with two metafounders and four loci:

::

  MF_1      MF_2
  0.468005  0.390000
  0.195520  0.219890
  0.733061  0.509849
  0.145847  0.090000


Input details
~~~~~~~~~~~~~

If provided,
values from this file initialises *alternative allele probabilities*.
If metafounders are defined,
the file should include all metafounders present in pedigrees.
Any missing values, or whole metafounder columns missing,
default to 0.5 for all loci.
If you do not have information for some loci, provide 0.5 for those loci.

Output details
~~~~~~~~~~~~~~

The file is saved as ``.alt_allele_prob.txt``.
