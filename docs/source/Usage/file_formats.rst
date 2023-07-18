============
File formats
============

Input file formats
------------------

Genotype file 
=============

Genotype files contain the input genotypes for each individual. The first value in each line is the individual's id. The remaining values are the genotypes of the individual at each locus, either 0, 1, or 2 (or 9 if missing). The following examples gives the genotypes for four individuals genotyped on four markers each.

Example: ::

  id1 0 2 9 0 
  id2 1 1 1 1 
  id3 2 0 2 0 
  id4 0 2 1 0

Sequence file
=============

The sequence data file is in a similar Sequence data is given in a similar format to the genotype data. For each individual there are two lines. The first line gives the individual's id and the read counts for the reference allele. The second line gives the individual's id and the read counts for the alternative allele.

Example: ::

  id1 4 0 0 7 # Reference allele for id1
  id1 0 3 0 0 # Alternative allele for id2
  id2 1 3 4 3
  id2 1 1 6 2
  id3 0 3 0 1
  id3 5 0 2 0
  id4 2 0 6 7
  id4 0 7 7 0

Pedigree file
=============

Each line of a pedigree file has three values, the individual's id, their father's id, and their mother's id. "0" represents an unknown id.

Example: ::

  id1 0 0
  id2 0 0
  id3 id1 id2
  id4 id1 id2

Binary plink file
=================

|Software| supports the use of binary plink files using the package ``AlphaPlinkPython``. |Software| will use the pedigree supplied by the ``.fam`` file if a pedigree file is not supplied. Otherwise the pedigree file will be used and the ``.fam`` file will be ignored. 


Map file 
========

The map file gives the chromosome number and the marker name and the base pair position for each marker in two columns. |Software| needs to be run with all of the markers on the same chromosome. 

Example: ::

  1 snp_a 12483939
  1 snp_b 192152913
  1 snp_c 65429279
  1 snp_d 107421759


Output file formats
-------------------

Phase file
==========

The phase file gives the phased haplotypes (either 0 or 1) for each individual in two lines. For individuals where we can determine the haplotype of origin, the first line will provide information on the paternal haplotype, and the second line will provide information on the maternal haplotype.

Example: ::

  id1 0 1 9 0 # Paternal haplotype
  id1 0 1 9 0 # Maternal haplotype
  id2 1 1 1 0
  id2 0 0 0 1
  id3 1 0 1 0
  id3 1 0 1 0 
  id4 0 1 0 0
  id4 0 1 1 0

Genotype probability file
=========================

The haplotype file (*.haps*) provides the (phased) allele probabilities for each locus. There are four lines per individual containing the allele probability for the (aa, aA, Aa, AA) alleles where the paternal allele is listed first, and where *a* is the reference (or major) allele and *A* is the alternative (or minor) allele. 

Example: ::

  id1    0.9998    0.0001    0.0001    1.0000
  id1    0.0000    0.4999    0.4999    0.0000
  id1    0.0000    0.4999    0.4999    0.0000
  id1    0.0001    0.0001    0.0001    0.0000
  id2    0.0000    1.0000    0.0000    1.0000
  id2    0.9601    0.0000    0.0455    0.0000
  id2    0.0399    0.0000    0.9545    0.0000
  id2    0.0000    0.0000    0.0000    0.0000
  id3    0.9998    0.0001    0.0001    1.0000
  id3    0.0000    0.4999    0.4999    0.0000
  id3    0.0000    0.4999    0.4999    0.0000
  id3    0.0001    0.0001    0.0001    0.0000
  id4    1.0000    1.0000    0.0000    1.0000
  id4    0.0000    0.0000    0.0000    0.0000
  id4    0.0000    0.0000    0.0000    0.0000
  id4    0.0000    0.0000    1.0000    0.0000

Dosage file
===========

The dosage file gives the expected allele dosage for the alternative (or minor) allele for each individual. The first value in each line is the individual ID. The remaining values are the allele dosages at each loci. These values will be between 0 and 2.

Example: ::

  1    0.0003    1.0000    1.0000    0.0001
  2    1.0000    0.0000    1.0000    0.0000
  3    0.0003    1.0000    1.0000    0.0001
  4    0.0000    0.0000    2.0000    0.0000


Segregation file
================

The segregation file gives the joint probability of each pattern of inheritance. There are four lines for each individual representing the probability of inheriting: 

  1. the grand **paternal** allele from the father and the grand **paternal** allele from the mother
  2. the grand **paternal** allele from the father and the grand **maternal** allele from the mother
  3. the grand **maternal** allele from the father and the grand **paternal** allele from the mother
  4. the grand **maternal** allele from the father and the grand **maternal** allele from the mother

Example: ::

  id1    1.0000    0.9288    0.9583    0.9834
  id1    0.0000    0.0149    0.0000    0.0000
  id1    0.0000    0.0554    0.0417    0.0166
  id1    0.0000    0.0009    0.0000    0.0000
  id2    0.9810    0.9842    1.0000    0.9971
  id2    0.0174    0.0158    0.0000    0.0013
  id2    0.0016    0.0000    0.0000    0.0016
  id2    0.0000    0.0000    0.0000    0.0000
  id3    0.0164    0.0149    0.0000    0.0065
  id3    0.9259    0.9288    0.9582    0.9769
  id3    0.0010    0.0009    0.0000    0.0001
  id3    0.0567    0.0554    0.0417    0.0165
  id4    0.0002    0.0000    0.0002    0.0004
  id4    0.0015    0.0000    0.0019    0.0041
  id4    0.1189    0.1179    0.1052    0.0834
  id4    0.8794    0.8821    0.8927    0.9122

Parameter files
===============

|Software| outputs three parameter files, ``.maf``, ``.seqError``, ``.genoError``. These give the minor allele frequency, sequencing error rates, and genotyping error rates used. All three files contain a single column with an entry for each marker. 

Example ``.maf`` file for four loci: 
::

  0.468005
  0.195520
  0.733061
  0.145847


.. |Software| replace:: AlphaPeel
