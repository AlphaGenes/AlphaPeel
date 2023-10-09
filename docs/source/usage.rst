-----
Usage
-----

===============
Program options
===============

|Software| takes in several command line arguments to control the program's behaviour. To view a list of arguments, run |Software| without any command line arguments, i.e. ``AlphaPeel`` or ``AlphaPeel -h``. 

Input Arguments 
---------------

::

    Input Options:
      -pedigree [PEDIGREE [PEDIGREE ...]]
                          Pedigree file(s) (see format below).
      -genotypes [GENOTYPES [GENOTYPES ...]]
                          Genotype File(s) (see format below).
      -seqfile [SEQFILE [SEQFILE ...]]
                          Sequence allele read count file(s) (see format below).
      -bfile [BFILE [BFILE ...]]
                          Plink (binary) file(s).
      -startsnp STARTSNP  The first marker to consider. The first marker is "1".
      -stopsnp STOPSNP    The last marker to consider.

|Software| requires a pedigree file (``-pedigree``) and one or more genomic data files to run the analysis.

|Software| supports the following genomic data files: genotype files in the AlphaGenes format (``-genotypes``), sequence allele read in the AlphaGenes format (``-seqfile``), and binary Plink files (``-bfile``). Use of binary Plink files requires the package ``alphaplinkpython``, which  can be installed via ``pip``, but is only stable for Linux. There are known issues with this package, so we do not advocate its use at the moment.

Use the ``-startsnp`` and ``-stopsnp`` to run the analysis only on a subset of markers.

The input options in the form of ``[xxx [xxx ...]]`` can take in more than one input file seperated by space.

Output Arguments 
----------------

::

    Output options:
      -out prefix           The output file prefix. All file outputs will be stored
                            as "prefix.dosage" and so on.
      -writekey WRITEKEY    Determines the order in which individuals are ordered
                            in the output file based on their order in the
                            corresponding input file. Individuals not in the input
                            file are placed at the end of the file and sorted in
                            alphanumeric order. These individuals can be suppressed
                            with the "-onlykeyed" option. Options: id, pedigree,
                            genotypes, sequence, segregation. Default: id.
      -onlykeyed            Flag to suppress the individuals not present in
                            the file used with "-outputkey". It also suppresses "dummy"
                            individuals.
      -iothreads IOTHREADS  Number of threads to use for input/output. Default: 1.


    Peeling output options:
      -no_dosages           Flag to suppress the dosage files.
      -no_seg               Flag to suppress the segregation files (e.g. when
                            running for chip imputation and not hybrid peeling).
      -no_params            Flag to suppress writing the parameter files.
      -haps                 Flag to enable writing out the genotype probabilities.
      -calling_threshold [CALLING_THRESHOLD [CALLING_THRESHOLD ...]]
                            Genotype calling threshold(s). Multiple space
                            separated values allowed. Use .3 for best guess
                            genotype.
      -binary_call_files    Flag to write out the called genotype files as a
                            binary plink output [Not yet implemented].

By default |Software| produces a dosages file, a segregation files and two parameter files (genotyping error and recombination rate). Creation of these files can be suppressed with the ``-no_dosages``, ``-no_seg``, and ``-no_params`` options. |Software| can also write out the genotype probability file (.haps) with the `-haps` argument.

The ``-calling_threshold`` arguments controls which genotypes (and phased haplotypes) are called as part of the algorithm. A calling threshold of 0.9 indicates that genotypes are only called if greater than 90% of the final probability mass is on that genotype. Using a higher-value will increase the accuracy of called genotypes, but will result in fewer genotypes being called. Since there are three genotypes states,  "best-guess" genotypes are produced with a calling threshold less than ``0.33``. ``-calling_threshThe ``-binary_call_files`` option can be used to change the output to a plink binary format. 

The order in which individuals are output can be changed by using the ``writekey`` option. This option changes the order in which individuals are written out to the order in which they were observed in the corresponding file. The ```-onlykeyed`` option suppresses the output of dummy individuals (not recommended for hybrid peeling). 

The argument ``-iothreads`` controls the number of threads/processes used by |Software|. |Software| uses additional threads to parse and format input and output files. Setting this option to a value greater than 1 is only recommended for very large files (i.e. >10,000 individuals).

Peeling arguments 
------------------

::

    Mandatory peeling arguments:
      -runtype RUNTYPE      Program run type. Either "single" or "multi".
    
    Optional peeling arguments:
      -ncycles NCYCLES      Number of peeling cycles. Default: 5.
      -maxthreads MAXTHREADS
                            Number of threads to use. Default: 1.
      -length LENGTH        Estimated length of the chromosome in Morgans.
                            [Default 1.00]

    Peeling control arguments:
      -esterrors            Flag to re-estimate the genotyping error rates after
                            each peeling cycle.
      -estmaf               Flag to re-estimate the minor allele frequency after
                            each peeling cycle.
      -nophasefounders      A flag phase a heterozygous allele in one of the
                            founders (if such an allele can be found).
      -sexchrom             A flag to that this is a sex chromosome. Sex needs to
                            be given in the pedigree file. This is currently an
                            experimental option.

    Genotype probability arguments:
      -error ERROR          Genotyping error rate. [Default 0.01]
      -seqerror SEQERROR    Assumed sequencing error rate. [Default 0.001]

``-runtype`` controls whether the program is run in "single-locus" or "multi-locus" model. Single locus mode does not use linkage information to perform imputation. It is fast, but not very accurate. Multi-locus mode runs multi-locus iterative peeling which uses linkage information to increase accuracy and calculate segregation values.

For hybrid peeling, where a large amount (millions of segregating sites) of sequence allele read counts needs to be imputed, first run the program in multi-locus mode to generate a segregation file, and then run the program in single-locus mode with a known segregation file.

The ``-error``, ``-seqerror`` and ``-length`` arguments control some of the parameters used in the model. ``-seqerror`` must not be zero. |Software| is robust to deviations in genotyping error rate and sequencing error rate so it is not recommended to use these options unless large deviations from the default are known. Changing the ``-length`` argument to match the genetic map length can increase accuracy in some situations.

The ``-esterrors`` option estimated the genotyping error rate based on observed information, this option is generally not necessary and can increase runtime. ``-estmaf`` estimates the minor allele frequency after each peeling cycle. This option can be useful if there are a large number of non-genotyped founders. 

Hybrid peeling arguments 
------------------------

::

    Single locus arguments:
      -mapfile MAPFILE      A map file (chr marker_name position) for genotype data.
      -segmapfile SEGMAPFILE
                            a map file for the segregation estimates for hybrid
                            peeling.
      -segfile SEGFILE      A segregation file for hybrid peeling.

In order to run hybrid peeling the user needs to supply a ``-mapfile`` which gives the genetic positions for the SNPs in the sequence allele read counts data supplied, a ``-segmapfile`` which gives the genetic position for the SNPs in the segregation file, and a ``-segfile`` which gives the segregation values generated via multi-locus iterative peeling. These arguments are not required for running in multi-locus mode.

============
File formats
============

Input file formats
------------------

Pedigree file
=============

Each line of a pedigree file has three values, the individual's id, their father's id, and their mother's id. "0" represents an unknown id.

Example:

::

  id1 0 0
  id2 0 0
  id3 id1 id2
  id4 id1 id2

Genotype file 
=============

Genotype files contain the input genotypes for each individual. The first value in each line is the individual's id. The remaining values are the genotypes of the individual at each locus, either 0, 1, or 2 (or 9 if missing). The following examples gives the genotypes for four individuals genotyped on four markers each.

Example:

::

  id1 0 2 9 0 
  id2 1 1 1 1 
  id3 2 0 2 0 
  id4 0 2 1 0

Sequence allele read counts file
================================

The sequence allele read counts file has two lines for each individual. The first line gives the individual's id and read counts for the reference allele. The second line gives the individual's id and allele read counts for the alternative allele.

Example:

::

  id1 4 0 0 7 # Reference allele for id1
  id1 0 3 0 0 # Alternative allele for id1
  id2 1 3 4 3
  id2 1 1 6 2
  id3 0 3 0 1
  id3 5 0 2 0
  id4 2 0 6 7
  id4 0 7 7 0

Binary plink file
=================

Binary Plink files are supported using the package ``AlphaPlinkPython``. The pedigree supplied by the ``.fam`` file will be used if a pedigree file is not supplied. Otherwise, the pedigree file will be used and the ``.fam`` file will be ignored. 

Map file 
========

The map file gives the chromosome number, the marker name, and the base pair position for each marker in two columns. Only markers on one chromosome should be provided! 

Example:

::

  1 snp_a 12483939
  1 snp_b 192152913
  1 snp_c 65429279
  1 snp_d 107421759


Output file formats
-------------------

Phase file
==========

The phase file gives the phased haplotypes (either 0 or 1) for each individual in two lines. For individuals where we can determine the haplotype of origin, the first line will provide information on the paternal haplotype, and the second line will provide information on the maternal haplotype.

Example:

::

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

Example:

::

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

Example:

::

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

Example:

::

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
