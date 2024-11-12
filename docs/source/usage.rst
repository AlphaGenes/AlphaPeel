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
      -ped_file [PEDIGREE ...]
                          Pedigree file(s) (see format below).
      -geno_file [GENOTYPES ...]
                          Genotype file(s) (see format below).
      -seq_file [SEQFILE ...]
                          Sequence allele read count file(s) (see format below).
      -plink_file [BFILE ...]
                          Plink (binary) file(s).
      -alt_allele_prob_file [ALT_ALLELE_PROB_FILE...]
                          The alternative allele probabilities per metafounder(s). Default: 0.5 per marker.
      -start_snp START_SNP
                          The first marker to consider. The first marker is "1". Default: 1.
      -stop_snp STOP_SNP
                          The last marker to consider. Default: all markers considered.
      -main_metafounder MAIN_METAFOUNDER
                          The metafounder to use where parents are unknown with input "0". Default: MF_1.

|Software| requires a pedigree file (``-ped_file``) and one or more genomic data files to run the analysis.

|Software| supports the following genomic data files: genotype files in the AlphaGenes format (``-geno_file``), sequence allele read in the AlphaGenes format (``-seq_file``), and binary Plink files (``-plink_file``). Use of binary Plink files requires the package ``alphaplinkpython``, which  can be installed via ``pip``, but is only stable for Linux. There are known issues with this package, so we do not advocate its use at the moment.

Use the ``-start_snp`` and ``-stop_snp`` to run the analysis only on a subset of markers.

The input options in the form of ``[xxx ...]`` can take in more than one input file separated by space.

Output Arguments 
----------------

::

    Output options:
      -out_file PREFIX      The output file prefix. All file outputs will be stored
                            as "PREFIX.dosage.txt" and so on.
      -out_id_order OUT_ID_ORDER
                            Determines the order in which individuals are ordered
                            in the output file based on their order in the
                            corresponding input file. Individuals not in the input
                            file are placed at the end of the file and sorted in
                            alphanumeric order. These individuals can be suppressed
                            with the "-out_id_only" option. Options: id, pedigree,
                            genotypes, sequence, segregation. Default: id.
      -out_id_only          Flag to suppress the individuals not present in
                            the file used with "-out_id_order". It also suppresses "dummy"
                            individuals.
      -n_io_thread N_IO_THREAD
                            Number of threads to use for input/output. Default: 1.


    Peeling output options:
      -no_dosage            Flag to suppress the dosage files.
      -no_param             Flag to suppress writing the model parameter files.
      -alt_allele_prob      Flag to write out the alternative allele frequencies for each metafounder.
      -seg_prob             Flag to enable writing out the segregation probabilities.
      -phased_geno_prob     Flag to enable writing out the phased genotype probabilities.
      -geno_prob            Flag to enable writing out the genotype probabilities.
      -hap                  Flag to call and write out the haplotypes.
      -geno                 Flag to call and write out the genotypes.
      -geno_threshold [GENO_THRESHOLD ...]
                            Genotype calling threshold(s). Multiple space separated values allowed.
                            Value less than 1/3 will be replaced by 1/3.
      -hap_threshold [HAP_THRESHOLD ...]
                            Haplotype calling threshold(s). Multiple space separated values allowed.
                            Value less than 1/2 will be replaced by 1/2.
      -binary_call_file    Flag to write out the called genotype files as a
                            binary plink output [Not yet implemented].

By default |Software| produces a dosage file and two model parameter files (genotype error rate and recombination rate). Creation of these files can be suppressed with the ``-no_dosage``, and ``-no_param`` options. |Software| can also write out the alternative allele frequencies per metafounder (*.alt_allele_prob.txt*) with ``-alt_allele_prob`` argument, the phased genotype probability file (*.phased_geno_prob.txt*) with the ``-phased_geno_prob`` argument and the segregation probability file (*.seg_prob.txt*) with the ``-seg_prob`` argument.

The ``-geno_threshold`` and ``-hap_threshold`` arguments respectively control control which genotypes and haplotypes are called. A threshold of 0.9 will give calls only if the probability mass for one genotype (or haplotype) is higher than 0.9. Using a higher-value will increase the accuracy of called genotypes (or haplotypes), but will result in fewer called genotypes (or haplotypes). Since there are three genotypes states and two haplotype states, "best-guess" genotypes and haplotypes are respectively called with a threshold less than ``1/3`` and ``1/2``.

``-binary_call_file`` option can be used to change the output to a plink binary format.

The order in which individuals are output can be changed by using the ``out_id_order`` option. This option changes the order in which individuals are written out to the order in which they were observed in the corresponding file. The ```-out_id_only`` option suppresses the output of dummy individuals (not recommended for hybrid peeling).

The argument ``-n_io_thread`` controls the number of threads/processes used by |Software|. |Software| uses additional threads to parse and format input and output files. Setting this option to a value greater than 1 is only recommended for very large files (i.e. >10,000 individuals).

Peeling arguments 
------------------

::

    Mandatory peeling arguments:
      -method METHOD        Program run type. Either "single" or "multi".
    
    Optional peeling arguments:
      -n_cycle N_CYCLE    Number of peeling cycles. Default: 5.
      -n_thread N_THREAD
                            Number of threads to use. Default: 1.
      -rec_length REC_LENGTH
                            Estimated recombination length of the chromosome in Morgans.
                            [Default 1.00]

    Peeling control arguments:
      -est_geno_error_prob  Flag to re-estimate the genotyping error rates after
                            each peeling cycle.
      -est_seq_error_prob   Flag to re-estimate the sequencing error rates after
                            each peeling cycle.
      -est_rec_prob         Flag to re-estimate the recombination rates after
                            each peeling cycle.
      -est_alt_allele_prob  Flag to estimate the alternative allele frequencies using
                            all observed genotypes prior to peeling.
      -update_alt_allele_prob
                            Flag to re-estimate the alternative allele frequencies
                            for each metafounder after each peeling cycle.
      -no_phase_founder     A flag phase a heterozygous allele in one of the
                            founders (if such an allele can be found).
      -sex_chrom            A flag to indicate that input data is for a sex chromosome. Sex needs to
                            be given in the pedigree file. This is currently an
                            experimental option.

    Genotype probability arguments:
      -geno_error_prob GENO_ERROR_PROB
                            Genotyping error rate. [Default 0.0001]
      -seq_error_prob SEQ_ERROR_PROB
                            Sequencing error rate. [Default 0.001]

``-method`` controls whether the program is run in "single-locus" or "multi-locus" model. Single locus mode does not use linkage information to perform imputation. It is fast, but not very accurate. Multi-locus mode runs multi-locus iterative peeling which uses linkage information to increase accuracy and calculate segregation values.

For hybrid peeling, where a large amount (millions of segregating sites) of sequence allele read counts needs to be imputed, first run the program in multi-locus mode to generate a segregation file, and then run the program in single-locus mode with a known segregation file.

The ``-geno_error_prob``, ``-seq_error_prob`` and ``-rec_length`` arguments control some of the model parameters used in the model. ``-seq_error_prob`` must not be zero. |Software| is robust to deviations in genotyping error rate and sequencing error rate so it is not recommended to use these options unless large deviations from the default are known. Changing the ``-length`` argument to match the genetic map length can increase accuracy in some situations.

The ``-est_geno_error_prob`` and ``-est_seq_error_prob`` options estimate the genotyping error rate and the sequencing error rate based on miss-match between observed and inferred states. This option is generally not necessary and can increase runtime. ``-est_alt_allele_prob`` estimates the alternative allele frequencies before peeling using all available observed genotypes. This option can be useful if there are a large number of non-genotyped founders. ``-update_alt_allele_prob`` re-estimates the alternative allele frequencies per metafounder after each peeling cycle using the inferred genotype probabilities of the founders. For implementation of metafounders (**without** ``-alt_allele_prob_file``), both ``-est_alt_allele_prob`` and ``-update_alt_allele_prob`` should be used.

Hybrid peeling arguments 
------------------------

::

    Single locus arguments:
      -seg_file SEG_FILE    A segregation probabilities file for hybrid peeling.
      -seg_map_file SEG_MAP_FILE
                            A map file for loci in the segregation probabilities file.
      -map_file MAP_FILE    A map file for all loci in hybrid peeling.

In order to run hybrid peeling the user needs to supply a ``-map_file`` which gives the genetic positions for the SNPs in the sequence allele read counts data supplied, a ``-seg_map_file`` which gives the genetic position for the SNPs in the segregation file, and a ``-seg_file`` which gives the segregation values generated via multi-locus iterative peeling. These arguments are not required for running in multi-locus mode.

============
File formats
============

Input file formats
------------------

Pedigree file
=============

Each line of a pedigree file has three values, the individual's id, their father's id, and their mother's id. "0" represents an unknown id. Individuals with one unknown parent get internally assigned a dummy/unknown parent. Hence all individuals have both or none parents known. Individuals with two unknown parents are considered as founders and are internally allocated to a metafounder (unknown parent group) ``"MF_1"`` (or defined by the user through ``-main_metafounder``). Users can provide additional metafounders as shown below - these must start with ``"MF_"``.

Example:

::

  id1 0 0
  id2 0 0
  id3 id1 id2
  id4 id1 id2

or

::

  id1 MF_1 MF_1
  id2 MF_2 MF_2
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

Binary Plink files are supported using the package ``AlphaPlinkPython``. The pedigree supplied by the *.fam* file will be used if a pedigree file is not supplied. Otherwise, the pedigree file will be used and the *.fam* file will be ignored.

Map file 
========

The map file gives the chromosome number, the marker name, and the base pair position for each marker in two columns. Only markers on one chromosome should be provided! 

Example:

::

  1 snp_a 12483939
  1 snp_b 192152913
  1 snp_c 65429279
  1 snp_d 107421759

Alternative Allele Probability File
===================================

The alternative allele probability file allows for user-defined population alternative allele probabilities. This file contains the metafounder group denoted MF_x, where x is by default "1" but see ``-main_metafounder``, followed by alternative allele probabilities for all the markers. In case of multiple metafounders, provide multiple rows in the file. The default starting alternative allele probabilities are 0.5 for each marker. If you don't have information for some markers, provide 0.5 for these in the file.

Example:

::

  MF_1 0.30 0.21 0.44 0.24

Or

::

  MF_1 0.30 0.21 0.44 0.24
  MF_2 0.40 0.34 0.25 0.40

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

The haplotype file (*.phased_geno_prob.txt*) provides the (phased) allele probabilities for each locus. There are four lines per individual containing the allele probability for the (aa, aA, Aa, AA) alleles where the paternal allele is listed first, and where *a* is the reference (or major) allele and *A* is the alternative (or minor) allele.

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

Model parameter files
=====================

|Software| outputs four model parameter files: *.alt_allele_prob.txt*, *.seq_error_prob.txt*, *.geno_error_prob.txt*, *.rec_prob.txt*. These give the alternative allele frequency, sequencing error rates, genotyping error rates and the recombination rates used. In the *.alt_allele_prob.txt*, there is a column per metafounder with an alternative allele frequency for each marker. Here, all values will range from 0.01 to 0.99. The other three files contain a single column with an entry for each marker. By default, |Software| will output *.seq_error_prob.txt*, *.geno_error_prob.txt* and *.rec_prob.txt*. The *.alt_allele_prob.txt* will only be outputted with the argument ``-alt_allele_prob``.

Example ``.alt_allele_prob.txt`` file for two metafounders and four loci:

::

  MF_1      MF_2
  0.468005  0.390000
  0.195520  0.219890
  0.733061  0.509849
  0.145847  0.090000


.. |Software| replace:: AlphaPeel
