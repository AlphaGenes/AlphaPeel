===============
Program options
===============

|Software| takes in a number of command line arguments to control the program's behavior. To view a list of arguments, run |Software| without any command line arguments, i.e. ``AlphaPeel`` or ``AlphaPeel -h``. 


Core Arguments 
--------------

::
  
  Core arguments
    -out prefix              The output file prefix.

The ``-out`` argument gives the output file prefix for where the outputs of |Software| should be stored. By default, |Software| outputs a file with imputed genotypes, ``prefix.genotypes``, phased haplotypes ``prefix.phase``, and genotype dosages ``prefix.dosages``. For more information on which files are created, see "Output Arguments", below.


Input Arguments 
----------------

::

    Input Options:
      -bfile [BFILE [BFILE ...]]
                          File(s) in plink (binary) format. Only stable on
                          Linux).
      -genotypes [GENOTYPES [GENOTYPES ...]]
                          File(s) in AlphaGenes format.
      -seqfile [SEQFILE [SEQFILE ...]]
                          Sequence data file(s).
      -pedigree [PEDIGREE [PEDIGREE ...]]
                          Pedigree file(s) in AlphaGenes format.
      -startsnp STARTSNP    The first marker to consider. The first marker in the
                          file is marker "1".
      -stopsnp STOPSNP      The last marker to consider.

|Software| requires a pedigree file and one or more genotype files to run the analysis.

|Software| supports binary plink files, ``-bfile``, genotype files in the AlphaGenesFormat, ``-genotypes``, and sequence data read counts in the AlphaGenes format, ``-seqfile``. A pedigree file must be supplied using the ``-pedigree`` option. 

Use the ``-startsnp`` and ``-stopsnp`` comands to run the analysis only on a subset of markers.

The input options in the form of ``[xxx [xxx ...]]`` can take in more than one input file that are seperated by space.

Binary plink files require the package ``alphaplinkpython``. This can be installed via ``pip`` but is only stable for Linux.


Output Arguments 
----------------
::

    Output options:
      -writekey WRITEKEY    Determines the order in which individuals are ordered
                            in the output file based on their order in the
                            corresponding input file. Animals not in the input
                            file are placed at the end of the file and sorted in
                            alphanumeric order. These animals can be suppressed
                            with the "-onlykeyed" option. Options: id, pedigree,
                            genotypes, sequence, segregation. Defualt: id.
      -onlykeyed            Flag to suppress the animals who are not present in
                            the file used with -outputkey. Also suppresses "dummy"
                            animals.
      -iothreads IOTHREADS  Number of threads to use for io. Default: 1.


    Peeling output options:
      -no_dosages           Flag to suppress the dosage files.
      -no_seg               Flag to suppress the segregation files (e.g. when
                            running for chip imputation and not hybrid peeling).
      -no_params            Flag to suppress writing the parameter files.
      -haps                 Flag to enable writing out the genotype probabilities.
      -calling_threshold [CALLING_THRESHOLD [CALLING_THRESHOLD ...]]
                            Genotype calling threshold(s). Multiple space
                            separated values allowed. Use. .3 for best guess
                            genotype.
      -binary_call_files    Flag to write out the called genotype files as a
                            binary plink output [Not yet implemented].

By default |Software| produces a dosages file, a segregation files and two parameter files (genotyping error and recombination rate). Creation of each of these files can be suppressed with the ``-no_dosages``, ``-no_seg``, and ``-no_params`` options. |Software| can also write out the genotype probability file (.haps) with the `-haps` argument.

The ``-calling_threshold`` arguments controls which genotypes (and phased haplotypes) are called as part of the algorithm. A calling threshold of 0.9 indicates that genotypes are only called if greater than 90% of the final probability mass is on that genotype. Using a higher-value will increase the accuracy of called genotypes, but will result in fewer genotypes being called. Since there are three genotypes states,  "best-guess" genotypes are produced with a calling threshold less than ``0.33``. ``-calling_threshThe ``-binary_call_files`` option can be used to change the output to a plink binary format. 

The order in which individuals are output can be changed by using the ``writekey`` option. This option changes the order in which individuals are written out to the order in which they were observed in the corresponding file. The ```-onlykeyed`` option suppresses the output of dummy individuals (not recommended for hybrid peeling). 

The parameter ``-iothreads`` controls the number of threads/processes used by |Software|. |Software| uses additional threads to parse and format input and output files. Setting this option to a value greater than 1 is only recommended for very large files (i.e. >10,000 individuals).


Peeling arguments: 
------------------------
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

For hybrid peeling, where a large amount (millions of segregating sites) of sequence data needs to be imputed, first run the program in multi-locus mode to generate a segregation file, and then run the program in single-locus mode with a known segregation file.

The ``-error``, ``-seqerror`` and ``-length`` arguments control some of the parameters used in the model. |Software| is robust to deviations in genotyping error rate and sequencing error rate so it is not recommended to use these options unless large deviations from the default are known. Changing the ``-length`` argument to match the genetic map length can increase accuracy in some situations.

The ``-esterrors`` option estimated the genotyping error rate based on observed information, this option is generally not necessary and can increase runtime. ``-estmaf`` estimates the minor allele frequency after each peeling cycle. This option can be useful if there are a large number of non-genotyped founders. 


Hybrid peeling arguments 
-----------------------------
::

    Single locus arguments:
      -mapfile MAPFILE      A map file (chr marker_name position) for genotype data.
      -segmapfile SEGMAPFILE
                            a map file for the segregation estimates for hybrid
                            peeling.
      -segfile SEGFILE      A segregation file for hybrid peeling.

In order to run hybrid peeling the user needs to supply a ``-mapfile`` which gives the genetic positions for the SNPs in the sequence data supplied, a ``-segmapfile`` which gives the genetic position for the SNPs in the segregation file, and a ``-segfile`` which gives the segregation values generated via multi-locus iterative peeling. These arguments are not required for running in multi-locus mode.

.. |Software| replace:: AlphaPeel
