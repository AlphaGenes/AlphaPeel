.. _usage:

-----
Usage
-----

.. _program_options:

===============
Program options
===============

|Software| accepts several command-line options to control the program's behaviour.
To view a list of all supported options, run |Software| like this:
``AlphaPeel`` or ``AlphaPeel -h``.

.. _input_options:

Input options
-------------

.. parsed-literal::

    Individuals:
      -ped_file [PED_FILE ...]
                          Pedigree file(s)
                          (:ref:`see format details <ped_file_format>`).
      -geno_file [GENO_FILE ...]
                          Genotype file(s)
                          (:ref:`see format details <geno_file_format>`).
      -seq_file [SEQ_FILE ...]
                          Sequence allele read count file(s)
                          (:ref:`see format details <seq_file_format>`).
      -x_chr              Indicate that input data is for the :ref:`X chromosome <zero_one_two_etc>`.
      -pheno_file [PHENO_FILE ...]
                          Phenotype file(s)
                          (:ref:`see format details <pheno_file_format>`).
      -plink_file [PLINK_FILE ...]
                          Plink (binary) file(s)
                          (:ref:`see format details <plink_file_format>`).

    :ref:`Markers <markers>`:
      -map_file MAP_FILE  Map file for :ref:`loci <markers>` in genomic data files
                          (:ref:`see format details <map_file_format>`).
      -start_snp START_SNP
                          The first :ref:`locus <markers>` to consider.
                          Counting starts at 1.
                          Default: 1.
      -stop_snp STOP_SNP  The last :ref:`locus <markers>` to consider.
                          Default: all loci in input files.

    :ref:`Model parameters <prob_freq_rate>` and other:
      -alt_allele_prob_file [ALT_ALLELE_PROB_FILE ...]
                          Alternative allele :ref:`probability <prob_freq_rate>` file
                          (:ref:`see format details <alt_allele_prob_file_format>`).
                          Default: 0.5 for each locus.
      -main_metafounder MAIN_METAFOUNDER
                          Metafounder name for unknown parents
                          (:ref:`see format details <ped_file_format>`).
                          Default: MF_1.
      TODO: how is rec length used?
      -rec_length REC_LENGTH
                          Recombination length of the chromosome in Morgans.
                          Default: 1.00.
      TODO: rename mutation_rate to mut_prob
      -mutation_rate MUTATION_RATE
                          Mutation :ref:`probability <prob_freq_rate>`.
                          Default: 1e-8.
      -geno_error_prob GENO_ERROR_PROB
                          Genotype error :ref:`probability <prob_freq_rate>`.
                          Default: 0.0001.
      -seq_error_prob SEQ_ERROR_PROB
                          Sequence error :ref:`probability <prob_freq_rate>`.
                          Must not be 0.
                          Default: 0.001.
      -pheno_penetrance_prob_file
                          [PHENO_PENETRANCE_PROB_FILE ...]
                          Phenotype penetrance :ref:`probability <prob_freq_rate>` file
                          (:ref:`see format details <pheno_penetrance_prob_file_format>`).

|Software| requires a pedigree file (``-ped_file``) and
one or more genomic data files to run the analysis.

|Software| supports the following genomic data files:
genotype file in the AlphaGenes format (``-geno_file``),
sequence allele read count file in the AlphaGenes format (``-seq_file``), and
binary Plink file (``-plink_file``).

When ``-x_chr`` is used the :ref:`pedigree file <ped_file_format>` and
:ref:`genotype file <geno_file_format>` have specific requirements.
Follow the links to respective file formats for more details.

|Software| also supports phenotype files (``-pheno_file``) and
corresponding phenotype penetrance files (``-pheno_penetrance_prob_file``).
Both must be provided for work with phenotypes.

The input options in the form of ``-opt [XYZ ...]``
can accept more than one argument separated by spaces.

Use ``-start_snp`` and ``-stop_snp``
to run the analysis only on a subset of :ref:`markers <markers>`.

.. _markers:

.. note::

    We use interchangeably the terms "marker(s)", "locus/loci", "SNP(s)",
    or "site(s)", to refer to a specific position in the genome,
    where we typically observe polymorphism in the population.

|Software| supports specifying a number of model parameters,
which are :ref:`probabilities <prob_freq_rate>` of different events or outcomes:
alternative allele probabilities in the base population(s) through a file (``-alt_allele_prob_file``),
recombination length of the chromosome (``-rec_length``),
TODO: rename mutation_rate to mut_prob
mutation rate (``-mutation_rate``),
genotype error probability (``-geno_error_prob``),
sequence error probability (``-seq_error_prob``), and
phenotype penetrance probability through a file (``-pheno_penetrance_prob_file``).

.. _robust_parameters:

.. note::

    The accuracy of |Software| results has been shown
    to be quite robust to deviations in most of the model parameters,
    so it might not be needed to change or estimate them from the input data;
    at least unless large deviations from the defaults are known or expected.
    Having said this, do explore what works best for your data and your aims!

.. _prob_freq_rate:

.. note::

    We use the term "probability" to also represent the commonly used terms
    "frequency" and "rate", to refer to the same concept of a value
    between 0 and 1 that quantifies the likelihood or proportion
    of a certain event or outcome.

.. _output_options:

Output options
--------------

.. parsed-literal::

    Individuals:
      -no_dosage            Suppress default output of :ref:`allele dosages <dosage_file_format>`.
      -geno                 Call and output :ref:`genotypes <genotype_file_format>`.
      -geno_threshold [GENO_THRESHOLD ...]
                            Genotype calling threshold(s) from the genotype probabilities.
                            Multiple space separated values allowed.
                            Value(s) less than 1/3 are replaced by 1/3.
                            Default: 1/3.
      -geno_prob            Output :ref:`genotype probabilities <geno_prob_file_format>`.
      -phased_geno_prob     Output :ref:`phased genotype probabilities <phased_geno_prob_file_format>`.
      -hap                  Call and output :ref:`haplotypes <hap_file_format>`.
                            Default: 1/2.
      -hap_threshold [HAP_THRESHOLD ...]
                            Haplotype calling threshold(s) from the phased genotype probabilities
                            Multiple space separated values allowed.
                            Value(s) less than 1/2 are replaced by 1/2.
      -seg_prob             Output :ref:`segregation probabilities <seg_prob_file_format>`.
      -pheno_prob           Output :ref:`phenotype probabilities <pheno_prob_file_format>`.

    Prefix, order, and IO:
      -out_file PREFIX      The output file prefix. All file outputs will be named
                            as "PREFIX.OUTPUT.txt", where "OUTPUT" is the type of output
                            (for example, "dosage" and "geno_prob").
      -out_id_order OUT_ID_ORDER
                            Determines the order of individuals in the output
                            file based on their order in the
                            corresponding input file. Individuals not in the input
                            file are placed at the end of the file and sorted in
                            alphanumeric order. These individuals can be suppressed
                            with the -out_id_only option. Accepted arguments for
                            this option are: id, pedigree, genotypes, sequence, and
                            segregation.
                            Default: id.
      -out_id_only          Suppress output for individuals not present in
                            the file specified with -out_id_order. It also suppresses
                            "dummy" individuals.
      -n_io_thread N_IO_THREAD
                            Number of threads to use for input/output (IO).
                            Default: 1.

|Software| by default produces a :ref:`dosage file <dosage_file_format>`.
Additional individual-level outputs can be requested with the options described above.
The ``-geno_threshold`` and ``-hap_threshold``
respectively control which genotypes and haplotypes are called.
A threshold of 0.9 gives calls only if
the probability for one genotype (or haplotype allele) is higher than 0.9.
When the probability is lower than the threshold,
the output is set to value :ref:`9 (missing) <zero_one_two_etc>`.
Using a higher value will increase the accuracy of called genotypes (or haplotypes),
but will result in fewer calls.
Since there are three genotype states and two haplotype states,
"best-guess" genotypes and haplotypes are
respectively called with a threshold less than ``1/3`` and ``1/2``.

The output order of individuals can be changed
using the ``-out_id_order`` option,
with additional control provided with the ``-out_id_only`` option.
The latter option is not recommended for hybrid peeling or
any combination of different input files.

The option ``-n_io_thread`` controls the number of threads/processes
used by |Software|.
|Software| uses additional threads to parse and format input and output files.
Setting this option to a value greater than 1 is only recommended for very large files
(i.e. >10,000 individuals).

.. _peeling_methods:

Peeling methods
---------------

.. parsed-literal::

    Strategy:
      -method METHOD        Peeling method: single or multi.
                            Default: multi.

    Single-locus options for the second stage of hybrid peeling:
      -seg_file SEG_FILE    :ref:`Segregation probabilities file <seg_prob_file_format>`.
      -seg_map_file SEG_MAP_FILE
                            Map file for loci in the segregation probabilities file.

|Software| supports three peeling strategies: single-locus, multi-locus, and hybrid.

Single-locus peeling method does not use linkage information between loci in iterative peeling.
It is fast, but not very accurate.

Multi-locus peeling method runs multi-locus iterative peeling,
which uses linkage information to increase accuracy and calculate segregation probabilities,
but it is much slower than single-locus method.

Hybrid peeling is useful in settings with
a SNP genotypes with a limited number of markers and
a sequence allele read counts from a large number of loci.
In this setting, you can
first run the multi-locus peeling method on
SNP genotypes to estimate segregation probabilities, and
then run the single-locus peeling method on
sequence allele read counts and the segregation probabilities.
In this second stage of hybrid peeling, provide:
a ``-map_file`` with positions for loci in the sequence allele read count data,
a ``-seg_file`` with segregation probabilities generated via multi-locus method, and
a ``-seg_map_file`` with genetic positions for loci in the segregation probabilities file.
This combination of options is not required in the standard multi-locus mode.

.. _peeling_parameters:

Peeling parameters
------------------

.. parsed-literal::

    Computational parameters:
      -n_cycle N_CYCLE      Number of peeling cycles.
                            Default: 5.
      -n_thread N_THREAD
                            Number of threads to use.
                            Default: 1.

    Estimation of model parameters:
      TODO: rename to est_start_alt_allele_prob
      -est_alt_allele_prob  Estimate starting alternative allele probabilities using
                            all inputted genomic data prior to peeling.
      TODO: rename to est_alt_allele_prob
      -update_alt_allele_prob
                            Estimate :ref:`alternative allele probabilities <alt_allele_prob_file_format>`
                            for each metafounder after each peeling cycle.
      TODO: -est_rec_prob is not recognised
      -est_rec_prob         Estimate :ref:`recombination probabilities <rec_prob_file_format>`
                            after each peeling cycle.
      TODO: no_phase_founder is not documented - talking to Ros she noticed that
      the default behaviour is to phase heterozygous genotypes in founders
      in some way (TODO: what way - she noticed it is not driven by data possibly?!)
      and this flag suppresses this behaviour.
      -no_phase_founder     Phase a heterozygous allele in one of the
                            founders (if such an allele can be found).
      -est_geno_error_prob  Estimate :ref:`genotype error probability <geno_error_prob_file_format>`
                            after each peeling cycle.
      -est_seq_error_prob   Estimate :ref:`sequence error probability <seq_error_prob_file_format>`
                            after each peeling cycle.
      -est_pheno_penetrance_prob
                            Estimate :ref:`phenotype penetrance probabilities <pheno_penetrance_prob_file_format>`
                            after each peeling cycle.

Computational effort and speed of |Software| can be controlled with
the number of peeling cycles (``-n_cycle``,
increasing the number will marginally increase accuracy, but also runtime) and
the number of threads (``-n_thread``, to reduce runtime on large datasets).

TODO: delete these options as they are just making a mess - below is a clear and simple behaviour
      -no_param             Suppress output of model parameter files.
      -alt_allele_prob      Output :ref:`alternative allele probabilities <alt_allele_prob_file_format>`.

|Software| can estimate the model parameters from the input data.
The :ref:`default or user provided input values<input_options>`
are used as a starting point for estimation.
See a note on :ref:`robustness of results <robust_parameters>`
to these parameters.
TODO: check for this behaviour for each parameter (output only when estimated)
When estimation options are used,
the respective parameters are estimated after each peeling cycle and
output to a file at the end.
This process usually increases running times and
might require additional peeling cycles to converge.
The estimates are based on inferred states of the modelled events and
their match/mismatch between observed and inferred states.

Alternative allele probabilities in the founders are estimated
(using ``-est_alt_allele_prob``)
as half of the mean of estimated :ref:`allele dosage <dosage_file_format>`
in the founders (potentially grouped into multiple populations via
metafounders).
The estimates are constrained to be between 0.01 and 0.99
to avoid TODO: discuss with Evie how to word this.
This estimation can be initiated with an estimate from inputted genomic data
(using ``-est_init_alt_allele_prob``).
This option uses Newton optimisation, which also requires starting values.
These starting values are by default 0.5,
but can also be provided by the user using ``-alt_allele_prob_file``.
Note that this estimation is not taking the pedigree structure into account,
so it is a naive population estimate and does not pertain to founders of the pedigree
and ignores metafounders.

For a pedigree with multiple metafounders,
the user has three options to obtain metafounder-specific alternative allele probabilities:
TODO: rename est_alt_allele_prob to est_start_alt_allele_prob
TODO: rename to update_alt_allele_prob to est_alt_allele_prob
(1) use the default starting value of 0.5 for all loci and
``-update_alt_allele_prob``,
(2) use ``-est_alt_allele_prob`` to get more informed starting value and
``-update_alt_allele_prob``, or
(3) provide starting values using ``-alt_allele_prob_file`` and
then potentially ``-update_alt_allele_prob`` from the inputted genomic data.

TODO: Say something about -est_rec_prob

TODO: no_phase_founder is not documented / see above

Error probabilities (using ``-est_geno_error_prob`` and ``-est_seq_error_prob``)
are estimated as the proportion of mismatches between observed and inferred states.

TODO: Say something about penetrance probabilities estimation

.. _file_formats:

============
File formats
============

.. _input_file_formats:

Input file formats
------------------

.. _ped_file_format:

Pedigree file
=============

This file has one line with
*recorded parents* for each individual.
The file(s) should include all individuals present in other files
and their known relatives present only in pedigree.
Each line of a *pedigree* file has three values,
the individual's ID,
their father's ID, and
their mother's ID.
Unknown parent ID is encoded with ``0``.
Individuals with one unknown parent are internally assigned a "dummy" parent.
Therefore, all individuals have either both parents known or neither parent known.
Individuals with two unknown parents are considered as founders and
are internally assigned to a metafounder (unknown parent group) ``MF_1``
(or defined by the user through ``-main_metafounder``).
Users can provide additional metafounders as shown below;
these must start with ``MF_``.

When working with the X chromosome,
a fourth value is needed for individual's sex:
``0`` for males and ``1`` for females.

Example with four individuals across two generations:

::

  id1 0 0
  id2 0 0
  id3 id1 id2
  id4 id1 id2

Example with two metafounders:

::

  id1 MF_1 MF_1
  id2 MF_2 MF_2
  id3 id1 id2
  id4 id1 id2

Example with sex information;
id1 and id3 are males, while id2 and id4 are females:

::

  id1 0 0 0
  id2 0 0 1
  id3 id1 id2 0
  id4 id1 id2 1

.. _geno_file_format:

Genotype file
=============

This file has one line with
*observed genotypes* for each genotyped individual.
The file does not need to include all individuals present in other files.
The first value in each line is the individual's ID.
The remaining values are observed genotypes at each locus.
Only loci on one chromosome should be provided!

.. _zero_one_two_etc:

.. note::

    |Software| works with bi-allelic loci with two *alleles*;
    *reference (major)* allele ``a`` and
    *alternative (minor)* allele ``A``.
    The terms major and minor indicates their frequency in a population,
    but this is not a requirement for |Software|.
    These alleles are numerically encoded as ``0`` and ``1``, respectively.
    Combining two alleles in a diploid individual gives three possible *genotypes*:
    reference homozygote ``a/a``,
    heterozygote ``a/A`` or ``A/a``, and
    alternative homozygote ``A/A``.
    When the origin of alleles that an individual inherited is known,
    we have four *phased genotypes*: ``aa``, ``aA``, ``Aa``, and ``AA``,
    where the *paternal allele* is listed first and
    the *maternal allele* is listed second.
    These genotypes are numerically encoded as ``0``, ``1``, and ``2``, respectively.
    Missing alleles and genotypes are numerically encoded as ``9``.
    The numerical codes are called *allele dosages*, because
    they represent the number (dose) of alternative alleles.

    When working with the X chromosome:
    (1) *heterogametic genotypes* (for males in mammals (XY) and for females in birds (ZW))
    should be coded as:
    ``0`` (reference allele ``a`` on the X chromosome of the XY genotype),
    ``1`` (alternative allele ``A`` on the X chromosome of the XY genotype), or
    ``9`` (missing) and
    (2) *Homogametic genotypes* (females in mammals (XX) and males in birds (ZZ))
    should be coded as described above for autosomes
    (since they have the XX genotype).

Example with four individuals and their genotypes at four loci:

::

  id1 0 2 9 0
  id2 1 1 1 1
  id3 2 0 2 0
  id4 0 2 1 0

Example with four individuals and their X chromosome genotypes at four loci;
id1 and id3 are males, while id2 and id4 are females:

::

  id1 0 1 9 0
  id2 1 1 1 1
  id3 1 0 1 0
  id4 0 2 1 0

.. _seq_file_format:

Sequence allele read counts file
================================

This file has two lines with
*observed allele read counts* for each sequenced individual.
The file does not need to include all individuals present in other files.
The first line gives the individual's ID and read counts for the :ref:`reference allele <zero_one_two_etc>` at each locus.
The second line gives the individual's ID and allele read counts for the :ref:`alternative allele <zero_one_two_etc>` at each locus.
TODO: work with Yuni on X seq data input
The same format works for the autosomes and the X chromosome.
Only loci on one chromosome should be provided!

Example with four individuals and their allele read counts at four loci:

::

  id1 4 0 0 7 # Reference allele read counts for id1
  id1 0 3 0 0 # Alternative allele read counts for id1
  id2 1 3 4 3
  id2 1 1 6 2
  id3 0 3 0 1
  id3 5 0 2 0
  id4 2 0 6 7
  id4 0 7 7 0

.. _pheno_file_format:

Phenotype file
==============

This file has one or multiple lines with
*observed phenotypes* for each phenotyped individual.
Multiple lines per individual support repeated phenotyping.
The file does not need to include all individuals present in other files.
The first value in each line is the individual's ID.
The remaining values are the phenotypes of the individual for a specific trait,
coded from ``0`` onwards.
For example, a binary trait should be coded as ``0`` and ``1``,
while a trait with three states should be coded as ``0``, ``1``, and ``2``
(current phenotype functionality works only with one locus and one trait).

Example with four individuals and their phenotypes for a binary trait:

::

  id1 0
  id2 1
  id3 1
  id4 0

.. _plink_file_format:

Binary Plink file
=================

This file is supported through the package ``AlphaPlinkPython``,
which can be installed via ``pip``, but is only stable for Linux.
There are known issues with this package,
so we do not advocate its use at the moment.
The pedigree supplied by the ``.fam`` file will be used if a pedigree file is not supplied.
Otherwise, the pedigree file will be used and the ``.fam`` file will be ignored.

.. _map_file_format:

Map file
========

This file has one line with
*genome information* for each :ref:`locus <markers>`.
The file should include all loci present in other files.
Each line of a *map* file has three values,
the chromosome number,
the locus name, and
the base-pair position.
Only loci on one chromosome should be provided!

TODO: mention somewhere how we convert base-pair positions to genetic positions (Morgans) and interaction? with -rec_length

Example with four SNP loci on chromosome 1:

::

  1 snp_a 12483939
  1 snp_b 192152913
  1 snp_c 65429279
  1 snp_d 107421759

TODO: how does seg_map_file look like?

.. _output_file_formats:

Output file formats
-------------------

.. _dosage_file_format:

Dosage file
===========

The ``.dosage.txt`` file contains *expected dosages for the alternative allele*
for each individual.
These expected dosages are obtained by
(1) taking the :ref:`allele dosages <zero_one_two_etc>` (``0``, ``1``, and ``2``)
of each individual at each locus,
(2) weighting them with corresponding :ref:`genotype probabilities <geno_prob_file_format>`, and
(3) summing the weighted dosages.
There is one line per individual.
The first value in each line is the individual ID.
The remaining values are allele dosages at each locus.
These values will be between ``0`` and ``2``, inclusive
(:ref:`see the note on encoding alleles and genotypes <zero_one_two_etc>`).

Example with four individuals and four loci:

::

  id1 0.2089 1.7820 1.4725 0.0000
  id2 1.0000 1.0000 1.0000 0.9999
  id3 0.8293 1.1329 1.9999 0.0000
  id4 0.0001 1.9999 1.0000 0.0000

.. _genotype_file_format:

Genotype file
=============

The ``.geno_THRESHOLD.txt`` file contains *called genotypes* for each individual
for each specified threshold (``-geno_threshold THRESHOLD``).
There is one line per individual.
The first value in each line is the individual ID.
The remaining values are called genotypes at each locus,
encoded as ``0``, ``1``, and ``2``, or ``9`` when
the probability is to low to make the call
(:ref:`see the note on encoding alleles and genotypes <zero_one_two_etc>`).

Example with four individuals and four loci:

::

  id1 0 2 1 0
  id2 1 1 1 1
  id3 1 1 2 0
  id4 0 2 1 0

.. _geno_prob_file_format:

Genotype probability file
=========================

The ``.geno_prob.txt`` file contains *genotype probabilities* for each individual.
TODO: remove the empty lines in the output files? (4 will be 3 in text!)
There are four lines per individual
(an empty line and three lines with probabilities for
``aa``, ``aA`` or ``Aa``, and ``AA`` genotypes).
The first value in each line is the individual ID.
The remaining values are genotype probabilities at each locus.

Example with four individuals and four loci:

::

  id1 0.7912 0.0000 0.0000 1.0000
  id1 0.2088 0.2179 0.5274 0.0000
  id1 0.0000 0.7820 0.4725 0.0000
  TODO: remove the empty lines in the output files?
  id2 0.0000 0.0000 0.0000 0.0001
  id2 1.0000 1.0000 1.0000 0.9999
  id2 0.0000 0.0000 0.0000 0.0000

  id3 0.3784 0.2171 0.0000 1.0000
  id3 0.4140 0.4329 0.0001 0.0000
  id3 0.2076 0.3500 0.9999 0.0000

  id4 0.9999 0.0000 0.0000 1.0000
  id4 0.0001 0.0001 0.9999 0.0000
  id4 0.0000 0.9999 0.0000 0.0000

.. _phased_geno_prob_file_format:

Phased genotype probability file
================================

The ``.phased_geno_prob.txt`` file contains *phased genotype probabilities* for each individual.
TODO: remove the empty lines in the output files? (5 will be 4 in text!)
There are five lines per individual
(an empty line and four lines with probabilities for
``aa``, ``aA``, ``Aa``, and ``AA`` phased genotypes).
The first value in each line is the individual ID.
The remaining values are phased genotype probabilities at each locus.

Example with four individuals and four loci:

::

  id1 0.7912 0.0000 0.0000 1.0000
  id1 0.1044 0.1090 0.2637 0.0000
  id1 0.1044 0.1090 0.2637 0.0000
  id1 0.0000 0.7820 0.4725 0.0000
  TODO: remove the empty lines in the output files?
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


.. _hap_file_format:

Phase / haplotype file
======================

The ``.hap_THRESHOLD.txt`` file contains *called haplotypes* for each individual.
There are two lines per individual,
one for each of the paternal and maternal haplotypes of the individual.
The first value in each line is the individual ID.
The remaining values are called alleles at each locus,
encoded as ``0`` and ``1``, or ``9`` when
the probability is to low to make the call
(:ref:`see the note on encoding alleles and genotypes <zero_one_two_etc>`).
For individuals for whom the allele/haplotype of origin can be determined,
the first line provides information on the paternal haplotype, and
the second line provides information on the maternal haplotype.

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

.. _seg_prob_file_format:

Segregation file
================

The ``.seg_prob.txt`` file contains *segregation probabilities* for each individual.
There are four lines per individual,
corresponding to the four possible patterns of segregation (inheritance):

  1. The grand *paternal* allele from the father and
     the grand *paternal* allele from the mother (``Pr(F_p,- & M_p,-)``),
  2. The grand *paternal* allele from the father and
     the grand *maternal* allele from the mother (``Pr(F_p,- & M_-,m)``),
  3. The grand *maternal* allele from the father and
     the grand *paternal* allele from the mother (``Pr(F_-,m & M_p,-)``), and
  4. The grand *maternal* allele from the father and
     the grand *maternal* allele from the mother (``Pr(F_-,m & M_-,m)``).

The symbols above mean the following:
``Pr()`` probability,
``F`` father,
``M`` mother,
``p`` *paternal* allele,
``m`` *maternal* allele, and
``-`` the allele that was not inherited.

TODO: Do we need a diagram for this maybe even show the 4 cases so we bring home the message!?

The first value in each line is the individual ID.
The remaining values are segregation probabilities at each locus.

Example with four individuals and four loci:

::

  id1 0.2500 0.2500 0.2500 0.2500
  id1 0.2500 0.2500 0.2500 0.2500
  id1 0.2500 0.2500 0.2500 0.2500
  id1 0.2500 0.2500 0.2500 0.2500
  id2 0.2500 0.2500 0.2500 0.2500
  id2 0.2500 0.2500 0.2500 0.2500
  id2 0.2500 0.2500 0.2500 0.2500
  id2 0.2500 0.2500 0.2500 0.2500
  id3 0.3356 0.3894 0.5000 0.4390
  id3 0.1644 0.1106 0.0000 0.0610
  id3 0.3356 0.3894 0.5000 0.4390
  id3 0.1644 0.1106 0.0000 0.0610
  id4 0.2046 0.1953 0.2760 0.3914
  id4 0.2954 0.3047 0.2240 0.1086
  id4 0.2046 0.1953 0.2760 0.3914
  id4 0.2954 0.3047 0.2240 0.1086

.. _pheno_prob_file_format:

Phenotype probability file
==========================

The ``.pheno_prob.txt`` file contains *phenotype probabilities* for each individual.
TODO: remove the empty lines in the output files? (1+k will be k in text!)
There are ``1+k`` lines per individual
(an empty line and ``k`` lines with probabilities for ``0:k`` phenotypes).
The first value in each line is the individual ID.
The remaining values are phenotype probabilities
(current phenotype functionality works only with one locus and one trait).

Example with four individuals and two phenotype states for a binary trait:

::

  id1 0.9000
  id1 0.1000
  TODO: remove the empty lines in the output files?
  id2 0.9000
  id2 0.1000

  id3 0.2453
  id3 0.7547

  id4 0.9000
  id4 0.1000

.. _param_file_format:

Parameter file formats
----------------------

.. _alt_allele_prob_file_format:

Alternative allele probability file
===================================

TODO: by row or by col? We think we will do it such that
      metafounders are in columns and loci are in rows.
      So, the task is to use this text below and modify it
      accordingly and then merge the two alt allele prob file format
      section below into one. Ros will lead on this with code
      and documentation. The conclusion from looking across all the files
      was that we will have loci by rows here and most of The
      model parameter files.

This file has one line with
*alternative allele probabilities* for each metafounder.
It provides a way to include known information on allele probabilities
in the pedigree founders (base population).
The file should include all metafounders present in pedigrees.
The first value in each line is the metafounder's ID.
The remaining values are alternative allele probabilities for all loci.
The default alternative allele probabilities are 0.5 for each locus.
If you do not have information for some loci, provide 0.5 for those loci.

Example with one metafounder and four loci:

::

  MF_1 0.30 0.21 0.44 0.24

Example with two metafounders and four loci:

::

  MF_1 0.30 0.21 0.44 0.24
  MF_2 0.40 0.34 0.25 0.40

TODO: See the above duplicate section

Alternative allele probability file
===================================

The ``.alt_allele_prob.txt`` file TODO

TODO: is this the same as the input version? Review with Ros
      :ref:`here <alt_allele_prob_file_format>`.
      Hmm, we should transpose one of these!?
      Best to make all of the inputs & outputs consistent - row-wise or column-wise!

Example with two metafounders and four loci:

::

  MF_1      MF_2
  0.468005  0.390000
  0.195520  0.219890
  0.733061  0.509849
  0.145847  0.090000

.. _rec_prob_file_format:

Recombination probability file
==============================

The ``.rec_prob.txt`` file contains *recombination probabilities*
for each pair of neighbouring loci.
There is one line per one pair of neighbouring loci with one value.

Example with four loci:

::
  0.0000

TODO: -est_rec_prob is not recognised - once it is, fill the example

.. _geno_error_prob_file_format:

Genotype error probability file
===============================

The ``.geno_error_prob.txt`` file contains *genotype error probabilities* for each locus.
There is one line per locus with one value.

Example with four loci:

::

  0.050000
  0.050000
  0.000100
  0.000100

.. _seq_error_prob_file_format:

Sequence error probability file
===============================

The ``.seq_error_prob.txt`` file contains *sequence error probabilities* for each locus.
There is one line per locus with one value.

Example with four loci:

::

  0.000137
  0.000100
  0.000362
  0.000100

.. _pheno_penetrance_prob_file_format:

Phenotype penetrance probability file
=====================================

This file has one line with
*phenotype penetrance probabilities* for each distinct genotype impacting a phenotype.
It provides a way to include known relationship between genotypes and phenotypes.
The file should include all phenotypes present in phenotype files
(current phenotype functionality works only with one locus and one trait).
Specifically, one line contains conditional probabilities of
possible phenotypes for each true *phased genotype*,
including phenotyping errors or other deviations.
Each column represents a different phenotype state
(ordered from ``0`` to the total number of states), and
each row represents the :ref:`phased genotypes <zero_one_two_etc>`
(``aa``, ``aA``, ``Aa``, and ``AA``).

Example for a monogenic recessive binary trait:

::

  0.9 0.1
  0.9 0.1
  0.9 0.1
  0.1 0.9
