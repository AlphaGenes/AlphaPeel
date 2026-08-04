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
``AlphaPeel``, ``AlphaPeel -h``, ``AlphaPeel -help`` or ``AlphaPeel --help``.

User can check the version of the program with ``AlphaPeel -version``. 
Remember to use the correct version of the documentation for the version of the program you are using.
For example, the link to the documentation for version ``v1.3.0`` is https://alphapeel.readthedocs.io/en/v1.3.0/index.html.

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
      -hap_file [HAP_FILE ...]
                          Haplotype file 
                          (:ref:`see format details <hap_file_format>`)
      -x_chr              Indicate that input data is for the :ref:`X chromosome <zero_one_two_etc>`.
      -pheno_file [PHENO_FILE ...]
                          Phenotype file(s)
                          (:ref:`see format details <pheno_file_format>`).
      -phased_geno_prob_file [PHASED_GENO_PROB_FILE ...]
                          Optional external phased_geno_prob_file file(s) 
                          (:ref:`see format details <phased_geno_prob_file_format>`).
                          This will provide the starting internal genotype probability state. 

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
                          ID used to represent the base population of the pedigree
                          (metafounder / unknown parent group)
                          (:ref:`see format details <ped_file_format>`).
                          Default: MF_1.
      -rec_length REC_LENGTH
                          Recombination length of the chromosome in Morgans.
                          Default: 1.00.
      -mut_prob MUT_PROB
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
sequence allele read count file in the AlphaGenes format (``-seq_file``), 
and haplotype file in the AlphaGenes format (``-hap_file``).

.. note::

  Internally, the haplotype input by ``-hap_file`` is used to fill up the missing genotype information,
  and the phasing information is not preserved. If you want to input phased data, 
  consider the ``-phased_geno_prob_file`` input.

|Software| supports genotype probabilities input in the AlphaGenes format (``-phased_geno_prob_file``). 
The ``-phased_geno_prob_file`` is assumed to contain error-free phased genotype probabilities. 
Hence, when using ``-phased_geno_prob_file`` as the input, you should not use with options ``-est_geno_error_prob``, 
``-est_seq_error_prob``. 
The option does not support the ``-x_chr`` option for now. 

When ``-x_chr`` is used the :ref:`pedigree file <ped_file_format>` and
:ref:`genotype file <geno_file_format>` have specific requirements.
Follow the links to file formats for more details.

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
alternative allele probabilities in the base population(s) (``-alt_allele_prob_file``),
recombination length of the chromosome (``-rec_length``),
mutation rate (``-mut_prob``),
genotype error probability (``-geno_error_prob``),
sequence error probability (``-seq_error_prob``), and
phenotype penetrance probabilities (``-pheno_penetrance_prob_file``).

The option ``-rec_length`` provides the software with the 
recombination length of the input chromosome. |Software|
assumes the recombination happens equally likely across
the chromosome. Therefore, the recombination probability
between two loci would be calculated by dividing the 
values of ``-rec_length`` by the distance of the physical
positions of the two loci relative to the total length of the 
chromosome.

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
      -geno                 Call and output :ref:`genotypes <geno_file_format>`.
                            The default genotype calling threshold is set to 1/3.
      -geno_threshold [GENO_THRESHOLD ...]
                            Custom genotype calling threshold(s) from the genotype probabilities.
                            Multiple space separated values allowed.
                            Value(s) less than 1/3 are replaced by 1/3.
      -geno_prob            Output :ref:`genotype probabilities <geno_prob_file_format>`.
      -phased_geno_prob     Output :ref:`phased genotype probabilities <phased_geno_prob_file_format>`.
      -hap                  Call and output :ref:`haplotypes <hap_file_format>`.
                            The default haplotype calling threshold is set to 1/2.
      -hap_threshold [HAP_THRESHOLD ...]
                            Custom haplotype calling threshold(s) from the 
                            phased genotype probabilities.
                            Multiple space separated values allowed.
                            Value(s) less than 1/2 are replaced by 1/2.
      -seg_prob             Output :ref:`segregation probabilities <seg_prob_file_format>`.
      -pheno_prob           Output :ref:`phenotype probabilities <pheno_prob_file_format>`.
      -alt_allele_prob      Output :ref:`alternative allele probabilities <alt_allele_prob_file_format>`. 
                            Output initial value ``0.5`` if none of ``est_start_alt_allele_prob``, 
                            ``est_alt_allele_prob``, or ``alt_allele_prob_file`` is used.
      -pheno_penetrance_prob
                            Output :ref:`phenotype penetrance probabilities <pheno_penetrance_prob_file_format>`.

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
      -out_digits           Specify the number of digits to round the outputs. 
                            Does not apply to outputs from ``alt_allele_prob``, 
                            ``geno_error_prob``, ``seq_error_prob``, 
                            and ``pheno_penetrance``. 
                            Default: 4.

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

We round the thresholds in the output filenames to three digits. 

The output order of individuals can be changed
using the ``-out_id_order`` option,
with additional control provided with the ``-out_id_only`` option.
The latter option is not recommended for hybrid peeling or
any combination of different input files.

.. _peeling_methods:

Peeling methods
---------------

.. parsed-literal::

    Strategy:
      -method METHOD        Peeling method: single or multi. Default: multi.

    Single-locus options for the second stage of hybrid peeling:
      -seg_file SEG_FILE    :ref:`Segregation probabilities file <seg_prob_file_format>`.
      -seg_map_file SEG_MAP_FILE
                            Map file for loci in the segregation probabilities file
                             (:ref:`see format details <map_file_format>`).

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

  * a ``-map_file`` with positions for loci in the sequence allele read count data,
  * a ``-seg_file`` with segregation probabilities generated via multi-locus method, and
  * a ``-seg_map_file`` with genetic positions for loci in the segregation probabilities file.

This combination of options is not required in the standard multi-locus mode.

.. _peeling_parameters:

Peeling parameters
------------------

.. parsed-literal::

    Computational parameters:
      -n_cycle N_CYCLE      Number of peeling cycles.
                            Default: 5.
      -n_thread_fam N_THREAD_FAM
                            Number of threads to parallelise computation across families.
                            Default: 1.
      -n_thread_loci N_THREAD_LOCI
                            Number of threads to parallelise computation across loci.
                            Default: 1.

    Estimation of model parameters:
      -est_start_alt_allele_prob
                            Estimate from all inputted genomic data prior to peeling and 
                            output :ref:`alternative allele probabilities <alt_allele_prob_file_format>`.
      -est_alt_allele_prob  Estimate after each peeling cycle and output 
                            :ref:`alternative allele probabilities <alt_allele_prob_file_format>`.
      -est_geno_error_prob  Estimate after each peeling cycle and 
                            output :ref:`genotype error probabilities <geno_error_prob_file_format>`. 
      -est_seq_error_prob   Estimate after each peeling cycle and 
                            output :ref:`sequence error probabilities <seq_error_prob_file_format>`. 
      -est_pheno_penetrance_prob
                            Estimate after each peeling cycle 
                            and output :ref:`phenotype penetrance probabilities <pheno_penetrance_prob_file_format>`.
      -no_phase_founder     Suppress phasing a heterozygous allele 
                            (if such an allele can be found) in
                            genotyped individuals without genotyped parents.

The option ``-no_phase_founder`` can be used to suppress the default 
behaviour of phasing a midpoint heterozygote locus in a genotyped 
individual with ungenotyped parents,
where the alternative allele is set as paternally inherited. 
Calling this option, we expect an equal probability
of the alternative allele as paternally or maternally inherited.

Computational effort and speed of |Software| can be controlled with
the number of peeling cycles (``-n_cycle``,
increasing the number will marginally increase accuracy, but also runtime) and
the number of threads to parallelise computations across families and loci 
(``-n_thread_fam`` and ``-n_thread_loci``, to reduce runtime on large datasets).

Both multithreading options, ``-n_thread_fam`` and
``-n_thread_loci``, are used for parallelisation of the peeling process.
``-n_thread_fam`` controls the number of threads to parallelise accross
families in a generation, while ``-n_thread_loci`` controls the number of threads 
to parallelise across loci. Both options can be used together 
to speed up the analysis, but ``-n_thread_fam`` is more impactful than 
``-n_thread_loci`` when the number of loci is small (i.e., less than 3000
per chromosome) and computing resources are limited.
The ``-n_thread_loci`` option might be more useful if number of CPU cores 
is abundant and the number of loci is large 
(e.g., more than 5 cores and 4000 loci).

You can try the script described in :ref:`multi-threading benchmarking <multi-threading-benchmarking>`
to see how the number of threads affects the runtime of the 2000 loci simulation data on your device.

You are encouraged to find a combination that works best for your use case.

|Software| can estimate the model parameters from the input data.
The :ref:`default or user provided input values<input_options>`
are used as a starting point for estimation.
See a note on :ref:`robustness of results <robust_parameters>`
to these parameters.

When estimation options are used,
the respective parameters are estimated after each peeling cycle and
output to a file at the end.
This process usually increases running times and
might require additional peeling cycles to converge.
The estimates are based on inferred states of the modelled events and
their match/mismatch between observed and inferred states.

Alternative allele probabilities are estimated (using ``-est_alt_allele_prob``)
as half of the mean of estimated :ref:`allele dosage <dosage_file_format>`
in the base population(s) (metafounders).

This estimation can be warm-started with a sample estimate from inputted genomic data
(using ``-est_start_alt_allele_prob``).
Note that this sample estimate is not taking the pedigree structure into account,
so it is a naive estimate and does not pertain to the pedigree base population(s)
and hence ignores metafounders.
Option ``-est_start_alt_allele_prob`` uses Newton optimisation,
which also requires starting values.
These starting values are by default 0.5,
but can also be provided using ``-alt_allele_prob_file``.

For a pedigree with :ref:`multiple metafounders <ped_file_format>`,
there are three options to obtain metafounder-specific alternative allele probabilities:

  (1) use the default starting values of 0.5 for all loci and then
  ``-est_alt_allele_prob``,

  (2) use ``-est_start_alt_allele_prob`` to use sample-based starting values and
  then ``-est_alt_allele_prob``, or
  
  (3) provide starting values using ``-alt_allele_prob_file`` and
  then ``-est_alt_allele_prob``.
  
In all three cases,
``-est_alt_allele_prob`` is optional.

The estimates (using ``-est_start_alt_allele_prob`` or ``-est_alt_allele_prob``) or 
inputs (using ``-alt_allele_prob_file``) of the alternative allele probabilities 
are constrained to be between 0.001 and 0.999 to ensure valid probabilities and 
avoid getting trapped in boundary values, 0 or 1.

Error probabilities (using ``-est_geno_error_prob`` and ``-est_seq_error_prob``)
are estimated as the proportion of mismatches between observed and inferred states.

Phenotype penetrance probabilities (using ``-est_pheno_penetrance_prob``) are estimated
from the average conditional probabilities of phenotype states for each phased genotype
across phenotyped individuals. This follows Kinghorn (2003)
A simple method to detect a single gene that determines a categorical trait with incomplete penetrance.
Assoc. Advmt. Anim. Breed. Genet. 15:103-106.
