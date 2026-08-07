=========
Changelog
=========

.. 
  If you are add a new entry to the changelog, please follow the format of the existing entries and 
  add your changes to the top of the file under the [Unreleased] section. 
  Please also include the pull request number and the GitHub username(s) of the contributor(s) who made the change.

[1.4.0] - 2026-08-07
====================

New features
------------
* Allow a founding individual to have two different metafounders, e.g., one for the paternal lineage and one for the maternal lineage
  (:pr:`222`, :user:`RosCraddock`, :user:`gregorgorjanc`, :user:`XingerTang`).

* Add a version option to the command line interface
  (:pr:`237`, :user:`XingerTang`).

* Add support for direct in-process AlphaPeel calls for accuracy benchmarking,
  profiling, and development workflows (:pr:`292`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Add locus-based multi-threading options and benchmarking support for comparing
  family-wise and locus-wise thread settings (:pr:`292`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Add switch-error rate and phasing-error rate calculation utilities for accuracy
  reporting (:issue:`255`, :pr:`256`, :user:`XingerTang`).

* Add controls for output rounding precision (:issue:`271`, :pr:`272`, :user:`XingerTang`).

* Add explicit output options for estimated genotype error probabilities,
  estimated sequence error probabilities, and recombination probabilities
  (:issue:`252`, :pr:`268`, :user:`XingerTang`).


Bug fixes
---------
* Fix bug for subsetting snps with map file input (:pr:`241`, :user:`XingerTang`).

* Fix X chromosome peeling behaviour, including penetrance handling and default
  segregation probabilities (:issue:`225`, :pr:`226`, :user:`AprilYUZhang`).

* Fix the accuracy report test configuration file names
  (:issue:`259`, :pr:`262`, :user:`XingerTang`).

* Fix estimation of alternative allele probabilities when using up to two
  metafounders per individual (:issue:`260`, :pr:`263`, :user:`RosCraddock`).

Maintenance
-----------
* Updated the documentation for ``alt_allele_prob_file`` and estimation of
  phenotype penetrance probabilities (:pr:`220`, :user:`RosCraddock`).

* Update the alternative allele probability limits from [0.01, 0.99] to [0.001, 0.999]
  (:issue:`185`, :issue:`276`, :pr:`277`, :user:`RosCraddock`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Renaming of command line arguments, changing of input format, and corresponding documentation and test updates
  (:issue:`221`, :pr:`219` ,:pr:`222`, :user:`RosCraddock`, :user:`gregorgorjanc`, :user:`XingerTang`).

    - ``est_alt_allele_prob`` to ``est_start_alt_allele_prob``.

    - ``update_alt_allele_prob`` to ``est_alt_allele_prob``.

    - ``pheno_penetrance_file`` to ``phenotype_penetrance_prob_file``.

    - ``update_pheno_penetrance`` to ``est_pheno_penetrance_file``.

    - Reformatted the ``alt_allele_prob_file`` input to match the outputted ``alt_allele_prob_file``.

* Updated the documentation to clarify inputs, outputs, and parameters
  (:pr:`219`, :user:`gregorgorjanc`).

* Add collaboration guidelines to the documentation
  (:pr:`224`, :user:`XingerTang`, :user:`RosCraddock`, :user:`AprilYuZhang`, :user:`gregorgorjanc`).

* Update the ``tinyhouse`` submodule reference after metafounder updates
  (:issue:`228`, :pr:`229`, :user:`RosCraddock`).

* Remove obsolete repository files and make minor documentation and packaging
  refinements (:issue:`233`, :pr:`234`, :user:`XingerTang`).

* Updated the documentation to add more techincal instructions for developers and contributors
  (:pr:`237`, :user:`XingerTang`).

* Round up threshold for output file names (:pr:`246`, :user:`XingerTang`).

* Update help messages for command line arguments (:pr:`246`, :user:`XingerTang`).

* Rename option ``mutation_rate`` to ``mut_prob`` for consistency with other options
  (:pr:`246`, :user:`XingerTang`).

* Add explanation of ``no_phase_founder`` option in the documentation
  (:pr:`246`, :user:`XingerTang`).

* Add explanation of ``rec_length`` option in the documentation
  (:pr:`246`, :user:`XingerTang`).

* Add explanation of file format for the ``seg_map_file`` in the documentation
  (:pr:`246`, :user:`XingerTang`).

* Remove extra empty lines in the genotype probabilities output file
  (:pr:`246`, :user:`XingerTang`).

* Remove the ineffective IO multi-threading option
  (:issue:`210`, :issue:`250`, :pr:`251`, :user:`XingerTang`).

* Remove the sequence file requirement from the accuracy test runner
  (:issue:`257`, :pr:`261`, :user:`AprilYUZhang`).

* Refactor warning handling for clearer runtime messages
  (:issue:`266`, :pr:`267`, :user:`XingerTang`).

* Clarify documentation and help text for method selection, alternative allele
  probability output, and male X chromosome genotype probability formatting
  (:pr:`275`, :user:`XingerTang`).

* Clarify haplotype and penetrance input handling, remove the unsupported
  reference option, and update the ``tinyhouse`` submodule for renamed haplotype
  input messages (:issue:`73`, :issue:`74`, :issue:`127`, :issue:`167`,
  :issue:`231`, :pr:`278`, :user:`XingerTang`).

* Restructure the documentation to make file format documentation easier to navigate and cross-reference
  (:pr:`283`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Remove unsupported PLINK IO functionality and its functional tests
  (:issue:`286`, :pr:`291`, :user:`XingerTang`).

* Rework the accuracy test and reporting scheme to produce clearer benchmark
  outputs, generation-wise metrics, switch-error metrics, runtime summaries,
  and visualization-ready reports (:pr:`292`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Optimize memory use and runtime in the peeling implementation, including
  lower-memory data structures, streamlined array handling, and faster
  alternative allele probability updates (:pr:`292`, :user:`XingerTang`, 
  :user:`gregorgorjanc`).

* Refactor code and documentation for clearer module-level reuse
  (:pr:`292`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Add developer documentation for testing, benchmarking, memory profiling,
  runtime profiling, multi-thread profiling, and code-quality
  checks (:pr:`292`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Update continuous-integration configuration and supported Python-version
  testing setup for the current development workflow (:pr:`292`, :user:`XingerTang`,
  :user:`gregorgorjanc`).

[1.3.0] - 2026-02-05
====================

New features
------------

* Allow phenotype input, currently only works for one phenotype and one genotype/locus,
  currently still experimental (:pr:`195`, :user:`RosCraddock`, :user:`XingerTang`).

    - Added ``pheno_file`` file with phenotypes for phenotyped individuals.

    - Added ``pheno_penetrance_file`` file with user-provided phenotype error rates
      for the phenotypes, that is, the conditional probability of each phenotype
      given the genotype.

    - Added ``pheno_prob`` to output the phenotype probabilities of individuals.

    - Added ``update_pheno_penetrance_file`` to re-estimate the phenotype penetrance
      after each peeling cycle following Kinghorn (2003).

* Allow X chromosome input,
  currently still experimental (:pr:`198`, :user:`AprilYuZhang`, :user:`XingerTang`).

    - Move X chromosome flag ``sex_chrom`` to ``x_chr``.

    - Update the X chromosome related peeling functions.

* Add mutation rate input (:pr:`198`, :user:`AprilYuZhang`, :user:`XingerTang`).

    - Added ``mutation_rate`` to allow user-provided mutation rate.

* Add map file input (:pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`).

    - Modified ``map_file`` to enable map file input for non-hybrid mode.

Bug fixes
---------

* Fix minor bug in simulation code for accuracy tests
  (:issue:`181`, :pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Fix bug in the example code and the corresponding accuracy check code
  (:issue:`205`, :issue:`206`, :pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Fix tinyhouse bug of incorrectly classified founders
  (:issue:`118`, :pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`, :user:`augustusgrant`).

* Fix map file bug for non-ascending order input
  (:issue:`155`, :pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Fix minor typo in the ``.gitattributes`` file name
  (:issue:`30`, :pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`).

* Fix the bug that ignores the first locus while calculating the accuracy in the accuracy test
  (:pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`).

Maintenance
-----------

* Update the documentation and tests for the phenotype file input
  (:pr:`195`, :user:`RosCraddock`, :user:`XingerTang`).

* Add docstrings for main peeling functions
  (:pr:`195`, :user:`RosCraddock`, :user:`XingerTang`).

* Add changelog (:pr:`195`, :user:`RosCraddock`, :user:`XingerTang`).

* Update the documentation and tests for the X chromosome peeling
  (:pr:`198`, :pr:`201`, :user:`AprilYuZhang`, :user:`XingerTang`).

* Update the documentation for mutation rate
  (:pr:`198`, :user:`AprilYuZhang`, :user:`XingerTang`).

* Move changelog to documentation; add algorithm section, add a simple example, and
  update installation instructions in the documentation
  (:pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`).

[1.2.0] - 2025-03-21
====================

New features
------------

* Allow metafounders, defined as “MF\_”, in the pedigree file input
  (:pr:`175`, :user:`RosCraddock`, :user:`XingerTang`)>

    - Added ``alt_allele_prob_file`` for user-inputted alternative allele frequencies
      for each metafounder and loci. For now, these are restricted to be between 0.01 and 0.99.

    - Added ``main_metafounder`` to allow user to assign the default metafounder to use where
      a metafounder has not been assigned to a founder in the pedigree.

    - Added ``update_alt_allele_prob`` to allow the base alternative allele frequencies
      to be updated after each peeling cycle based on the mean of the founders within the assigned metafounder

Bug fixes
---------

* Fixed bug due to setuptools package being updated for all wheel file naming
  to follow binary distribution specification (i.e., all lower case) and updated documentation
  (:pr:`182`, :user:`RosCraddock`).

Maintenance
-----------

* User-warnings and documentation updates for metafounder implementation and
  estimation of alternative allele frequency
  (:pr:`152`, :pr:`175`, :pr:`182`, :user:`RosCraddock`, :user:`XingerTang`).

* Functional and accuracy tests for metafounder implementation
  (:pr:`156`, :pr:`182`, :user:`XingerTang`, :user:`RosCraddock`).

* Updated reference to tinyhouse
  (:pr:`177`, :user:`XingerTang`).

[1.1.6] - 2024-10-22
====================

New features
------------

* Addition of map file input for non-hybrid mode
  (:pr:`154`, :user:`XingerTang`).

Bug fixes
---------

* Resolved bug to produce output file with ``-hap`` and ``-geno``
  (:pr:`157`, :user:`AprilYUZhang`).

Maintenance
-----------

* Set default hap and geno threshold as 1/3 when calling genotypes
  (:pr:`157`, :user:`AprilYUZhang`).

[1.1.5] - 2023-12-01
====================

New features
------------

* Addition of output options: ``geno``, and ``hap_threshold``
  (:pr:`119`, :user:`AprilYuZhang`).

Maintenance
-----------

* Updates in option and file names
  (:pr:`105`, :pr:`115`, :pr:`122`, :user:`XingerTang`, :user:`AprilYuZhang`),
  the major ones include:

    - ``no_dosages`` to ``no_dosage``,

    - ``calling_threshold`` to ``geno_threshold``,

    - ``call_phase`` to ``hap``,

    - ``haps`` to ``phased_geno_prob``,

    - ``pedigree`` to ``ped_file``,

    - ``genotypes`` to ``geno_file``,

    - and more, for all changes please visit: https://github.com/AlphaGenes/AlphaPeel/issues/113#issue-1935197000.

* Updates the documentation and help functions
  (:pr:`88`, :pr:`119`, :user:`XingerTang`, :user:`AprilYuZhang`).

* Updates to accuracy and functional tests for new option names
  (:pr:`126`, :pr:`130`, :pr:`131`, :user:`XingerTang`).

[1.1.4] - 2023-08-25
====================

New features / additions
------------------------

* Implementation of functional and accuracy testing with pytest
  (:pr:`53`, :user:`XingerTang`).

* Implementation of pre-commit code formatting with Black and Flake8
  (:pr:`40`, :user:`XingerTang`).

* Implementation of cross-platform tests workflow with GitHub actions
  (:pr:`59`, :user:`XingerTang`).

* Added instructions on how to contribute to ``AlphaPeel``
  (:pr:`27`, :pr:`31`, :user:`XingerTang`).

Bug fixes
---------

* Fixed bug on the loading of submodule
  (:pr:`22`, :user:`XingerTang`).

Maintenance
-----------

* Update theme for the HTML documentation
  (:pr:`38`, :user:`XingerTang`).

* Modified the URL for installation
  (:pr:`11`, :user:`XingerTang`).
