=========
Changelog
=========

[1.3.0] - 2026-02-05
====================

New features
------------

* Allow phenotype input, currently only works for one phenotype and one genotype/locus,
  currently still experimental (:pr:`195`, :user:`RosCraddock`, :user:`XingerTang`).

    - Added ``pheno_file`` file with phenotypes for phenotyped individuals.

    - Added ``pheno_penetrance_file`` file with user-provided phenotype error rates
      for the phenotypes, that is, the conditional probability of each phenotype
      given the genotype. TODO: Rename to `pheno_error_prob_file`

    - Added ``pheno_prob`` to output the phenotype probabilities of individuals.

    - Added ``update_pheno_penetrance_file`` to re-estimate the phenotype penetrance
      after each peeling cycle following Kinghorn (2003).

* Allow X chromosome input,
  currently still experimental (:pr:`198`, :user:`AprilYuZhang`, :user:`XingerTang`).

    - Move X chromosome flag ``sex_chrom`` to ``x_chr``.

    - Update the X chromosome related peeling functions.

* Add mutation rate input (:pr:`198`, :user:`AprilYuZhang`, :user:`XingerTang`).

    - Added ``mutation_rate`` to allow user-provided mutation rate.

Add map file input (:pr:`208`, :user:`XingerTang`, :user:`gregorgorjanc`).

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

* Updates to accuracy and functional tests for new argument names
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

* Added instructions on how to contribute to AlphaPeel
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
