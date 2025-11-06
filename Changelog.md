
# Changelog


## [1.2.0] - 2025-03-21
**New features**
* allow metafounders, defined as “MF_”, in the pedigree file input ({pr}`175`, {user}`RosCraddock`, {user}`XingerTang`)

* added ``alt_allele_prob_file`` to allow user-inputted alternative allele frequencies for each metafounder and loci. For now, these are restricted to be between 0.01 and 0.99 ({pr}`175`, {user}`RosCraddock`, {user}`XingerTang`)

* added ``main_metafounder`` to allow user to assign the default metafounder to use where a metafounder has not been assigned to a founder in the pedigree ({pr}`175`, {user}`RosCraddock`, {user}`XingerTang`)

* added ``update_alt_allele_prob`` to allow the base alternative allele frequencies to be updated after each peeling cycle based on the mean of the founders within the assigned metafounder ({pr}`182`, {user}`RosCraddock`)

**Bug fixes**
* Fixed bug due to setuptools package being updated for all wheel file naming to follow binary distribution specification (i.e all lower case) and updated documentation ({pr}`182`, {user}`RosCraddock`)

**Maintenance**
* User-warnings and documentation updates for metafounder implementation and estimation of alternative allele frequency ({pr} `152`, {pr}`175`, {pr}`182`, {user}`RosCraddock`, {user}`XingerTang`)

* Functional and accuracy tests for metafounder implementation ({pr}`156`, {pr}`182`, {user}`XingerTang`, {user}`RosCraddock`)

* Updated reference to tinyhouse ({pr}`177`, {user}`XingerTang`)


## [1.1.6] - 2024-10-22
**New features**
* Addition of map file input for non-hybrid mode ({pr}`154`, {user}`XingerTang`)

**Bug fixes**
* resolved bug to produce output file with ``-hap`` and ``-geno`` ({pr}`157`, {user}`AprilYUZhang`)

**Maintenance**
* set default hap and geno threshold as 1/3 when calling genotypes ({pr}`157`, {user}`AprilYUZhang`)

## [1.1.5] - 2023-12-01
**New features**
* Addition of output options: ``geno``, and ``hap_threshold`` ({pr}`119`, {user}``AprilYuZhang`)

**Maintenance**
* Updates in option and file names. These include:

  * ``no_dosages`` to ``no_dosage``
  * ``calling_threshold`` to ``geno_threshold``
  * ``call_phase`` to ``hap``
  * ``haps`` to ``phased_geno_prob``
  * ``pedigree`` to ``ped_file``
  * ``genotypes`` to ``geno_file``
  * and more, for all changes please visit: https://github.com/AlphaGenes/AlphaPeel/issues/113#issue-1935197000 ({pr}`105`, {pr}`115`, {pr}`122`, {user}`XingerTang`, {user}`AprilYuZhang`)

* Updates the documentation and help functions ({pr}`88`, {pr}`119`, {user}`XingerTang`, {user}`AprilYuZhang`).

* Updates to accuracy and functional tests for new argument names ({pr}`126`, {pr}`130`, {pr}`131`, {user}`XingerTang`).

## [1.1.4] - 2023-08-25
**New features**
* Implementation of functional and accuracy testing with pytest ({pr}`53`, {user}`XingerTang`)

* Implementation of pre-committ protocol through Black and Flake8 ({pr}`40`, {user}`XingerTang`)

* Implementation of cross-platform tests workflow through github actions ({pr}`59`, {user}`XingerTang`)

* Added instructions for how to contribute to AlphaPeel ({pr}`27`, {pr}`31`, {user}`XingerTang`)

**Bug Fixes**
* Fixed bug on the loading of submodule ({pr}`22`, {user}`XingerTang`)

**Maintenance**
* Update theme for the HTML documentation ({pr}`38`, {user}`XingerTang`)

* Modified the URL for installation ({pr}`11`, {user}`XingerTang`)
