.. _ped_file_format:

Pedigree
========

Format
~~~~~~

This file has one line with
*recorded parents* for each individual.
Each line has three values,
the individual's ID,
their father's ID, and
their mother's ID.
Unknown parent ID is encoded with ``0``.

Example with four individuals across two generations:

::

  id1 0 0
  id2 0 0
  id3 id1 id2
  id4 id1 id2

Users can define multiple base populations (metafounders) as shown below.
The IDs of these metafounders must start with ``MF_``.

Example with two metafounders:

::

  id1 MF_1 MF_1
  id2 MF_2 MF_2
  id3 id1 id2
  id4 id1 id2

When working with the X chromosome,
a fourth value is needed for individual's sex:
``0`` for males and ``1`` for females.

Example with sex information:

id1 and id3 are males, while id2 and id4 are females:

::

  id1 0 0 0
  id2 0 0 1
  id3 id1 id2 0
  id4 id1 id2 1

Input details
~~~~~~~~~~~~~

The file(s) should include all individuals present in other files
and their known relatives present only in pedigree.
Individuals with one unknown parent are internally assigned a "dummy" parent.
Therefore, all individuals have either both parents known or neither parent known.
Individuals with two unknown parents are considered as founders and
as such represent the base population of the pedigree.
Internally, they are assigned to a metafounder (unknown parent group) ``MF_1``
(or whatever is specified through ``-main_metafounder``).
