.. _map_file_format:

Map
===

Format
~~~~~~

This file has one line with
*genome information* for each :ref:`locus <markers>`.
Each line has three values, 
the chromosome number,
the locus name, and
the physical position (in base-pairs).

Example with four SNP loci on chromosome 1:

::

  1 snp_a 12483939
  1 snp_b 192152913
  1 snp_c 65429279
  1 snp_d 107421759

Input details
~~~~~~~~~~~~~

If provided, values from this file determines the 
*recombination rate* between loci.
The file should include all loci present in other files.
Only loci on one chromosome should be provided!

The input file for hybrid method, ``seg_map_file``, 
has the same format as the above.
It should include only the subsets of the loci that 
have been chosen to run the single peeling stage of 
the hybrid peeling.
