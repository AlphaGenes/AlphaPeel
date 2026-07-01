.. _seg_prob_file_format:

Segregation probability
=======================

Format
~~~~~~

There are four lines per individual,
corresponding to the four possible patterns of segregation (inheritance):

  1. The grand *paternal* allele from the father and
     the grand *paternal* allele from the mother: :math:`Pr(F_{p,-} \& M_{p,-})`,
  2. The grand *paternal* allele from the father and
     the grand *maternal* allele from the mother: :math:`Pr(F_{p,-} \& M_{-,m})`,
  3. The grand *maternal* allele from the father and
     the grand *paternal* allele from the mother: :math:`Pr(F_{-,m} \& M_{p,-})`, and
  4. The grand *maternal* allele from the father and
     the grand *maternal* allele from the mother: :math:`Pr(F_{-,m} \& M_{-,m})`.

The symbols above mean the following:
``Pr()`` probability,
``F`` father,
``M`` mother,
``p`` *paternal* allele,
``m`` *maternal* allele, and
``-`` the allele that was not inherited.

.. TODO: Do we need a diagram for this maybe even show the 4 cases so we bring home the message!?
.. https://github.com/AlphaGenes/AlphaPeel/issues/227

The first value in each line is the individual ID.
The remaining values are *segregation probabilities* at each locus.

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


Input details
~~~~~~~~~~~~~

This input file is obtained from the output of the multi-peeling cycle of the hybrid peeling,
and is only required and allowed for the single-peeling cycle of the hybrid peeling.
If provided, values from this file initializes the *segregation probabilities* in the current analysis. 
See details of ``segregationTensor`` in :ref:`function explanation <function_explanation>` for more explanation of 
how segregation probabilities are used in the peeling cycle. 


Output details
~~~~~~~~~~~~~~

The file is saved as ``.seg_prob.txt``.
The file contains all individuals from all inputs.
