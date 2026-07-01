.. _pheno_penetrance_prob_file_format:

Phenotype penetrance probability
================================

Format
~~~~~~

There are four lines per locus,
containing *conditional probabilities of
all phenotypes for a given phased genotype*,
including phenotyping errors or other deviations.
These probabilities encode the relationship between genotypes and phenotypes.
Each column represents a different phenotype state
(ordered from ``0`` to the total number of states), and
each row represents the :ref:`phased genotypes <zero_one_two_etc>`
(``aa``, ``aA``, ``Aa``, and ``AA``).
That is, the values within a row
form a conditional probability distribution
of phenotypes for the corresponding phased genotype.
All phenotypes present in phenotype files should be included in the file.
Current functionality works only with one locus and one trait.

Example for a monogenic recessive binary trait with 10% phenotyping error:

::

  0.9 0.1
  0.9 0.1
  0.9 0.1
  0.1 0.9


Input details
~~~~~~~~~~~~~

Values from this file parameterise how observed phenotypes  
translate to individual's :ref:`penetrance <penetrance>` 
component of genotype probabilities.

Output details
~~~~~~~~~~~~~~

The file is saved as ``.pheno_penetrance_prob.txt``.
