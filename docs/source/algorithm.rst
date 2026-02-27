==================
Algorithm Overview
==================

Here are some notes about how the peeling algorithm works.
For further details read:

.. [1] Whalen, A, Ros-Freixedes, R, Wilson, DL, Gorjanc, G, Hickey, JM. (2018). *Hybrid peeling for fast and accurate calling, phasing, and imputation with sequence data of any coverage in pedigrees*. Genetics Selection Evolution; doi: https://doi.org/10.1186/s12711-018-0438-2

Peeling basics
==============

The idea behind multi-locus iterative peeling is
to estimate an individual's genotype based on information from their own SNP gentoype data,
the genotypes of their parents, and the genotypes of their offspring.
This can be challenging because we have to combine these different sources of information and
there can be feedback loops that need to be managed when doing this.
For example, if a child's genotypes were imputed based on the genotypes of their parents,
we would not want to use those genotypes directly to re-impute the parent.
Instead we would want to use only the "new" information from the genotypes of the child
(and potentially their offspring).

To keep these sources of information separate,
we explicitly calculate and store three different values.
These are,

- anterior: genotype probabilities for the individual based on their parent's genotypes,
- penetrance: genotype probabilities for an individual based on their own genotypes, and
- posterior: genotype probabilities for an individual based on the offspring's genotypes.
  For the penetrance term you will often need the mate's genotypes as well.

For each of these values we store genotype probabilities over the four, phased genotype states,
that is, ``aa``, ``aA``, ``Aa``, and ``AA``.
In all cases the first allele is the allele inherited from the father,
and the second is the allele that was inherited from the mother.

One of the key things about inheritance and associated probability calculations
is that individual's inherit their chromosomes from their parents in large blocks.
If we knew which blocks an individual inherited,
we could use that information to help us determine what alleles the individual carried,
and also to determine which alleles their parents carried.
For each individual, we estimate these values using a segregation probabilities.

We code segregation probabilities in two forms,

- as the joint probability for both the paternal and maternal segregation
  (given by four values, ordered ``pp``, ``pm``, ``mp``, and ``mm``, where the ``p`` denotes the paternal allele,
  and the ``m`` denotes the maternal allele, while position is as before - the first allele inherited from father,
  meaning that ``pp`` denotes that individual inherited paternal allele from the father and mother),
- or as the probability for either the father or the mother
  (a pair of single values giving the probability that the individual inherits the maternal haplotype,
  that is, a seg=1 means that the individual inherits the maternal haplotype for that parent).

In a lot of places, we will rely heavily on calculating the "inheritance" or "transmission" probabilities.
These give the probability that a child inherits an allele based on the genotypes of their parents.
For some cases this is simple, for example, if the father is ``AA``, and the mother is ``aa``
then their offspring will be ``Aa`` since we know that the father will  transmit a copy of ``A`` and
the mother will transmit a copy of ``a`` so the resulting offspring will be ``a``.

Other cases may be more challenging and involve multiple possible outcomes
with associated probabilities.
If the father is heterozygous, ``aA``, and the segregation is unknown,
then there will be
a 50% probability that the offspring will inherit an ``a`` and
a 50% probability they inherit a ``A``.
If we take such a father and the mother is ``AA``
then the genotype probabilities of the offspring are:

::

    p(aa): 0.0
    p(aA): 0.5
    p(Aa): 0.0
    p(AA): 0.5

However we may know which haplotype the individual inherits.
If the segregation of the individual is ``mm``
(it inherited the maternal allele from both parents)
then the resulting genotype probabilities of an offspring
from an ``aA`` father and an ``AA`` mother are:

::

    p(aa): 0.0
    p(aA): 0.0
    p(Aa): 0.0
    p(AA): 1.0

Segregation probabilities are really helpful for determining which alleles an individual inherits
from their parents and are used for both the peel down (anterior) and peel up (posterior) steps.
The following sections outline how the penetrance, anterior, posterior, and segregation terms are calculated.

Penetrance
----------

The penetrance term gives the probability of each of the four phased genotypes
based on the direct information we have about an individual's genotype.
It is calculated via the helper function ``tinyhouse.ProbMath.getGenotypeProbabilities_ind``.

The basic idea is that we want to take a genotype and turn it into genotype probabilities.

- If the individual is genotyped as a ``0``, that is ``aa``, we want to set ``p(aa)=1.0``.
- If the individual is genotyped ``1``, that is ``aA`` or ``Aa``, we want to set ``p(aA) = p(Aa) = 0.5``.
- If the individual is genotyped as a ``2``, that is ``AA``, we want to set ``p(AA) = 1.0``.

The difference between the heterozygous and homozygous states is that
with the heterozygous state there are two possible allele combinations that produce the same genotype.
In all of the cases we also want to add a small amount of noise to the estimate to allow for genotyping errors.
This error is often called as "penetrance" error.

Anterior
--------

The anterior term is calculated differently for the founders and non-founders.

The anterior term for founders in the population (individuals without any parents)
are set using the population minor allele frequency.
To set the term, we use:

::

    pedigree.setMaf()
    founder_anterior = ProbMath.getGenotypesFromMaf(pedigree.maf)
    founder_anterior = founder_anterior*(1-0.1) + 0.1/4
    for ind in pedigree:
        if ind.isFounder():
            ind.peeling_view.setAnterior(founder_anterior.copy())

For non-founders the anterior value is based on the genotypes of their parents.
The probability that the inherit an ``A`` allele from their father if their segregation is unknown is:

::

    P(a|father) = 0.5 p(aA) + 0.5 p(Aa) + 1 p(AA)

where ``p(aA)`` relates to the genotype probabilities of their father.
If the segregation value is known, this is instead:

::

    P(a|father) = seg p(aA) + (1-seg) p(Aa) + 1 p(AA)

We can generate genotype probabilities for the four phased genotype
as the product of the allele probabilities for each parent:

::

    p(aA) = p(a|father) p(A|mother)

Posterior
---------

The posterior term is a bit trickier than the anterior term.
The idea is that for each full-sub family (that is, each father-mother mate pair)
we will use their children to estimate their genotype  probabilities.
To do this, we construct the join probabilities of their parent's genotypes,
a 4x4 matrix of values.
We then calculate the posterior estimate for a single parent by
marginalizing the joint genotype probabilities by the genotypes of the other parents.

The joint genotypes are estimated by

.. math::
    p(g_{father}, g_{mother}|children) = \prod(p(g_{father}, g_{mother}|child)).

Simulation results using AlphaPeel have suggested that accuracy may be increased
by using the called genotype probabilities.
Because of this we call the individual's genotypes, haplotypes, and segregation values.
This has the added benefit of allowing us to use a look-up table to
produce the joint genotype probabilities at a given locus.

We then marginalize the genotypes for each parent by,

.. math::
    p(g_{father}|children) = \sum(p(g_{father}, g_{mother}|children)p(g_{mother})).

We assume that the posterior values for each family are independent.
This lets them calculate them separately for each family group, and
then sum them together to produce the final called genotype.
Because some individuals have a large number of offspring,
we calculate the joint genotypes on a log scale and then convert back in the marginalization step.

Segregation
-----------

We use a standard forward-backward algorithm.
The transmission rate determines the uncertainty when moving from one loci to the next.

Specific code comments
======================

General math
------------

**Exponential + normalization:**
When working with the posterior values it is important to handle overflow/underflow.
We do this by treating most of the posterior values on a log scale.
To convert back from a log scale to a normal scale requires an exponential operation.
Most of the cases where we do this, we also need to normalize the resulting values so that they sum to one.
In order to prevent issues with underflows,
we first calculate the maximum value in the slice of the array that we are dealing with.
We then subtract the maximum value and take the exponential.
This means that the greatest value in the array will be set to :math:`exp(0)=1`.
Very small values may be set to zero, but this just indicates that they have vanishingly small probability.

Parallel
--------

In order to parallelize the code, we take this approach.
The overall idea is to find blocks of individuals
who can be updated at the same time (in parallel) and perform those updates.
In the context of peeling, we can break up individuals into discrete generations and
perform updates on all of the individuals in the same generation at the same time.
We use a ``ThreadPoolExecutor`` to perform tasks in parallel, and
use ``numba``'s ``nogil`` flag to get around python's global interpreter lock.

Because of the overhead in crating threads, and calling ``numba`` functions,
we split out the tasks in groups of full-sub families,
which will be updated at the same time.
Because a given parent can be a parent of multiple families
(but an offspring can only be an offspring of one family),
we set the genotypes for the parents separately in a non-parallel mode.
A complete breakdown of the main steps (and parallelization) is:

- **father and mother genotypes** (no parallel):
  We set the genotypes of the fathers and mothers in a non-parallel mode.
  Individual father and mothers may be used in multiple families within the same generation,
  and so setting them within a family block may lead to conflicts.
- **Child genotypes** (parallel):
  We set the genotypes of the children in parallel.
  This is safe because children are only members of a single full sib family
  (the family which contains their father and mother).
  We also collapse the posterior terms for individuals with offspring in this step
  when estimating the posterior term for their parents.
  This is done because all of the posterior terms for the individual
  should have been calculated by this point
  (the offspring's generation is always greater than their parent's).
- **Child segregation** (parallel):
  Individual child segregation values are re-estimated in parallel.
  These values only depend on the genotypes of the child and their parents.
  Their parents genotypes are fixed for the parallel blocks, and
  their genotypes are set on a family-by-family basis.
- **Child anterior term** (parallel):
  The anterior terms are estimated in parallel.
  These depend only on the parent genotype probabilities which fixed for a given generation.
- **Parent posterior** (parallel, but outside of ``numba``):
  We add the family posterior elements to each of the parents in parallel but outside of ``numba``.
  I think this should be safe because the code is run in multi-threaded mode,
  and so the GIL should make the append operations threadsafe.

Speed
-----

- Time seems to be split pretty evenly between the anterior and posterior computation terms.
- On the whole, the posterior term seems to dominate compared to the anterior term (usually by a factor of ~2).
- Estimating the segregation appears to be fairly low cost.
- There do not seem to be any obvious speed bottlenecks.

Memory
------

The storage of values takes place inside ``jit_Peeling_Individual`` in ``ImputationIndividual``. The main costs are:

- three ``4xnLoci float32`` s:
    - anterior
    - penetrance
    - genotypeProbabilities
- For individuals with offspring there are an additional two ``4xnloci float32`` s:
    - posterior
    - newPosterior
- one ``2xnLoci float32``:
    - segregation

There are a lot of possible places to obtain substantial memory savings.

- For all individuals
    - Because the anterior terms segregate independently, this could be reduced down to ``2xnLoci float32``,
      and re-calculated on the fly.
    - All of the ``4xnLoci float32`` s could potentially be stored as ``int8`` s with some conversion from ``int8->float32``.
      We don't actually need these terms to be very accurate (as long as we can still accurately call values).
- Individuals without offspring
    - For individuals without offspring we only ever use their called genotypes
      (with just the penetrance term) and
      segregation estimates in the peeling.
    - These terms come to play in the context of the posterior term for their parents.
    - This means we could potentially just save their genotypes, haplotypes, and segregation estimate.
    - If we need to call the individual in the future, we can re-calculate their anterior term on the fly, and
      recombine with the individual's genotype.
- For individuals with offspring:
    - We only ever use the penetrance+posterior or penetrance+anterior+posterior.
      We could calculate and store these values explicitly instead of storing them independently and
      re-calculating the genotype probabilities.
    - We currently store the posterior estimates as a list and re-add.
      We could instead store the values as a single matrix and just add each time.
      We need to be careful with the parallel updates on this term though.

Function explanation
====================

The main peeling function of AlphaPeel is given by ``tinypeel.Peeling.Peeling.peel()`` function:

.. autofunction:: tinypeel.Peeling.Peeling.peel

The peeling process consists of two parts:

1. The first part performs a Baum-Welch-like algorithm over each family of each generation of the pedigree.
    The details of the HMM are the following:

   - Hidden states: phased genotype.
   - Time dimension: generation (parent -> child).
   - Observed states: input genotype and sequence data.
   - Emission function: combined with the observed states, are introduced via the ``penetrance`` of the peeling information container,
     which is generated from ``tinypeel.tinyhouse.ProbMath.getGenotypeProbabilities()``.

     - For genotype input: suppose genotype error is :math:`e`, the following table represents
       how the input genotype data is encoded to be the probabilities of the phased genotype data,
       with rows representing the phased genotypes and columns representing the input genotype data:

       .. list-table::
          :header-rows: 1

          * -
            - aa
            - aA
            - Aa
            - AA
          * - 0
            - 1 - e
            - e/2
            - e/2
            - e/2
          * - 1
            - e/2
            - 1 - e
            - 1 - e
            - e/2
          * - 2
            - e/2
            - e/2
            - e/2
            - 1 - e

     - For sequence input: suppose sequence error is :math:`e`, the following table represents
       how the input sequence data is encoded to be the phased genotype data,
       with :math:`r` representing the number of reference alleles
       and :math:`a` representing the number of alternative alleles:

       .. list-table::
          :header-rows: 1

          * - phased genotype
            - corresponding probability
          * - aa
            - :math:`(1 - e)^r e^a`
          * - aA
            - :math:`\left(\frac{1}{2}\right)^{r + a}`
          * - Aa
            - :math:`\left(\frac{1}{2}\right)^{r + a}`
          * - AA
            - :math:`e^r (1 - e)^a`

      All the values are normalised before use.
      In the case that both sequence and genotype inputs are provided,
      the two probabilities are multiplied together.

   - Transmission function: implemented via ``segregationTensor`` defined in ``tinypeel.tinyhouse.ProbMath()``.
     The ``segregationTensor`` is a 4D numpy array of float32 of size :math:`4 \times 4 \times 4 \times 4`, with

     * the first dimension represents the paternal genotype (``p``),

     * the second dimension represents the maternal genotype (``m``),

     * the third dimension represents the child's genotype (``allele``), and

     * the last dimension represents the child's segregation (``seg``).

     Mathematically, ``segregationTensor`` represents the probability of each combination of
     the father's genotype, the mother's genotype,
     child's genotype and segregation without any other information (:math:`P(p, m, allele, seg)`).
     The following are two examples:

     * Example 1: Suppose ``allele = 0`` (genotype = aa) and ``seg = 0`` (segregation = pp),
       the values of ``segregationTensor[:, :, 0, 0]`` are:

     .. list-table::
          :header-rows: 1

          * -
            - m = aa
            - m = aA
            - m = Aa
            - m = AA
          * - p = aa
            - 1
            - 1
            - 0
            - 0
          * - p = aA
            - 1
            - 1
            - 0
            - 0
          * - p = Aa
            - 0
            - 0
            - 0
            - 0
          * - p = AA
            - 0
            - 0
            - 0
            - 0

     * Example 2: Suppose ``allele = 2`` (genotype = Aa) and ``seg = 1`` (segregation = pm),
       the values of ``segregationTensor[:, :, 2, 1]`` are:

     .. list-table::
          :header-rows: 1

          * -
            - m = aa
            - m = aA
            - m = Aa
            - m = AA
          * - p = aa
            - 0
            - 0
            - 0
            - 0
          * - p = aA
            - 0
            - 0
            - 0
            - 0
          * - p = Aa
            - 1
            - 0
            - 1
            - 0
          * - p = AA
            - 1
            - 0
            - 1
            - 0

     A pre-defined error :math:`e` is used here with ``segregationTensor = segregationTensor*(1-e) + e/4``
     (mutation error).

     This matrix can be used to generate ``ChildSegs``, which is the matrix controls
     how the information is passed across generations.
     By ``tinypeel.Peeling.Peeling.createChildSegs()``, the probabilities of each combination of
     father's genotype, mother's genotype and child's genotype can be calculated via
     summing over the child's segregation states.

     .. autofunction:: tinypeel.Peeling.Peeling.createChildSegs

     The information are passed with functions ``tinypeel.Peeling.Peeling.projectChildGenotypes()`` and
     ``tinypeel.Peeling.Peeling.projectParentGenotypes()``.

     .. autofunction:: tinypeel.Peeling.Peeling.projectChildGenotypes

     .. autofunction:: tinypeel.Peeling.Peeling.projectParentGenotypes

     The pre-defined error :math:`e` is used in a similar way as the ``segregationTensor`` on the following variables:

     - ``JointParents``,
     - ``probSire`` and ``probDam``,
     - ``childValues``, and
     - ``sirePosterior`` and ``damPosterior``,

     which are all intermediate values or the results of the transmission across generation.

   - Probabilities update: The phased genotype probabilities are calculatd via ``anterior * penetrance * posterior``,
     where:

     * ``anterior``: is updated when peeling down, every time a new value is calculated.
     * ``penetrance``: depends only on input.
     * ``posterior``: is updated when peeling up, but only when a peeling cycle is finished.

2. The second part performs a Baum-Welch-like algorithm over each locus of each individual.
   This part is implemented only when the multi-locus peeling mode is used.
   The details of the HMM is the following:

   - Hidden states: segregation states.
   - Time dimension: locus (locus_i -> locus_i+1).
   - Observed states: phased genotype probabilities from part 1.
   - Emission function: A segregation estimate ``pointSeg`` is obtained by ``estimateSegregationWithNorm()``,

     .. autofunction:: tinypeel.Peeling.Peeling.estimateSegregationWithNorm

     then same usage of :math:`e` as the ``segregationTensor`` is applied,
     but now on the resulting output of function for the Baum-Welch algorithm implementation:

     * matched segregation: :math:`1 - \frac{3}{4}e`,
     * unmatched segregation: :math:`\frac{1}{4}e`.

   - Transmission function: first generate equal recombination rates ``r`` across all loci,
     which the values are calculated by assuming there is exactly 1 recombination happened per snippet
     of input and the distances between each locus are equal.
     The transmission function is setting up as the following:

     * if the segregation states are same at locus i and locus i + 1: :math:`(1 - r)^2`,
     * if the segregation states are different by one parent at locus i and locus i + 1: :math:`r \times (1 - r)`,
     * if the segregation states are different by both parents at locus i and locus i + 1: :math:`r^2`.
