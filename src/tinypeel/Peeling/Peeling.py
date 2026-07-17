from numba import jit, float32
import numpy as np


# Defining variables for peel up and peel down. Ideally these would be characters, but numba does not support characters.
PEEL_UP = 0
PEEL_DOWN = 1


@jit(
    nopython=True,
    nogil=True,
)
def peel(family, operation, peelingInfo, singleLocusMode):
    """Compatibility dispatcher for the direction-specific peeling kernels."""

    if operation == PEEL_DOWN:
        peelDown(family, peelingInfo, singleLocusMode)
    else:
        peelUp(family, peelingInfo, singleLocusMode)


# This is the main peeling down function.
@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e16": float32, "e1e": float32},
)
def peelDown(family, peelingInfo, singleLocusMode):
    """Peel information down from parents to offspring.

    :param family: The family object that the peeling is performed on.
    :type family: class:`tinyhouse.Pedigree.Family`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :param singleLocusMode: A flag to indicate the mode of peeling.
        `False` if using multi-locus peeling, and `True` if using single-locus peeling.
    :type singleLocusMode: bool
    :return: None. The function modifies the peelingInfo object in place.
    """
    isXChr = peelingInfo.isXChr

    e = 0.000001
    e1e = 1 - e
    e4 = e / 4
    e16 = e / 16

    # Setup local variables from the peeling information container.
    anterior = peelingInfo.anterior
    penetrance = peelingInfo.penetrance
    posterior = peelingInfo.posterior
    segregation = peelingInfo.segregation

    segregationTensor = peelingInfo.segregationTensor
    segregationTensor_norm = peelingInfo.segregationTensor_norm

    nLoci = peelingInfo.nLoci
    nOffspring = len(family.offspring)
    sire = family.sire
    dam = family.dam
    fam = family.idn

    # Creating variables here:
    # childSegs: The segregation estimates for a particular child (These are re-used? so need to be stored)
    # allToParents: The projection of each child onto the parental genotypes.
    # parentsMinustChild: The estimate of the parent's genotypes minus the contribution from a specific individual.

    childToParentsCurrent = np.full((4, 4, nLoci), 0, dtype=np.float32)
    allToParents = np.full((4, 4, nLoci), 0, dtype=np.float32)
    probSire = np.full((4, nLoci), 0, dtype=np.float32)
    probDam = np.full((4, nLoci), 0, dtype=np.float32)
    childValuesCurrent = np.full((4, nLoci), 0, dtype=np.float32)

    # Some local variables. currentSeg is the segregation estimate of a child (but may be modified).
    currentSeg = np.full((4, nLoci), 1, dtype=np.float32)
    forwardSeg = np.full((4, nLoci), 1, dtype=np.float32)

    needsSegregationUpdate = not singleLocusMode
    childSegTensor = np.full((nOffspring, 4, 4, 4, nLoci), 0, dtype=np.float32)
    parentsMinusChild = np.full((nOffspring, 4, 4, nLoci), 1, dtype=np.float32)
    if needsSegregationUpdate:
        childValuesTensor = np.full((nOffspring, 4, nLoci), 0, dtype=np.float32)

    # Construct the joint parent genotypes based on the parent's anterior, penetrance, and posterior terms minus this family.
    setupParentGenotypeProbs(
        anterior,
        penetrance,
        posterior,
        peelingInfo.posteriorSireContribution[fam, :, :],
        peelingInfo.posteriorDamContribution[fam, :, :],
        sire,
        dam,
        probSire,
        probDam,
        nLoci,
    )

    # Einstien sum notation 1: create the joint parental genotypes based on the probabilities for each parent.
    # jointParents = np.einsum("ai, bi -> abi", probSire, probDam)

    jointParents = getJointParents(probSire, probDam, nLoci)
    smooth16ByLocus(jointParents, e1e, e16, nLoci)

    # Now construct the parental genotypes based on within-family information.

    for index in range(nOffspring):
        child = family.offspring[index]

        # Einstien sum notation 2: Create the child-specific segregation tensor using the child's currrent segregation estimate.
        # childSegTensor[index,:,:,:,:] = np.einsum("abcd, di -> abci", segregationTensor, currentSeg)
        childSegs = childSegTensor[index, :, :, :, :]

        # Einstien sum notation 3: Estimate the parental genotypes based on the child's genotypes and their segregation tensor.
        # childToParents[index,:,:,:] = np.einsum("abci, ci -> abi", childSegTensor[index,:,:,:,:], childValues)
        projectChildForPeel(
            posterior,
            penetrance,
            segregation,
            peelingInfo.sex,
            peelingInfo.segregationTensor,
            peelingInfo.segregationTensorXY,
            peelingInfo.segregationTensorXX,
            child,
            childValuesCurrent,
            currentSeg,
            childSegs,
            childToParentsCurrent,
            isXChr,
            nLoci,
            e1e,
            e4,
        )
        if needsSegregationUpdate:
            childValuesTensor[index, :, :] = childValuesCurrent

        addLogChildToParentsAndSubtractCurrent(
            childToParentsCurrent,
            allToParents,
            parentsMinusChild[index, :, :, :],
            nLoci,
        )

    # Estimate the parents genotype and the child-specific posterior terms using a slightly smarter log scale.
    addJointParentsAndAllToMinus(
        parentsMinusChild, jointParents, allToParents, nOffspring, nLoci
    )

    # Move from a log-scale to a non-log scale and re-normalize.
    allToParents = expNorm2D(allToParents, nLoci)
    for i in range(nOffspring):
        parentsMinusChild[i, :, :, :] = expNorm2D(parentsMinusChild[i, :, :, :], nLoci)

    for i in range(nOffspring):
        child = family.offspring[i]

        # Einstien sum notation 4: Project the parent genotypes down onto the child genotypes.
        # anterior[child,:,:] = np.einsum("abci, abi -> ci", childSegTensor[i,:,:,:,:], parentsMinusChild[i,:,:,:])
        childAnterior = anterior[child, :, :]
        projectParentGenotypes(
            childSegTensor[i, :, :, :, :],
            parentsMinusChild[i, :, :, :],
            childAnterior,
            nLoci,
        )
        normalize4ByLocus(childAnterior, nLoci)

    if needsSegregationUpdate:
        # Estimate the segregation probabilities for each child.

        for i in range(nOffspring):
            # Child values is the same as in the posterior estimation step above.
            child = family.offspring[i]
            childValues = childValuesTensor[i, :, :]

            if isXChr:
                if peelingInfo.sex[child] == 0:  # 0 for male, 1 for female.
                    segregationTensor = peelingInfo.segregationTensorXY
                    segregationTensor_norm = peelingInfo.segregationTensorXY_norm
                elif peelingInfo.sex[child] == 1:  # 0 for male, 1 for female.
                    segregationTensor = peelingInfo.segregationTensorXX
                    segregationTensor_norm = peelingInfo.segregationTensorXX_norm

            # Einstien sum notation 5:
            # segregation[child,:,:] = np.einsum("abcd, abi, ci-> di", segregationTensor, parentsMinusChild[i,:,:,:], childValues)
            # Estimate with normalizing.
            estimateSegregationWithNorm(
                segregationTensor,
                segregationTensor_norm,
                parentsMinusChild[i, :, :, :],
                childValues,
                segregation[child, :, :],
                nLoci,
            )

            collapseSegregationInPlace(
                segregation[child, :, :],
                peelingInfo.transmissionRate,
                forwardSeg,
                nLoci,
            )
            for locus in range(nLoci):
                for state in range(4):
                    segregation[child, state, locus] = (
                        e1e * segregation[child, state, locus] + e4
                    )


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e1e": float32},
)
def peelUp(family, peelingInfo, singleLocusMode):
    """Peel information up from offspring to parents.

    :param family: The family object that the peeling is performed on.
    :type family: class:`tinyhouse.Pedigree.Family`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :param singleLocusMode: A flag to indicate the mode of peeling.
        `False` if using multi-locus peeling, and `True` if using single-locus peeling.
    :type singleLocusMode: bool
    :return: None. The function modifies the peelingInfo object in place.
    """
    isXChr = peelingInfo.isXChr

    e = 0.000001
    e1e = 1 - e
    e4 = e / 4

    anterior = peelingInfo.anterior
    penetrance = peelingInfo.penetrance
    posterior = peelingInfo.posterior
    segregation = peelingInfo.segregation

    nLoci = peelingInfo.nLoci
    nOffspring = len(family.offspring)
    sire = family.sire
    dam = family.dam
    fam = family.idn

    childToParentsCurrent = np.full((4, 4, nLoci), 0, dtype=np.float32)
    allToParents = np.full((4, 4, nLoci), 0, dtype=np.float32)
    probSire = np.full((4, nLoci), 0, dtype=np.float32)
    probDam = np.full((4, nLoci), 0, dtype=np.float32)
    childValuesCurrent = np.full((4, nLoci), 0, dtype=np.float32)
    currentSeg = np.full((4, nLoci), 1, dtype=np.float32)
    childSegTensor = np.full((1, 4, 4, 4, nLoci), 0, dtype=np.float32)
    childSegsCurrent = childSegTensor[0, :, :, :, :]

    setupParentGenotypeProbs(
        anterior,
        penetrance,
        posterior,
        peelingInfo.posteriorSireContribution[fam, :, :],
        peelingInfo.posteriorDamContribution[fam, :, :],
        sire,
        dam,
        probSire,
        probDam,
        nLoci,
    )
    smooth4ByLocus(probSire, e1e, e4, nLoci)
    smooth4ByLocus(probDam, e1e, e4, nLoci)

    for index in range(nOffspring):
        child = family.offspring[index]

        projectChildForPeel(
            posterior,
            penetrance,
            segregation,
            peelingInfo.sex,
            peelingInfo.segregationTensor,
            peelingInfo.segregationTensorXY,
            peelingInfo.segregationTensorXX,
            child,
            childValuesCurrent,
            currentSeg,
            childSegsCurrent,
            childToParentsCurrent,
            isXChr,
            nLoci,
            e1e,
            e4,
        )
        addLogChildToParents(childToParentsCurrent, allToParents, nLoci)

    allToParents = expNorm2D(allToParents, nLoci)

    sirePosterior = peelingInfo.posteriorSireContribution[fam, :, :]
    combineAndReduceAxis1(allToParents, probDam, sirePosterior, nLoci)
    normalize4ByLocus(sirePosterior, nLoci)
    smooth4ByLocus(sirePosterior, e1e, e4, nLoci)

    damPosterior = peelingInfo.posteriorDamContribution[fam, :, :]
    combineAndReduceAxis0(allToParents, probSire, damPosterior, nLoci)
    normalize4ByLocus(damPosterior, nLoci)
    smooth4ByLocus(damPosterior, e1e, e4, nLoci)


#
# The following are a large number of "helper" jit functions that replace the einstien sums in the original scripts.
#


@jit(nopython=True, nogil=True)
def setupParentGenotypeProbs(
    anterior,
    penetrance,
    posterior,
    sireContribution,
    damContribution,
    sire,
    dam,
    probSire,
    probDam,
    nLoci,
):
    """Build the normalized sire and dam genotype probabilities for a family."""

    for state in range(4):
        for locus in range(nLoci):
            probSire[state, locus] = np.log(posterior[sire, state, locus]) - np.log(
                sireContribution[state, locus]
            )
            probDam[state, locus] = np.log(posterior[dam, state, locus]) - np.log(
                damContribution[state, locus]
            )

    expNorm1DInPlace(probSire, nLoci)
    expNorm1DInPlace(probDam, nLoci)

    for state in range(4):
        for locus in range(nLoci):
            probSire[state, locus] *= (
                anterior[sire, state, locus] * penetrance[sire, state, locus]
            )
            probDam[state, locus] *= (
                anterior[dam, state, locus] * penetrance[dam, state, locus]
            )

    normalize4ByLocus(probSire, nLoci)
    normalize4ByLocus(probDam, nLoci)


@jit(
    nopython=True,
    nogil=True,
    locals={"e4": float32, "e1e": float32},
)
def projectChildForPeel(
    posterior,
    penetrance,
    segregation,
    sex,
    segregationTensorDefault,
    segregationTensorXY,
    segregationTensorXX,
    child,
    childValues,
    currentSeg,
    childSegs,
    childToParents,
    isXChr,
    nLoci,
    e1e,
    e4,
):
    """Project one child onto parental genotypes for either peeling direction."""

    for state in range(4):
        for locus in range(nLoci):
            childValues[state, locus] = (
                posterior[child, state, locus] * penetrance[child, state, locus]
            )
            currentSeg[state, locus] = segregation[child, state, locus]

    normalize4ByLocus(childValues, nLoci)
    smooth4ByLocus(childValues, e1e, e4, nLoci)
    normalize4ByLocus(currentSeg, nLoci)

    segregationTensor = segregationTensorDefault
    if isXChr:
        if sex[child] == 0:  # 0 for male, 1 for female.
            segregationTensor = segregationTensorXY
        elif sex[child] == 1:  # 0 for male, 1 for female.
            segregationTensor = segregationTensorXX

    createChildSegs(segregationTensor, currentSeg, childSegs, nLoci)
    projectChildGenotypes(childSegs, childValues, childToParents, nLoci)


@jit(nopython=True, nogil=True)
def getJointParents(probSire, probDam, nLoci):
    """Creates the joint parental genotypes based on the probabilities for each parent.

    :param probSire: the probability of each genotype of each locus of the sire
        with information of the current child from previous peeling cycle (P(p))
    :type probSire: 2D numpy array of float32 with size 4 x nLoci
    :param probDam: the probability of each genotype of each locus of the dam
        with information of the current child from previous peeling cycle (P(m))
    :type probDam: 2D numpy array of float32 with size 4 x nLoci
    :return: the probabilities of all the combinations of the sire's genotype and the dam's genotype
        with information from previous peeling cycle (P(p, m))
    :rtype: 3D numpy array of float32 with size 4 x 4 x nLoci
    """
    # jointParents = np.einsum("ai, bi -> abi", probSire, probDam)
    output = np.full(shape=(4, 4, nLoci), fill_value=0, dtype=np.float32)
    for a in range(4):
        for b in range(4):
            for i in range(nLoci):
                output[a, b, i] = probSire[a, i] * probDam[b, i]
    return output


@jit(nopython=True, nogil=True)
def createChildSegs(segregationTensor, currentSeg, output, nLoci):
    """Creates the child-specific segregation tensor using the child's current segregation estimate.

    :param segregationTensor: the probability of each combination of the sire's genotype, the dam's genotype,
        child's genotype and segregation without any other information (P(p, m, allele, seg))
    :type segregationTensor: 4D numpy array of float32 with size 4 x 4 x 4 x 4
    :param currentSeg: The probability of each segregation of each locus of the child
        with information from previous peeling cycle (P(seg))
    :type currentSeg: 2D numpy array of float32 with size 4 x nLoci
    :param output: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle (P(p, m, allele))
    :type output: 4D numpy array of float32 with size 4 x 4 x 4 x nLoci
    """
    # childSegs[index,:,:,:,:] = np.einsum("abcd, di -> abci", segregationTensor, currentSeg)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for i in range(nLoci):
                    total = 0
                    for d in range(4):
                        total += segregationTensor[a, b, c, d] * currentSeg[d, i]
                    output[a, b, c, i] = total

    return output


@jit(nopython=True, nogil=True)
def projectChildGenotypes(childSegs, childValues, output, nLoci):
    """Estimate the parental genotypes based on the child's genotypes and their segregation tensor.

    :param childSegs: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle
        (P(p, m, allele))
    :type childSegs: 4D numpy array of float32 with size 4 x 4 x 4 x nLoci
    :param childValues: the probability of each genotype of each locus of the current child
        given the information of itself and its later generations from previous peeling cycle
        (P(allele))
    :type childValues: 2D numpy array of float32 with size 4 x nLoci
    :param output: the probability of each combination of the sire's genotype, the dam's genotype
        given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type output: 3D numpy array of float32 with size 4 x 4 x nLoci
    """
    # childToParents[index,:,:,:] = np.einsum("abci, ci -> abi", childSegs[index,:,:,:,:], childValues)
    for a in range(4):
        for b in range(4):
            for i in range(nLoci):
                total = 0
                for c in range(4):
                    total += childSegs[a, b, c, i] * childValues[c, i]
                output[a, b, i] = total
    return output


@jit(nopython=True, nogil=True)
def projectParentGenotypes(childSegs, parentValues, output, nLoci):
    """Project the parent genotypes down onto the child genotypes.

    :param childSegs: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle
        (P(p, m, allele))
    :type childSegs: 4D numpy array of float32 with size 4 x 4 x 4 x nLoci
    :param parentValues: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        without the information of the current child (P(p, m))
    :type parentValues: 3D numpy array of float32 with size 4 x 4 x nLoci
    :param output: the probability of child's genotype of each locus given the later and current generations
        from previous peeling cycle without the information of the current child (P(allele))
    :type output: 2D numpy array of float32 with size 4 x nLoci
    """
    # anterior[child,:,:] = np.einsum("abci, abi -> ci", childSegs[i,:,:,:,:], parentsMinusChild[i,:,:,:])

    for c in range(4):
        for i in range(nLoci):
            total = 0
            for a in range(4):
                for b in range(4):
                    total += childSegs[a, b, c, i] * parentValues[a, b, i]
            output[c, i] = total

    return output


@jit(nopython=True, nogil=True)
def addLogChildToParents(childToParents, allToParents, nLoci):
    """Add one child's parent projection to the log-scale family total."""

    for a in range(4):
        for b in range(4):
            for i in range(nLoci):
                allToParents[a, b, i] += np.log(childToParents[a, b, i])


@jit(nopython=True, nogil=True)
def addLogChildToParentsAndSubtractCurrent(
    childToParents, allToParents, parentsMinusCurrentChild, nLoci
):
    """Add one child to the family total and store its negative log contribution."""

    for a in range(4):
        for b in range(4):
            for i in range(nLoci):
                logValue = np.log(childToParents[a, b, i])
                allToParents[a, b, i] += logValue
                parentsMinusCurrentChild[a, b, i] = -logValue


@jit(nopython=True, nogil=True)
def addJointParentsAndAllToMinus(
    parentsMinusChild, jointParents, allToParents, nOffspring, nLoci
):
    """Complete parent-minus-child log terms after all children are accumulated."""

    for childIndex in range(nOffspring):
        for a in range(4):
            for b in range(4):
                for i in range(nLoci):
                    parentsMinusChild[childIndex, a, b, i] += (
                        np.log(jointParents[a, b, i]) + allToParents[a, b, i]
                    )


@jit(nopython=True, nogil=True)
def estimateSegregationWithNorm(
    segregationTensor, segregationTensor_norm, parentValues, childValues, output, nLoci
):
    """Estimate with normalizing.

    :param segregationTensor: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype and segregation without any other information (P(p, m, allele, seg))
    :type segregationTensor: 4D numpy array of float32 with size 4 x 4 x 4 x 4
    :param segregationTensor_norm: the mean probability of each combination of the sire's genotype, the dam's genotype
        and the child's genotype across the child's segregation without any other information
        (1 / 4 x (P(p, m, allele)))
    :type segregationTensor_norm: 3D numpy array of float32 with size 4 x 4 x 4
    :param parentValues: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        without the information of the current child (P(p, m))
    :type parentValues: 3D numpy array of float32 with size 4 x 4 x nLoci
    :param childValues: the probability of each genotype of each locus of the current child
        given the information of itself and its later generations from previous peeling cycle
        (P(allele))
    :type childValues: 2D numpy array of float32 with size 4 x nLoci
    :param output: the probability of each segregation states of each locus of the current child
        given the information of later and current generations from previous peeling cycle
        (P(seg))
    :type output: 2D numpy array of float32 with size 4 x nLoci
    """
    # output = np.einsum("abcd, abi, ci-> di", segregationTensor, parentsMinusChild, childValues)
    for d in range(4):
        for i in range(nLoci):
            total = 0
            for a in range(4):
                for b in range(4):
                    for c in range(4):
                        # Check if norm is 0. Otherwise use norm to normalize.
                        if segregationTensor_norm[a, b, c] != 0:
                            total += (
                                segregationTensor[a, b, c, d]
                                * parentValues[a, b, i]
                                * childValues[c, i]
                                / segregationTensor_norm[a, b, c]
                            )
            output[d, i] = total
    return output


@jit(nopython=True, nogil=True)
def combineAndReduceAxis1(jointEstimate, parentEstimate, output, nLoci):
    """Summing over axis 1 of jointEstimate with weights given by parentEstimate

    :param jointEstimate: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type jointEstimate: 3D numpy array of float32 with size 4 x 4 x nLoci
    :param parentEstimate: the probability of each genotype of each locus of the dam
        with information of the current child from previous peeling cycle (P(m))
    :type parentEstimate: 2D numpy array of float32 with size 4 x nLoci
    :return: the probability of each genotype of each locus of the sire
        given the information of later and current generations from previous peeling cycle (P(p))
    :rtype: 2D numpy array of float32 with size 4 x nLoci
    """
    # output = np.einsum("abi, bi-> ai", jointEstimate, parentEstimate)
    for a in range(4):
        for i in range(nLoci):
            total = 0
            for b in range(4):
                total += jointEstimate[a, b, i] * parentEstimate[b, i]
            output[a, i] = total
    return output


@jit(nopython=True, nogil=True)
def combineAndReduceAxis0(jointEstimate, parentEstimate, output, nLoci):
    """Summing over axis 0 of jointEstimate with weights given by parentEstimate

    :param jointEstimate: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type jointEstimate: 3D numpy array of float32 with size 4 x 4 x nLoci
    :param parentEstimate: the probability of each genotype of each locus of the sire
        with information of the current child from previous peeling cycle (P(p))
    :type parentEstimate: 2D numpy array of float32 with size 4 x nLoci
    :return: the probability of each genotype of each locus of the dam
        given the information of later and current generations from previous peeling cycle (P(m))
    :rtype: 2D numpy array of float32 with size 4 x nLoci
    """
    # output = np.einsum("abi, ai-> bi", jointEstimate, parentEstimate)
    for b in range(4):
        for i in range(nLoci):
            total = 0
            for a in range(4):
                total += jointEstimate[a, b, i] * parentEstimate[a, i]
            output[b, i] = total
    return output


@jit(nopython=True, nogil=True)
def normalize4ByLocus(mat, nLoci):
    """Normalize a 4 x nLoci probability matrix in place."""

    for i in range(nLoci):
        total = 0
        for a in range(4):
            total += mat[a, i]
        for a in range(4):
            mat[a, i] /= total
    return mat


@jit(nopython=True, nogil=True)
def normalize16ByLocus(mat, nLoci):
    """Normalize a 4 x 4 x nLoci probability tensor in place."""

    for i in range(nLoci):
        total = 0
        for a in range(4):
            for b in range(4):
                total += mat[a, b, i]
        for a in range(4):
            for b in range(4):
                mat[a, b, i] /= total
    return mat


@jit(nopython=True, nogil=True)
def smooth4ByLocus(mat, scale, floor, nLoci):
    """Apply a small in-place probability floor to a 4 x nLoci matrix."""

    for a in range(4):
        for i in range(nLoci):
            mat[a, i] = mat[a, i] * scale + floor
    return mat


@jit(nopython=True, nogil=True)
def smooth16ByLocus(mat, scale, floor, nLoci):
    """Apply a small in-place probability floor to a 4 x 4 x nLoci tensor."""

    for a in range(4):
        for b in range(4):
            for i in range(nLoci):
                mat[a, b, i] = mat[a, b, i] * scale + floor
    return mat


@jit(nopython=True, nogil=True)
def expNorm2D(mat, nLoci):
    """Output is to take the exponential of the matrix and normalize each locus.

    :param mat: a 3D matrix with last axis represents the locus
    :type mat: 3D numpy array of float32 with size 4 x 4 x nLoci
    :return: the normalized exponential of the `mat`
    :rtype: 3D numpy array of float32 with size 4 x 4 x nLoci
    """
    # Matrix is 4x4xnLoci: Output is to take the exponential of the matrix and normalize each locus. We need to make sure that there are not any overflow values.
    for i in range(nLoci):
        maxVal = (
            1  # Log of anything between 0-1 will be less than 0. Using 1 as a default.
        )
        for a in range(4):
            for b in range(4):
                if mat[a, b, i] > maxVal or maxVal == 1:
                    maxVal = mat[a, b, i]
        for a in range(4):
            for b in range(4):
                mat[a, b, i] -= maxVal
    tmp = np.exp(mat)
    normalize16ByLocus(tmp, nLoci)
    return tmp


@jit(nopython=True, nogil=True)
def expNorm1D(mat, nLoci):
    """Output is to take the exponential of the matrix and normalize each locus.

    :param mat: a 2D matrix with last axis represents the locus
    :type mat: 2D numpy array of float32 with size 4 x nLoci
    :return: the normalized exponential of the `mat`
    :rtype: 2D numpy array of float32 with size 4 x nLoci
    """
    # Matrix is 4x4xnLoci: Output is to take the exponential of the matrix and normalize each locus. We need to make sure that there are not any overflow values.
    for i in range(nLoci):
        maxVal = (
            1  # Log of anything between 0-1 will be less than 0. Using 1 as a default.
        )
        for a in range(4):
            if mat[a, i] > maxVal or maxVal == 1:
                maxVal = mat[a, i]
        for a in range(4):
            mat[a, i] -= maxVal
    tmp = np.exp(mat)
    normalize4ByLocus(tmp, nLoci)
    return tmp


@jit(nopython=True, nogil=True)
def expNorm1DInPlace(mat, nLoci):
    """Take the exponential of a 4 x nLoci log matrix and normalize in place."""

    for i in range(nLoci):
        maxVal = 1
        for a in range(4):
            if mat[a, i] > maxVal or maxVal == 1:
                maxVal = mat[a, i]
        total = 0
        for a in range(4):
            mat[a, i] = np.exp(mat[a, i] - maxVal)
            total += mat[a, i]
        for a in range(4):
            mat[a, i] /= total
    return mat


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e2": float32, "e1e": float32, "e2i": float32},
)
def collapseSegregationInPlace(segregation, transmission, forward, nLoci):
    """Using Baum-Welch algorithm to calculate the segregation probabilities.

    :param segregation: the probability of each segregation states of each locus of the current child
        given the information of later and current generations from previous peeling cycle
        (P(seg)). This is updated in place with the collapsed probabilities.
        the segregation state ordering: pp, pm, mp, mm
    :type segregation: 2D numpy array of float32 with size 4 x nLoci
    :param transmission: transmission function based on the distance between adjacent loci
    :type transmission: 1D numpy array of float32 with size (nLoci - 1)
    :param forward: work array for forward messages.
    :type forward: 2D numpy array of float32 with size 4 x nLoci
    :return: None. The function updates segregation in place with the collapsed probabilities.
    """
    # This is the forward backward algorithm.
    # Segregation estimate state ordering: pp, pm, mp, mm
    tmp = np.full((4), 0, dtype=np.float32)
    new = np.full((4), 0, dtype=np.float32)

    prev = np.full((4), 1, dtype=np.float32)
    for j in range(4):
        forward[j, 0] = 1

    for i in range(1, nLoci):
        e = transmission[i - 1]
        e2 = e**2
        e1e = e * (1 - e)
        e2i = (1.0 - e) ** 2
        for j in range(4):
            tmp[j] = prev[j] * segregation[j, i - 1]

        sum_j = 0
        for j in range(4):
            sum_j += tmp[j]
        for j in range(4):
            tmp[j] = tmp[j] / sum_j

        # !                  fm  fm  fm  fm
        # !segregationOrder: pp, pm, mp, mm

        new[0] = e2 * tmp[3] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[0]
        new[1] = e2 * tmp[2] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[1]
        new[2] = e2 * tmp[1] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[2]
        new[3] = e2 * tmp[0] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[3]

        for j in range(4):
            forward[j, i] = new[j]
            prev[j] = new[j]

    for j in range(4):
        prev[j] = 1

    for i in range(
        nLoci - 2, -1, -1
    ):  # zero indexed then minus one since we skip the boundary.
        e = transmission[i]
        e2 = e**2
        e1e = e * (1 - e)
        e2i = (1.0 - e) ** 2

        for j in range(4):
            tmp[j] = prev[j] * segregation[j, i + 1]

        sum_j = 0
        for j in range(4):
            sum_j += tmp[j]
        for j in range(4):
            tmp[j] = tmp[j] / sum_j

        new[0] = e2 * tmp[3] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[0]
        new[1] = e2 * tmp[2] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[1]
        new[2] = e2 * tmp[1] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[2]
        new[3] = e2 * tmp[0] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[3]

        sum_j = 0
        for j in range(4):
            segregation[j, i + 1] = segregation[j, i + 1] * forward[j, i + 1] * prev[j]
            sum_j += segregation[j, i + 1]
        for j in range(4):
            segregation[j, i + 1] = segregation[j, i + 1] / sum_j
            prev[j] = new[j]

    sum_j = 0
    for j in range(4):
        segregation[j, 0] = segregation[j, 0] * forward[j, 0] * prev[j]
        sum_j += segregation[j, 0]
    for j in range(4):
        segregation[j, 0] = segregation[j, 0] / sum_j

    return
