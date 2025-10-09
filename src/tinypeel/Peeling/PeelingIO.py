import numpy as np
from numba import jit
from ..tinyhouse import InputOutput


def readInSeg(pedigree, fileName, start=None, stop=None):
    """Reads in a segregation file and returns segregation probabilities.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param fileName: path to the external penetrance file
    :type fileName: str
    :param start: starting locus/marker, defaults to None
    :type start: int, optional
    :param stop: final locus/marker, defaults to None
    :type stop: int, optional
    :raises ValueError: If the length of the segregation line does not match the expected length.
    :raises ValueError: If the length of the line subsection does not match the expected length.
    :raises ValueError: If the individual in the segregation file is not found in the pedigree.
    :return: segregation probabilities for each locus and individual in the pedigree
    :rtype: 3D numpy array of float32 with shape nIndividuals x 4 x nLoci
    """
    print("Reading in seg file:", fileName)
    if start is None:
        start = 0
    if stop is None:
        stop = pedigree.nLoci
    nLoci = stop - start + 1  # Contains stop.

    seg = np.full((pedigree.maxIdn, 4, nLoci), 0.25, dtype=np.float32)

    index = 0
    fileNColumns = 0

    indHit = np.full(pedigree.maxIdn, 0, dtype=np.int64)

    with open(fileName) as f:
        e = 0
        currentInd = None
        for line in f:
            parts = line.split()
            idx = parts[0]

            if fileNColumns == 0:
                fileNColumns = len(parts)
            if fileNColumns != len(parts):
                raise ValueError(
                    f"The length of the line is not the expected length. Expected {fileNColumns} got {len(parts)} on individual {idx} and line {e}."
                )

            segLine = np.array(
                [float(val) for val in parts[(start + 1) : (stop + 2)]],
                dtype=np.float32,
            )
            if len(segLine) != nLoci:
                raise ValueError(
                    f"The length of the line subsection is not the expected length. Expected {nLoci} got {len(segLine)} on individual {idx} and line {e}."
                )

            if idx not in pedigree.individuals:
                print(f"Individual {idx} not found in pedigree. Individual ignored.")
            else:
                ind = pedigree.individuals[idx]
                if e == 0:
                    currentInd = ind.idx
                if currentInd != ind.idx:
                    raise ValueError(
                        f"Unexpected individual. Expecting individual {currentInd}, but got ind {ind.idx} on value {e}"
                    )
                seg[ind.idn, e, :] = segLine
                e = (e + 1) % 4
                ind.fileIndex["segregation"] = index
                index += 1
                indHit[ind.idn] += 1
        for ind in pedigree:
            if indHit[ind.idn] != 4:
                print(f"No segregation information found for individual {ind.idx}")

    return seg


def writeOutParamaters(peelingInfo):
    """Writes out the geno error rate, seq error rate, and recombination probabilities.

    :param peelingInfo: Peeling information container
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: None. Writes to files specified in the InputOutput.args.
    """
    args = InputOutput.args

    np.savetxt(args.out_file + ".geno_error_prob.txt", peelingInfo.genoError, fmt="%f")
    np.savetxt(args.out_file + ".seq_error_prob.txt", peelingInfo.seqError, fmt="%f")
    np.savetxt(
        args.out_file + ".rec_prob.txt", np.empty((1, 1)), fmt="%f"
    )  # not be realized, just as a placeholder
    # np.savetxt(args.out_file + ".trans", peelingInfo.transmissionRate, fmt = "%f")
    # np.savetxt(args.out_file + ".alt_allele_prob.txt", peelingInfo.maf, fmt="%f") - remove?


def writeOutAltAlleleProb(pedigree):
    """Writes out the alternative allele probabilities for each locus and metafounder in the pedigree.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :return: None. Writes to a file specified in the InputOutput.args.
    """
    args = InputOutput.args

    # Custom sorting key to extract numeric part of metafounder keys if it's an integer
    def sort_key(mf_key):
        part = mf_key.split("_")[1]
        return (0, int(part)) if part.isdigit() else (1, part)

    # Order the metafounders: MF_1, MF_2, ..., MF_11, etc., otherwise keep original order
    sorted_AAP = dict(sorted(pedigree.AAP.items(), key=lambda item: sort_key(item[0])))
    sorted_MF = list(sorted_AAP.keys())
    # Combine data into a single 2D array
    combined_AAP = np.hstack(
        [sorted_AAP[key].reshape(pedigree.nLoci, -1) for key in sorted_MF]
    )
    # Save into text file with metafounders heading columns
    np.savetxt(
        args.out_file + ".alt_allele_prob.txt",
        combined_AAP,
        delimiter="\t",
        fmt="%.2f",
        header="\t".join(sorted_MF),
        comments="",
    )


def writePhenoPenetrance(pedigree):
    """Writes out the phenotype penetrance for each individual in the pedigree.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :return: None. Writes to a file specified in the InputOutput.args.
    """
    args = InputOutput.args
    np.savetxt(
        args.out_file + ".pheno_penetrance.txt", pedigree.phenoPenetrance, fmt="%.2f"
    )


def writeGenotypes(pedigree, genoProbFunc, isSexChrom):
    """Writes out the genotypes for each individual. Format depends on user input and arguments.
    The output can include:
    - Dosages
    - Phased genotype probabilities
    - Genotype probabilities
    - Called genotypes (based on user inputted threshold or default 1/3)
    - Called haplotypes (based on user inputted threshold or default 1/2)

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param isSexChrom: flag whether inputted genotypes are sex chromosome
    :type isSexChrom: bool
    :return: None. Writes to files specified in the InputOutput.args.
    """
    args = InputOutput.args
    if not args.no_dosage:
        writeDosages(pedigree, genoProbFunc, isSexChrom, args.out_file + ".dosage.txt")
    if args.phased_geno_prob:
        writePhasedGenoProbs(
            pedigree, genoProbFunc, args.out_file + ".phased_geno_prob.txt"
        )
    if args.geno_prob:
        writeGenoProbs(pedigree, genoProbFunc, args.out_file + ".geno_prob.txt")
    if args.geno_threshold and args.geno:
        for thresh in args.geno_threshold:
            if thresh < 1 / 3:
                thresh = 1 / 3
            if args.binary_call_file:
                writeBinaryCalledGenotypes(
                    pedigree,
                    genoProbFunc,
                    isSexChrom,
                    args.out_file + ".called." + str(thresh),
                    thresh,
                )
            else:
                writeCalledGenotypes(
                    pedigree,
                    genoProbFunc,
                    isSexChrom,
                    args.out_file + ".geno_" + str(thresh) + ".txt",
                    thresh,
                )

    if args.hap_threshold and args.hap:
        for thresh in args.hap_threshold:
            if thresh < 1 / 2:
                thresh = 1 / 2
            if args.binary_call_file:
                pass  # this function is not applied
            else:
                writeCalledPhase(
                    pedigree,
                    genoProbFunc,
                    args.out_file + ".hap_" + str(thresh) + ".txt",
                    thresh,
                )


def writePhasedGenoProbs(pedigree, genoProbFunc, outputFile):
    """Writes the phased genotype probabilities to a file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param outputFile: name of output file to write to
    :type outputFile: str
    :return: None. Writes to the specified output file.
    """
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)
            f.write("\n")
            for i in range(matrix.shape[0]):
                f.write(
                    ind.idx + " " + " ".join(map("{:.4f}".format, matrix[i, :])) + "\n"
                )


def writeGenoProbs(pedigree, genoProbFunc, outputFile):
    """Writes out the non phased genotype probabilities to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param outputFile: name of output file to write to
    :type outputFile: str
    :return: None. Writes to the specified output file.
    """
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)
            f.write("\n")
            for i in range(matrix.shape[0]):
                if i == 1:  # Add up probabilities for aA and Aa
                    f.write(
                        ind.idx
                        + " "
                        + " ".join(
                            map("{:.4f}".format, matrix[i, :] + matrix[i + 1, :])
                        )
                        + "\n"
                    )
                elif i != 2:  # Print probabilities for aa and AA
                    f.write(
                        ind.idx
                        + " "
                        + " ".join(map("{:.4f}".format, matrix[i, :]))
                        + "\n"
                    )


def writePhenoProbs(pedigree, phenoProbFunc):
    """Writes out the phenotype probabilities to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :return: None. Writes to the specified output file.
    """
    args = InputOutput.args
    with open(args.out_file + ".pheno_prob.txt", "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = phenoProbFunc(ind.idn, pedigree.phenoPenetrance)
            f.write("\n")
            for i in range(matrix.shape[0]):
                f.write(
                    ind.idx + " " + " ".join(map("{:.4f}".format, matrix[i, :])) + "\n"
                )


def writeDosages(pedigree, genoProbFunc, isSexChrom, outputFile):
    """Writes out the allele dosages to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param isSexChrom: flag whether inputted genotypes are sex chromosome
    :type isSexChrom: bool
    :param outputFile: name of output file to write to
    :type outputFile: str
    :return: None. Writes to the specified output file.
    """
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = np.dot(np.array([0, 1, 1, 2]), genoProbFunc(ind.idn))

            if isSexChrom and ind.sex == 0:
                matrix *= 2

            f.write(ind.idx + " " + " ".join(map("{:.4f}".format, matrix)) + "\n")


def writeCalledGenotypes(pedigree, genoProbFunc, isSexChrom, outputFile, thresh):
    """Writes out the called genotypes to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param isSexChrom: flag whether inputted genotypes are sex chromosome
    :type isSexChrom: bool
    :param outputFile: name of output file to write to
    :type outputFile: str
    :param thresh: threshold for calling genotypes, defaults to 1/3
    :type thresh: float
    :return: None. Writes to the specified output file.
    """
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)

            matrixCollapsedHets = np.array(
                [matrix[0, :], matrix[1, :] + matrix[2, :], matrix[3, :]],
                dtype=np.float32,
            )
            calledGenotypes = np.argmax(matrixCollapsedHets, axis=0)
            setMissing(calledGenotypes, matrixCollapsedHets, thresh)
            if isSexChrom and ind.sex == 0:
                doubleIfNotMissing(calledGenotypes)

            f.write(ind.idx + " " + " ".join(map(str, calledGenotypes)) + "\n")


def writeCalledPhase(pedigree, genoProbFunc, outputFile, thresh):
    """Writes out the called haplotypes to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param outputFile: name of output file to write to
    :type outputFile: str
    :param thresh: threshold for calling genotypes, defaults to 1/3
    :type thresh: float
    :return: None. Writes to the specified output file.
    """
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)

            # Paternal
            paternal_probs = np.array(
                [matrix[0, :] + matrix[1, :], matrix[2, :] + matrix[3, :]],
                dtype=np.float32,
            )
            paternal_haplotype = np.argmax(paternal_probs, axis=0)
            setMissing(paternal_haplotype, paternal_probs, thresh)
            f.write(ind.idx + " " + " ".join(map(str, paternal_haplotype)) + "\n")

            # Maternal
            maternal_probs = np.array(
                [matrix[0, :] + matrix[2, :], matrix[1, :] + matrix[3, :]],
                dtype=np.float32,
            )
            maternal_haplotype = np.argmax(maternal_probs, axis=0)
            setMissing(maternal_haplotype, maternal_probs, thresh)
            f.write(ind.idx + " " + " ".join(map(str, maternal_haplotype)) + "\n")


def writeBinaryCalledGenotypes(pedigree, genoProbFunc, isSexChrom, outputFile, thresh):
    """Writes out the called genotypes to file in binary format.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param isSexChrom: flag whether inputted genotypes are sex chromosome
    :type isSexChrom: bool
    :param outputFile: name of output file to write to
    :type outputFile: str
    :param thresh: threshold for calling genotypes, defaults to 1/3
    :type thresh: float
    :return: None. Writes to the specified output file.
    """
    for idx, ind in pedigree.writeOrder():
        matrix = genoProbFunc(ind.idn)
        matrixCollapsedHets = np.array(
            [matrix[0, :], matrix[1, :] + matrix[2, :], matrix[3, :]], dtype=np.float32
        )
        calledGenotypes = np.argmax(matrixCollapsedHets, axis=0)
        setMissing(calledGenotypes, matrixCollapsedHets, thresh)
        if isSexChrom and ind.sex == 0:
            doubleIfNotMissing(calledGenotypes)
        ind.genotypes = calledGenotypes.astype(np.int8)

    InputOutput.writeOutGenotypesPlink(pedigree, outputFile)


@jit(nopython=True)
def doubleIfNotMissing(calledGenotypes):
    """Doubles the called genotype, used for (female) sex chromosomes.

    :param calledGenotypes: array of called genotypes
    :type calledGenotypes: 2D numpy array of float32 with shape 3 x nLoci
    :return: None. Modifies the calledGenotypes in place.
    """
    nLoci = len(calledGenotypes)
    for i in range(nLoci):
        if calledGenotypes[i] == 1:
            calledGenotypes[i] = 2


@jit(nopython=True)
def setMissing(calledGenotypes, matrix, thresh):
    """Sets the called genotypes to missing if the probability is below the threshold.

    :param calledGenotypes: array of called genotypes
    :type calledGenotypes: 2D numpy array of float32 with shape 3 x nLoci
    :param matrix: genotype probability matrix
    :type matrix: 2D numpy array of float32 with shape 3 x nLoci
    :param thresh: threshold for calling genotypes, defaults to 1/3
    :type thresh: float
    """
    nLoci = len(calledGenotypes)
    for i in range(nLoci):
        if matrix[calledGenotypes[i], i] <= thresh:
            calledGenotypes[i] = 9


def fullOutput(pedigree, peelingInfo, args):
    """Writes out the full output of the peeling process (for testing purposes).

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :param args: argument container with configuration options for peeling
    :type args: argparse.Namespace or similar object with attributes
    """
    InputOutput.writeIdnIndexedMatrix(
        pedigree, peelingInfo.penetrance, args.out_file + ".penetrance"
    )
    InputOutput.writeIdnIndexedMatrix(
        pedigree, peelingInfo.anterior, args.out_file + ".anterior"
    )
    InputOutput.writeIdnIndexedMatrix(
        pedigree, peelingInfo.posterior, args.out_file + ".posterior"
    )

    InputOutput.writeFamIndexedMatrix(
        pedigree,
        peelingInfo.posteriorSire_minusFam,
        args.out_file + ".posteriorSire_minusFam",
    )
    InputOutput.writeFamIndexedMatrix(
        pedigree,
        peelingInfo.posteriorDam_minusFam,
        args.out_file + ".posteriorDam_minusFam",
    )

    InputOutput.writeFamIndexedMatrix(
        pedigree, peelingInfo.posteriorSire_new, args.out_file + ".posteriorSire_new"
    )
    InputOutput.writeFamIndexedMatrix(
        pedigree, peelingInfo.posteriorDam_new, args.out_file + ".posteriorDam_new"
    )
