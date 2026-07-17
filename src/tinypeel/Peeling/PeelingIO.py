import numpy as np
from numba import jit
from contextlib import ExitStack
from ..tinyhouse import InputOutput


def writeOutParamaters(peelingInfo):
    """Writes out the geno error rate, seq error rate, and recombination probabilities.

    :param peelingInfo: Peeling information container
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: None. Writes to files specified in the InputOutput.args.
    """
    args = InputOutput.args
    if args.est_geno_error_prob:
        np.savetxt(
            args.out_file + ".geno_error_prob.txt",
            peelingInfo.genoError,
        )
    if args.est_seq_error_prob:
        np.savetxt(
            args.out_file + ".seq_error_prob.txt",
            peelingInfo.seqError,
        )
    if args.rec_prob:
        np.savetxt(
            args.out_file + ".rec_prob.txt",
            np.empty((1, 1)),
        )  # not implemented atm, just as a placeholder


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
        args.out_file + ".pheno_penetrance.txt",
        pedigree.phenoPenetrance,
    )


def writeGenotypes(pedigree, genoProbFunc, isXChr):
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
    :param isXChr: flag whether inputted genotypes are X chromosome
    :type isXChr: bool
    :return: None. Writes to files specified in the InputOutput.args.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    xChrMaleDosageWeights = np.array([0, 0, 0, 1])
    autosomeDosageWeights = np.array([0, 1, 1, 2])

    geno_threshold_list = []
    if args.geno:
        if args.geno_threshold:
            for thresh in args.geno_threshold:
                if thresh < 1 / 3:
                    geno_threshold_list.append(1 / 3)
                else:
                    geno_threshold_list.append(thresh)
        else:
            geno_threshold_list.append(1 / 3)

    hap_threshold_list = []
    if args.hap:
        if args.hap_threshold:
            for thresh in args.hap_threshold:
                if thresh < 1 / 2:
                    hap_threshold_list.append(1 / 2)
                else:
                    hap_threshold_list.append(thresh)
        else:
            hap_threshold_list.append(1 / 2)

    hasOutput = (
        (not args.no_dosage)
        or args.phased_geno_prob
        or args.geno_prob
        or len(geno_threshold_list) != 0
        or len(hap_threshold_list) != 0
    )
    if not hasOutput:
        return

    with ExitStack() as stack:
        dosageFile = None
        if not args.no_dosage:
            dosageFile = stack.enter_context(open(args.out_file + ".dosage.txt", "w+"))

        phasedGenoProbFile = None
        if args.phased_geno_prob:
            phasedGenoProbFile = stack.enter_context(
                open(args.out_file + ".phased_geno_prob.txt", "w+")
            )

        genoProbFile = None
        if args.geno_prob:
            genoProbFile = stack.enter_context(
                open(args.out_file + ".geno_prob.txt", "w+")
            )

        genoFiles = []
        for threshold in geno_threshold_list:
            outputFile = args.out_file + ".geno_" + str(round(threshold, 3)) + ".txt"
            print(
                f"Writing called genotypes with threshold {threshold} to {outputFile}"
            )
            genoFiles.append((threshold, stack.enter_context(open(outputFile, "w+"))))

        hapFiles = []
        for threshold in hap_threshold_list:
            outputFile = args.out_file + ".hap_" + str(round(threshold, 3)) + ".txt"
            print(
                f"Writing called haplotypes with threshold {threshold} to {outputFile}"
            )
            hapFiles.append((threshold, stack.enter_context(open(outputFile, "w+"))))

        if isXChr:
            for idx, ind in pedigree.writeOrder():
                matrix = genoProbFunc(ind.idn, ind.sex)

                if dosageFile is not None:
                    writeDosageFromMatrix(
                        dosageFile,
                        ind,
                        matrix,
                        True,
                        formatter,
                        xChrMaleDosageWeights,
                        autosomeDosageWeights,
                    )

                if phasedGenoProbFile is not None:
                    writePhasedGenoProbsFromMatrix(
                        phasedGenoProbFile, ind, matrix, formatter
                    )

                if genoProbFile is not None:
                    writeGenoProbsFromMatrix(genoProbFile, ind, matrix, formatter)

                if genoFiles:
                    matrixCollapsedHets = getCollapsedGenotypes(matrix, True, ind.sex)
                    for threshold, outputHandle in genoFiles:
                        writeCalledGenotypesFromCollapsed(
                            outputHandle, ind, matrixCollapsedHets, threshold
                        )

                if hapFiles:
                    for threshold, outputHandle in hapFiles:
                        writeCalledPhaseFromMatrix(
                            outputHandle, ind, matrix, True, threshold
                        )
        else:
            for idx, ind in pedigree.writeOrder():
                matrix = genoProbFunc(ind.idn, ind.sex)

                if dosageFile is not None:
                    writeAutosomeDosageFromMatrix(
                        dosageFile,
                        ind,
                        matrix,
                        formatter,
                        autosomeDosageWeights,
                    )

                if phasedGenoProbFile is not None:
                    writePhasedGenoProbsFromMatrix(
                        phasedGenoProbFile, ind, matrix, formatter
                    )

                if genoProbFile is not None:
                    writeGenoProbsFromMatrix(genoProbFile, ind, matrix, formatter)

                if genoFiles:
                    matrixCollapsedHets = getAutosomeCollapsedGenotypes(matrix)
                    for threshold, outputHandle in genoFiles:
                        writeCalledGenotypesFromCollapsed(
                            outputHandle, ind, matrixCollapsedHets, threshold
                        )

                if hapFiles:
                    for threshold, outputHandle in hapFiles:
                        writeAutosomeCalledPhaseFromMatrix(
                            outputHandle, ind, matrix, threshold
                        )


def writeGenotypesSeparatePasses(pedigree, genoProbFunc, isXChr):
    """Writes genotype outputs with one full pedigree pass per output file."""

    args = InputOutput.args
    if not args.no_dosage:
        writeDosages(pedigree, genoProbFunc, isXChr, args.out_file + ".dosage.txt")
    if args.phased_geno_prob:
        writePhasedGenoProbs(
            pedigree, genoProbFunc, args.out_file + ".phased_geno_prob.txt"
        )
    if args.geno_prob:
        writeGenoProbs(pedigree, genoProbFunc, args.out_file + ".geno_prob.txt")
    if args.geno:
        geno_threshold_list = []
        if args.geno_threshold:
            for thresh in args.geno_threshold:
                if thresh < 1 / 3:
                    geno_threshold_list.append(1 / 3)
                else:
                    geno_threshold_list.append(thresh)
        else:
            geno_threshold_list.append(1 / 3)

        for threshold in geno_threshold_list:
            print(
                f"Writing called genotypes with threshold {threshold} to {args.out_file + '.geno_' + str(round(threshold, 3)) + '.txt'}"
            )
            writeCalledGenotypes(
                pedigree,
                genoProbFunc,
                isXChr,
                args.out_file + ".geno_" + str(round(threshold, 3)) + ".txt",
                threshold,
            )

    if args.hap:
        hap_threshold_list = []
        if args.hap_threshold:
            for thresh in args.hap_threshold:
                if thresh < 1 / 2:
                    hap_threshold_list.append(1 / 2)
                else:
                    hap_threshold_list.append(thresh)
        else:
            hap_threshold_list.append(1 / 2)

        for threshold in hap_threshold_list:
            print(
                f"Writing called haplotypes with threshold {threshold} to {args.out_file + '.hap_' + str(round(threshold, 3)) + '.txt'}"
            )
            writeCalledPhase(
                pedigree,
                genoProbFunc,
                isXChr,
                args.out_file + ".hap_" + str(round(threshold, 3)) + ".txt",
                threshold,
            )


def writePhasedGenoProbsFromMatrix(f, ind, matrix, formatter):
    """Writes one individual's phased genotype probabilities."""

    for i in range(matrix.shape[0]):
        matrixRow = matrix[i, :]
        f.write(ind.idx + " " + " ".join(map(formatter, matrixRow)) + "\n")


def writeGenoProbsFromMatrix(f, ind, matrix, formatter):
    """Writes one individual's unphased genotype probabilities."""

    matrix0 = matrix[0, :]
    matrix1 = matrix[1, :]
    matrix2 = matrix[2, :]
    matrix3 = matrix[3, :]
    for i in range(matrix.shape[0]):
        if i == 1:  # Add up probabilities for aA and Aa
            f.write(ind.idx + " " + " ".join(map(formatter, matrix1 + matrix2)) + "\n")
        elif i != 2:  # Print probabilities for aa and AA
            matrixRow = matrix0
            if i == 3:
                matrixRow = matrix3
            f.write(ind.idx + " " + " ".join(map(formatter, matrixRow)) + "\n")


def writeDosageFromMatrix(
    f,
    ind,
    matrix,
    isXChr,
    formatter,
    xChrMaleDosageWeights,
    autosomeDosageWeights,
):
    """Writes one individual's allele dosage."""

    if isXChr and ind.sex == 0:
        tmp = xChrMaleDosageWeights
    else:
        tmp = autosomeDosageWeights
    dosage = np.dot(tmp, matrix)
    f.write(ind.idx + " " + " ".join(map(formatter, dosage)) + "\n")


def writeAutosomeDosageFromMatrix(f, ind, matrix, formatter, dosageWeights):
    """Writes one individual's autosomal allele dosage."""

    dosage = np.dot(dosageWeights, matrix)
    f.write(ind.idx + " " + " ".join(map(formatter, dosage)) + "\n")


def getCollapsedGenotypes(matrix, isXChr, sex):
    """Collapse phased genotype probabilities into called-genotype states."""

    matrix0 = matrix[0, :]
    matrix1 = matrix[1, :]
    matrix2 = matrix[2, :]
    matrix3 = matrix[3, :]
    if isXChr and sex == 0:
        matrixCollapsedHets = np.empty((2, matrix.shape[1]), dtype=np.float32)
        matrixCollapsedHets[0, :] = matrix0 + matrix2
        matrixCollapsedHets[1, :] = matrix1 + matrix3
    else:
        matrixCollapsedHets = np.empty((3, matrix.shape[1]), dtype=np.float32)
        matrixCollapsedHets[0, :] = matrix0
        matrixCollapsedHets[1, :] = matrix1 + matrix2
        matrixCollapsedHets[2, :] = matrix3
    return matrixCollapsedHets


def getAutosomeCollapsedGenotypes(matrix):
    """Collapse autosomal phased genotype probabilities into called-genotype states."""

    matrixCollapsedHets = np.empty((3, matrix.shape[1]), dtype=np.float32)
    matrixCollapsedHets[0, :] = matrix[0, :]
    matrixCollapsedHets[1, :] = matrix[1, :] + matrix[2, :]
    matrixCollapsedHets[2, :] = matrix[3, :]
    return matrixCollapsedHets


def writeCalledGenotypesFromCollapsed(f, ind, matrixCollapsedHets, thresh):
    """Writes one individual's called genotypes from collapsed probabilities."""

    calledGenotypes = np.argmax(matrixCollapsedHets, axis=0)
    setMissing(calledGenotypes, matrixCollapsedHets, thresh)
    f.write(ind.idx + " " + " ".join(map(str, calledGenotypes)) + "\n")


def writeCalledPhaseFromMatrix(f, ind, matrix, isXChr, thresh):
    """Writes one individual's called haplotypes."""

    matrix0 = matrix[0, :]
    matrix1 = matrix[1, :]
    matrix2 = matrix[2, :]
    matrix3 = matrix[3, :]

    if isXChr and ind.sex == 0:
        paternal_haplotype = np.full(matrix.shape[1], 9, dtype=np.int8)
    else:
        paternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
        paternal_probs[0, :] = matrix0 + matrix1
        paternal_probs[1, :] = matrix2 + matrix3
        paternal_haplotype = np.argmax(paternal_probs, axis=0)
        setMissing(paternal_haplotype, paternal_probs, thresh)
    f.write(ind.idx + " " + " ".join(map(str, paternal_haplotype)) + "\n")

    maternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
    maternal_probs[0, :] = matrix0 + matrix2
    maternal_probs[1, :] = matrix1 + matrix3
    maternal_haplotype = np.argmax(maternal_probs, axis=0)
    setMissing(maternal_haplotype, maternal_probs, thresh)
    f.write(ind.idx + " " + " ".join(map(str, maternal_haplotype)) + "\n")


def writeAutosomeCalledPhaseFromMatrix(f, ind, matrix, thresh):
    """Writes one individual's autosomal called haplotypes."""

    matrix0 = matrix[0, :]
    matrix1 = matrix[1, :]
    matrix2 = matrix[2, :]
    matrix3 = matrix[3, :]

    paternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
    paternal_probs[0, :] = matrix0 + matrix1
    paternal_probs[1, :] = matrix2 + matrix3
    paternal_haplotype = np.argmax(paternal_probs, axis=0)
    setMissing(paternal_haplotype, paternal_probs, thresh)
    f.write(ind.idx + " " + " ".join(map(str, paternal_haplotype)) + "\n")

    maternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
    maternal_probs[0, :] = matrix0 + matrix2
    maternal_probs[1, :] = matrix1 + matrix3
    maternal_haplotype = np.argmax(maternal_probs, axis=0)
    setMissing(maternal_haplotype, maternal_probs, thresh)
    f.write(ind.idx + " " + " ".join(map(str, maternal_haplotype)) + "\n")


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
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn, ind.sex)
            for i in range(matrix.shape[0]):
                matrixRow = matrix[i, :]
                f.write(ind.idx + " " + " ".join(map(formatter, matrixRow)) + "\n")


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
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn, ind.sex)
            matrix0 = matrix[0, :]
            matrix1 = matrix[1, :]
            matrix2 = matrix[2, :]
            matrix3 = matrix[3, :]
            for i in range(matrix.shape[0]):
                if i == 1:  # Add up probabilities for aA and Aa
                    f.write(
                        ind.idx
                        + " "
                        + " ".join(
                            map(
                                formatter,
                                matrix1 + matrix2,
                            )
                        )
                        + "\n"
                    )
                elif i != 2:  # Print probabilities for aa and AA
                    matrixRow = matrix0
                    if i == 3:
                        matrixRow = matrix3
                    f.write(ind.idx + " " + " ".join(map(formatter, matrixRow)) + "\n")


def writePhenoProbs(pedigree, phenoProbFunc):
    """Writes out the phenotype probabilities to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :return: None. Writes to the specified output file.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    with open(args.out_file + ".pheno_prob.txt", "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = phenoProbFunc(ind.idn, pedigree.phenoPenetrance)
            f.write("\n")
            for i in range(matrix.shape[0]):
                matrixRow = matrix[i, :]
                f.write(ind.idx + " " + " ".join(map(formatter, matrixRow)) + "\n")


def writeDosages(pedigree, genoProbFunc, isXChr, outputFile):
    """Writes out the allele dosages to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param isXChr: flag whether inputted genotypes are X chromosome
    :type isXChr: bool
    :param outputFile: name of output file to write to
    :type outputFile: str
    :return: None. Writes to the specified output file.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    xChrMaleDosageWeights = np.array([0, 0, 0, 1])
    autosomeDosageWeights = np.array([0, 1, 1, 2])
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            if isXChr and ind.sex == 0:
                tmp = xChrMaleDosageWeights
            else:
                tmp = autosomeDosageWeights
            matrix = np.dot(tmp, genoProbFunc(ind.idn, ind.sex))
            f.write(ind.idx + " " + " ".join(map(formatter, matrix)) + "\n")


def writeCalledGenotypes(pedigree, genoProbFunc, isXChr, outputFile, thresh):
    """Writes out the called genotypes to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param genoProbFunc: function to get genotype probabilities for an individual
    :type genoProbFunc: function
    :param isXChr: flag whether inputted genotypes are X chromosome
    :type isXChr: bool
    :param outputFile: name of output file to write to
    :type outputFile: str
    :param thresh: threshold for calling genotypes, defaults to 1/3
    :type thresh: float
    :return: None. Writes to the specified output file.
    """
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn, ind.sex)
            matrix0 = matrix[0, :]
            matrix1 = matrix[1, :]
            matrix2 = matrix[2, :]
            matrix3 = matrix[3, :]
            if isXChr and ind.sex == 0:
                matrixCollapsedHets = np.empty((2, matrix.shape[1]), dtype=np.float32)
                matrixCollapsedHets[0, :] = matrix0 + matrix2
                matrixCollapsedHets[1, :] = matrix1 + matrix3
            else:
                matrixCollapsedHets = np.empty((3, matrix.shape[1]), dtype=np.float32)
                matrixCollapsedHets[0, :] = matrix0
                matrixCollapsedHets[1, :] = matrix1 + matrix2
                matrixCollapsedHets[2, :] = matrix3

            calledGenotypes = np.argmax(matrixCollapsedHets, axis=0)
            setMissing(calledGenotypes, matrixCollapsedHets, thresh)
            f.write(ind.idx + " " + " ".join(map(str, calledGenotypes)) + "\n")


def writeCalledPhase(pedigree, genoProbFunc, isXChr, outputFile, thresh):
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
            matrix = genoProbFunc(ind.idn, ind.sex)
            matrix0 = matrix[0, :]
            matrix1 = matrix[1, :]
            matrix2 = matrix[2, :]
            matrix3 = matrix[3, :]

            # Paternal
            if isXChr and ind.sex == 0:
                paternal_haplotype = np.full(matrix.shape[1], 9, dtype=np.int8)
            else:
                paternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
                paternal_probs[0, :] = matrix0 + matrix1
                paternal_probs[1, :] = matrix2 + matrix3
                paternal_haplotype = np.argmax(paternal_probs, axis=0)
                setMissing(paternal_haplotype, paternal_probs, thresh)
            f.write(ind.idx + " " + " ".join(map(str, paternal_haplotype)) + "\n")

            # Maternal
            maternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
            maternal_probs[0, :] = matrix0 + matrix2
            maternal_probs[1, :] = matrix1 + matrix3
            maternal_haplotype = np.argmax(maternal_probs, axis=0)
            setMissing(maternal_haplotype, maternal_probs, thresh)
            f.write(ind.idx + " " + " ".join(map(str, maternal_haplotype)) + "\n")


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
