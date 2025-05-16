import numpy as np
from numba import jit
from ..tinyhouse import InputOutput


def readInSeg(pedigree, fileName, start=None, stop=None):
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
                print(f"Individual {idx} is not found in pedigree. Individual ignored.")
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
    args = InputOutput.args

    np.savetxt(args.out_file + ".geno_error_prob.txt", peelingInfo.genoError, fmt="%f")
    np.savetxt(args.out_file + ".seq_error_prob.txt", peelingInfo.seqError, fmt="%f")
    np.savetxt(
        args.out_file + ".rec_prob.txt", np.empty((1, 1)), fmt="%f"
    )  # not be realized, just as a placeholder
    # np.savetxt(args.out_file + ".trans", peelingInfo.transmissionRate, fmt = "%f")
    # np.savetxt(args.out_file + ".alt_allele_prob.txt", peelingInfo.maf, fmt="%f")


def writeOutAltAlleleProb(pedigree):
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


def writeGenotypes(pedigree, genoProbFunc, isXChr):
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
            if args.binary_call_file:
                writeBinaryCalledGenotypes(
                    pedigree,
                    genoProbFunc,
                    isXChr,
                    args.out_file + ".called." + str(threshold),
                    threshold,
                )
            else:
                writeCalledGenotypes(
                    pedigree,
                    genoProbFunc,
                    isXChr,
                    args.out_file + ".geno_" + str(threshold) + ".txt",
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
            if args.binary_call_file:
                pass  # this function is not applied
            else:
                writeCalledPhase(
                    pedigree,
                    genoProbFunc,
                    args.out_file + ".hap_" + str(threshold) + ".txt",
                    threshold,
                )


def writePhasedGenoProbs(pedigree, genoProbFunc, outputFile):
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)
            f.write("\n")
            for i in range(matrix.shape[0]):
                f.write(
                    ind.idx + " " + " ".join(map("{:.4f}".format, matrix[i, :])) + "\n"
                )


def writeGenoProbs(pedigree, genoProbFunc, outputFile):
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


def writeDosages(pedigree, genoProbFunc, isXChr, outputFile):
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            if isXChr and ind.sex == 0:
                tmp = np.array([0, 1, 0, 1])
            else:
                tmp = np.array([0, 1, 1, 2])
            matrix = np.dot(tmp, genoProbFunc(ind.idn))
            f.write(ind.idx + " " + " ".join(map("{:.4f}".format, matrix)) + "\n")


def writeCalledGenotypes(pedigree, genoProbFunc, isXChr, outputFile, thresh):
    with open(outputFile, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)
            if isXChr and ind.sex == 0:
                matrixCollapsedHets = np.array(
                    [matrix[0, :] + matrix[2, :], matrix[1, :] + matrix[3, :]],
                    dtype=np.float32,
                )
            else:
                matrixCollapsedHets = np.array(
                    [matrix[0, :], matrix[1, :] + matrix[2, :], matrix[3, :]],
                    dtype=np.float32,
                )

            calledGenotypes = np.argmax(matrixCollapsedHets, axis=0)
            setMissing(calledGenotypes, matrixCollapsedHets, thresh)
            f.write(ind.idx + " " + " ".join(map(str, calledGenotypes)) + "\n")


def writeCalledPhase(pedigree, genoProbFunc, outputFile, thresh):
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


def writeBinaryCalledGenotypes(pedigree, genoProbFunc, isXChr, outputFile, thresh):
    for idx, ind in pedigree.writeOrder():
        matrix = genoProbFunc(ind.idn)
        matrixCollapsedHets = np.array(
            [matrix[0, :], matrix[1, :] + matrix[2, :], matrix[3, :]], dtype=np.float32
        )
        calledGenotypes = np.argmax(matrixCollapsedHets, axis=0)
        setMissing(calledGenotypes, matrixCollapsedHets, thresh)
        ind.genotypes = calledGenotypes.astype(np.int8)

    InputOutput.writeOutGenotypesPlink(pedigree, outputFile)


@jit(nopython=True)
def setMissing(calledGenotypes, matrix, thresh):
    nLoci = len(calledGenotypes)
    for i in range(nLoci):
        if matrix[calledGenotypes[i], i] <= thresh:
            calledGenotypes[i] = 9


def fullOutput(pedigree, peelingInfo, args):
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
