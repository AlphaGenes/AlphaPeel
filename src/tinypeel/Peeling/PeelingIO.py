import numpy as np
from numba import jit
from ..tinyhouse import InputOutput

def writeOutParamaters(peelingInfo) :
    args = InputOutput.args

    np.savetxt(args.out + ".genoError", peelingInfo.genoError, fmt = "%f")
    np.savetxt(args.out + ".seqError", peelingInfo.seqError, fmt = "%f")
    # np.savetxt(args.out + ".trans", peelingInfo.transmissionRate, fmt = "%f")
    np.savetxt(args.out + ".maf", peelingInfo.maf, fmt = "%f")

def writeGenotypes(pedigree, genoProbFunc) :
    args = InputOutput.args
    if not args.no_dosages: writeDosages(pedigree, genoProbFunc, args.out + ".dosages")
    if args.haps: writeGenoProbs(pedigree, genoProbFunc, args.out + ".haps")

    if args.calling_threshold is not None:
        for thresh in args.calling_threshold:
            if args.binary_call_files : writeBinaryCalledGenotypes(pedigree, genoProbFunc, args.out + ".called." + str(thresh), thresh)
            if not args.binary_call_files : writeCalledGenotypes(pedigree, genoProbFunc, args.out + ".called." + str(thresh), thresh)


def writeGenoProbs(pedigree, genoProbFunc, outputFile):
    with open(outputFile, 'w+') as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)
            f.write('\n')
            for i in range(matrix.shape[0]) :
                f.write(ind.idx + ' ' + ' '.join(map("{:.4f}".format, matrix[i,:])) + '\n')

def writeDosages(pedigree, genoProbFunc, outputFile):
    with open(outputFile, 'w+') as f:
        for idx, ind in pedigree.writeOrder():
            matrix = np.dot(np.array([0,1,1,2]), genoProbFunc(ind.idn))
            
            if InputOutput.args.sexchrom and ind.sex == 0:
                matrix *= 2

            f.write(ind.idx + ' ' + ' '.join(map("{:.4f}".format, matrix)) + '\n')


def writeCalledGenotypes(pedigree, genoProbFunc, outputFile, thresh):
    with open(outputFile, 'w+') as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)
            
            matrixCollapsedHets = np.array([matrix[0,:], matrix[1,:] + matrix[2,:], matrix[3,:]], dtype=np.float32)
            calledGenotypes = np.argmax(matrixCollapsedHets, axis = 0)
            setMissing(calledGenotypes, matrixCollapsedHets, thresh)
            if InputOutput.args.sexchrom and ind.sex == 0:
                doubleIfNotMissing(calledGenotypes)

            f.write(ind.idx + ' ' + ' '.join(map(str, calledGenotypes)) + '\n')

def writeBinaryCalledGenotypes(pedigree, genoProbFunc, outputFile, thresh):
    for idx, ind in pedigree.writeOrder():
        matrix = genoProbFunc(ind.idn)
        matrixCollapsedHets = np.array([matrix[0,:], matrix[1,:] + matrix[2,:], matrix[3,:]], dtype=np.float32)
        calledGenotypes = np.argmax(matrixCollapsedHets, axis = 0)
        setMissing(calledGenotypes, matrixCollapsedHets, thresh)
        if InputOutput.args.sexchrom and ind.sex == 0:
            doubleIfNotMissing(calledGenotypes)
        ind.genotypes = calledGenotypes.astype(np.int8)

    InputOutput.writeOutGenotypesPlink(pedigree, outputFile)

@jit(nopython=True)
def doubleIfNotMissing(calledGenotypes):
    nLoci = len(calledGenotypes)
    for i in range(nLoci):
        if calledGenotypes[i] == 1:
            calledGenotypes[i] = 2

@jit(nopython=True)
def setMissing(calledGenotypes, matrix, thresh) :
    nLoci = len(calledGenotypes)
    for i in range(nLoci):
        if matrix[calledGenotypes[i],i] < thresh:
            calledGenotypes[i] = 9

