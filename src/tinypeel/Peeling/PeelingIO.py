import numpy as np
from numba import jit
from ..tinyhouse import InputOutput


def readInSeg(pedigree, fileName, start=None, stop = None):
    print("Reading in seg file:", fileName)
    if start is None: start = 0
    if stop is None: stop = pedigree.nLoci
    nLoci = stop - start + 1 #Contains stop.
    
    seg = np.full((pedigree.maxIdn, 4, nLoci), .25, dtype = np.float32)

    index = 0
    fileNColumns = 0

    indHit = np.full(pedigree.maxIdn, 0, dtype = np.int64)

    with open(fileName) as f:
        e = 0
        currentInd = None
        for line in f:
            parts = line.split();
            idx = parts[0]; 

            if fileNColumns == 0:
                fileNColumns = len(parts)
            if fileNColumns != len(parts):
                raise ValueError(f"The length of the line is not the expected length. Expected {fileNColumns} got {len(parts)} on individual {idx} and line {e}.")

            segLine=np.array([float(val) for val in parts[(start+1):(stop+2)]], dtype = np.float32)
            if len(segLine) != nLoci:
                raise ValueError(f"The length of the line subsection is not the expected length. Expected {nLoci} got {len(segLine)} on individual {idx} and line {e}.")

            if idx not in pedigree.individuals:
                print(f"Individual {idx} not found in pedigree. Individual ignored.")
            else:
                ind = pedigree.individuals[idx]
                if e == 0: 
                    currentInd = ind.idx
                if currentInd != ind.idx:
                    raise ValueError(f"Unexpected individual. Expecting individual {currentInd}, but got ind {ind.idx} on value {e}")
                seg[ind.idn,e,:] = segLine
                e = (e+1)%4
                ind.fileIndex['segregation'] = index; index += 1
                indHit[ind.idn] += 1
        for ind in pedigree:
            if indHit[ind.idn] != 4:
                print(f"No segregation information found for individual {ind.idx}")

    return seg

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

            if args.call_phase:
                writeCalledPhase(pedigree, genoProbFunc, args.out + ".called_phase." + str(thresh), thresh)

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

def writeCalledPhase(pedigree, genoProbFunc, outputFile, thresh):
    with open(outputFile, 'w+') as f:
        for idx, ind in pedigree.writeOrder():
            matrix = genoProbFunc(ind.idn)

            # Paternal            
            paternal_probs = np.array([matrix[0,:] + matrix[1,:], matrix[2,:] + matrix[3,:]], dtype=np.float32)
            paternal_haplotype = np.argmax(paternal_probs, axis = 0)
            setMissing(paternal_haplotype, paternal_probs, thresh)
            f.write(ind.idx + ' ' + ' '.join(map(str, paternal_haplotype)) + '\n')

            #Maternal
            maternal_probs = np.array([matrix[0,:] + matrix[2,:], matrix[1,:] + matrix[3,:]], dtype=np.float32)
            maternal_haplotype = np.argmax(maternal_probs, axis = 0)
            setMissing(maternal_haplotype, maternal_probs, thresh)
            f.write(ind.idx + ' ' + ' '.join(map(str, maternal_haplotype)) + '\n')

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

def fullOutput(pedigree, peelingInfo, args):
    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.penetrance, args.out + ".penetrance")
    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.anterior, args.out + ".anterior")
    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.posterior, args.out + ".posterior")

    InputOutput.writeFamIndexedMatrix(pedigree, peelingInfo.posteriorSire_minusFam, args.out + ".posteriorSire_minusFam")
    InputOutput.writeFamIndexedMatrix(pedigree, peelingInfo.posteriorDam_minusFam, args.out + ".posteriorDam_minusFam")

    InputOutput.writeFamIndexedMatrix(pedigree, peelingInfo.posteriorSire_new, args.out + ".posteriorSire_new")
    InputOutput.writeFamIndexedMatrix(pedigree, peelingInfo.posteriorDam_new, args.out + ".posteriorDam_new")
