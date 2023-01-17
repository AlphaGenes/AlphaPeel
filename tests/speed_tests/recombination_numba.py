from numba import jit
import numpy as np

@jit(nopython=True, nogil=True)
def estimateRecombinations(pointSeg, forwardSeg, backwardSeg):
    nLoci = pointSeg.shape[1]
    recomb = np.full(nLoci -1, 0, dtype = np.float32)
    recomb_mat = np.full(nLoci -1, 0, dtype = np.float32)
    recomb_pat = np.full(nLoci -1, 0, dtype = np.float32)
    
    e = 0.01
    score = np.array([[0, 1, 1, 2],
                      [1, 0, 2, 1],
                      [1, 2, 0, 1],
                      [2, 1, 1, 0]])
    
    
    score_pat = np.array([[0, 0, 1, 1],
                          [0, 0, 1, 1],
                          [1, 1, 0, 0],
                          [1, 1, 0, 0]])

    score_mat = np.array([[0, 1, 0, 1],
                          [1, 0, 1, 0],
                          [0, 1, 0, 1],
                          [1, 0, 1, 0]])

    base = np.array([[(1-e)**2,  (1-e)*e,    (1-e)*e,    e**2],
                    [(1-e)*e,   (1-e)**2,   e**2,       (1-e)*e],
                    [(1-e)*e,   e**2,       (1-e)**2,   (1-e)*e],
                    [e**2,      (1-e)*e,    (1-e)*e,    (1-e)**2]])


    for i in range(nLoci -1):
        recomb[i], recomb_mat[i], recomb_pat[i] = estimateLocusRecombination_for_loops(pointSeg, forwardSeg, backwardSeg, 0.01, i, score, score_pat, score_mat, base)

    return np.sum(recomb)/2

@jit(nopython=True, nogil=True)
def estimateLocusRecombination(pointSeg, forwardSeg, backwardSeg, transmissionRate, locus, score, score_pat, score_mat, base):
    # Estimates the transmission rate between locus and locus + 1.
    val_current = pointSeg[:,locus]*forwardSeg[:,locus]
    val_current = val_current/np.sum(val_current)

    val_next = pointSeg[:,locus+1]*backwardSeg[:,locus+1]
    val_next = val_next/np.sum(val_next)
    
    # norm1D(val_current)
    # norm1D(val_next)

    e = 0.01

    mat = np.full((4, 4), 0, dtype = np.float32)
    mat[:,:] = base[:,:]

    # Now create joint probabilities.
    for i in range(4):
        for j in range(4):
            mat[i,j] *= val_current[i]*val_next[j]

    mat = mat/np.sum(mat)
    error = np.sum(mat*score)
    error_mat = np.sum(mat*score_mat)
    error_pat = np.sum(mat*score_pat)
    # count = 0
    # for i in range(4):
    #     for j in range(4):
    #         count += mat[i,j]
    # for i in range(4):
    #     for j in rnage(4):
    #         mat[i,j]/=count
    # norm2D(mat)
    # error = 0
    # for i in range(4):
    #     for j in range(4):
    #         error += mat[i,j]*score[i,j]
    return(error, error_mat, error_pat)

@jit(nopython=True, nogil=True)
def estimateLocusRecombination_for_loops(pointSeg, forwardSeg, backwardSeg, transmissionRate, locus, score, score_pat, score_mat, base):
    # Estimates the transmission rate between locus and locus + 1.
    val_current = np.full(4, 0, dtype = np.float32)
    val_next = np.full(4, 0, dtype = np.float32)

    count = 0
    for i in range(4):
        val_current[i] = pointSeg[i,locus]*forwardSeg[i,locus]
    for i in range(4):
        count += val_current[i]
    for i in range(4):
        val_current[i]/=count

    count = 0
    for i in range(4):
        val_next[i] = pointSeg[i,locus+1]*backwardSeg[i,locus+1]
    for i in range(4):
        count += val_next[i]
    for i in range(4):
        val_next[i]/=count


    # norm1D(val_current)
    # norm1D(val_next)

    e = 0.01

    mat = np.full((4, 4), 0, dtype = np.float32)
    for i in range(4):
        for j in range(4):
            mat[i,j] = base[i,j]

    # Now create joint probabilities.
    for i in range(4):
        for j in range(4):
            mat[i,j] *= val_current[i]*val_next[j]


    # mat = mat/np.sum(mat)
    # count = 0
    # for i in range(4):
    #     for j in range(4):
    #         count += mat[i,j]
    # for i in range(4):
    #     for j in range(4):
    #         mat[i,j]/=count
    norm2D(mat)
    # error = np.sum(mat*score)
    # error_mat = np.sum(mat*score_mat)
    # error_pat = np.sum(mat*score_pat)

    error = 0
    error_mat = 0
    error_pat = 0
    for i in range(4):
        for j in range(4):
            error += mat[i,j]*score[i,j]
            error_mat += mat[i,j]*score_pat[i,j]
            error_pat += mat[i,j]*score_mat[i,j]
    return(error, error_mat, error_pat)

@jit(nopython=True, nogil=True)
def norm1D(values) :
    count = 0
    for i in range(4):
        count += values[i]
    for i in range(4):
        values[i] /= count

@jit(nopython=True, nogil=True)
def norm2D(values) :
    count = 0
    for i in range(4):
        for j in range(4):
            count += values[i, j]
    for i in range(4):
        for j in range(4):
            values[i,j] /= count



pointSeg = np.random.random((20000, 4, 1000)).astype(np.float32)
forwardSeg = np.random.random((20000, 4, 1000)).astype(np.float32)
backwardSeg = np.random.random((20000, 4, 1000)).astype(np.float32)


estimateRecombinations(pointSeg[0,:,:], forwardSeg[0,:,:], backwardSeg[0,:,:])
import time

start = time.time()
for i in range(20000):
    estimateRecombinations(pointSeg[i,:,:], forwardSeg[i,:,:], backwardSeg[i,:,:])

end = time.time()
print(end - start)






