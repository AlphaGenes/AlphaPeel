#!/bin/bash

mkdir -p outputs

# To check command line arguments run AlphaPeel without any arguments.
AlphaPeel

# Example 1: Performing multi-locus peeling with genotype data:
AlphaPeel -genotypes data/genotypes.txt \
         -pedigree data/pedigree.txt \
         -out outputs/multilocus \
         -nCycles 5 \
         -runType multi \
         -maxthreads 6

# Example 1b: Performing multi-locus peeling with genotype data and calling the values with a threshold of 0.98
AlphaPeel -genotypes data/genotypes.txt \
         -pedigree data/pedigree.txt \
         -out outputs/multilocus_with_phase \
         -nCycles 5 \
         -runType multi \
         -maxthreads 6 \
         -calling_threshold 0.98 \
         -hap

# Example 2: Performing single-locus "hybrid" peeling with sequence data and pre-computed segregation estimates (generated from Example 1).
AlphaPeel -seqfile data/sequence.txt \
         -pedigree data/pedigree.txt \
         -mapfile data/genotypes-map.txt\
         -out outputs/hybrid \
         -runType single \
         -segmapfile data/segregation-map.txt \
         -segfile outputs/multilocus.seg \
         -nCycles 5 \
         -maxthreads 6
