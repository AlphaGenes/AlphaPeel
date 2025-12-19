#!/bin/bash

# AlphaPeel is a command line package for imputation in pedigree populations
# Install AlphaPeel via pip using:
# pip install AlphaPeel

# To see command line arguments run AlphaPeel without any arguments or providing -h or --help
# AlphaPeel
# AlphaPeel -h
# AlphaPeel --help

# Below is a set of examples that give you a flavour on how to run and use AlphaPeel
mkdir -p outputs

# Example 1: Run the multi-locus peeling with genotype data
AlphaPeel -genotypes data/genotypes.txt \
          -pedigree data/pedigree.txt \
          -out outputs/multilocus \
          -nCycles 5 \
          -runType multi \
          -maxthreads 6

# Example 1b: Run the multi-locus peeling with genotype data and calling the values with high-confidence
AlphaPeel -genotypes data/genotypes.txt \
          -pedigree data/pedigree.txt \
          -out outputs/multilocus_with_phase \
          -nCycles 5 \
          -runType multi \
          -maxthreads 6 \
          -geno_threshold 0.98 \
          -hap_threshold 0.98 \
          -geno \
          -hap

# Example 2: Run the single-locus "hybrid" peeling with sequence data and pre-computed segregation estimates (generated from Example 1)
AlphaPeel -seqfile data/sequence.txt \
          -pedigree data/pedigree.txt \
          -mapfile data/genotypes-map.txt\
          -out outputs/hybrid \
          -runType single \
          -segmapfile data/segregation-map.txt \
          -segfile outputs/multilocus.seg \
          -nCycles 5 \
          -maxthreads 6
