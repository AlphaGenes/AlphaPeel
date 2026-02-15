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

# Example 1: Performing multi-locus peeling with genotype data:
AlphaPeel -geno_file data/genotypes.txt \
          -ped_file data/pedigree.txt \
          -out_file outputs/multilocus \
          -n_cycle 5 \
          -method multi \
          -n_thread 6 \
          -seg_prob

# Example 1b: Performing multi-locus peeling with genotype data and calling the values with a threshold of 0.98
AlphaPeel -geno_file data/genotypes.txt \
          -ped_file data/pedigree.txt \
          -out_file outputs/multilocus_with_phase \
          -n_cycle 5 \
          -method multi \
          -n_thread 6 \
          -geno_threshold 0.98 \
          -hap_threshold 0.98 \
          -geno \
          -hap

# Example 2: Performing single-locus "hybrid" peeling with sequence data and pre-computed segregation estimates (generated from Example 1).
AlphaPeel -seq_file data/sequence.txt \
         -ped_file data/pedigree.txt \
         -map_file data/genotypes-map.txt\
         -out_file outputs/hybrid \
         -method single \
         -seg_map_file data/segregation-map.txt \
         -seg_file outputs/multilocus.seg_prob.txt \
         -n_cycle 5 \
         -n_thread 6
