mkdir outputs
# AlphaPeel is a command line package for performing multilocus iterative peeling using genotype or sequence data.
# Install TinyPeel via pip using:
# pip install <wheel name>

#To check command line arguments run TinyPeel without any arguments.
TinyPeel

# Example 1: Performing multi-locus peeling with genotype data:

# TinyPeel -genotypes data/genotypes.txt \
python ../src/tinyPeel-script.py -genotypes data/genotypes.txt \
         -pedigree data/pedigree.txt \
         -out outputs/multilocus \
         -nCycles 5 \
         -runType multi \
         -maxthreads 6

# Example 2: Performing single-locus "hybrid" peeling with sequence data and pre-computed segregation estimates (generated from Example 1).

python ../src/tinyPeel-script.py -seqfile data/sequence.txt \
         -pedigree data/pedigree.txt \
         -mapfile data/genotypes-map.txt\
         -out outputs/hybrid \
         -runType single \
         -segmapfile data/segregation-map.txt \
         -segfile outputs/multilocus.seg \
         -nCycles 5 \
         -maxthreads 6
