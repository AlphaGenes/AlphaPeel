rm outputs/*

##AlphaPeel

python ../src/tinyPeel.py -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6

python ../tinyPeel.py -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.estErrorsAndMaf \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6 \
                        -esterrors \
                        -estmaf

#Multi locus with sequence
python ../tinyPeel.py -seqfile baseData/sequence.2 \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.multi.seq \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6 

#Single locus
python ../tinyPeel.py -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.single \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -estmaf


##AlphaPeel with sequence

Rscript extractSeg.r

python ../tinyPeel.py -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.hybrid \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -mapfile baseData/map.txt\
                        -segmapfile baseData/segmap.txt\
                        -segfile outputs/seg.subset.txt

python ../tinyPeel.py -seqfile baseData/sequence.2 \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.hybrid.seq \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -mapfile baseData/map.txt\
                        -segmapfile baseData/segmap.txt\
                        -segfile outputs/seg.subset.txt


###AlphaAssign

python ../tinyAssign.py -genotypes baseData/trueGenotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/assign \
                        -potentialsires baseData/sire.list

python ../tinyAssign.py -genotypes baseData/trueGenotypes.txt \
                        -pedigree baseData/pedigree.noDam \
                        -out outputs/assign.noDam \
                        -potentialsires baseData/sire.list

python ../tinyAssign.py -seqfile baseData/sequence.2 \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/assign.seq \
                        -potentialsires baseData/sire.list

python ../tinyAssign.py -genotypes baseData/trueGenotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/assign.check \
                        -checkpedigree

python ../tinyAssign.py -genotypes baseData/trueGenotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/assign.check.opp \
                        -checkpedigree \
                        -runtype opp 



#AlphaMGSAssign

python ../tinyMgsAssign.py -genotypes baseData/trueGenotypes.txt \
                        -pedigree baseData/pedigree.noDam \
                        -out outputs/mgsAssign \
                        -potentialgrandsires baseData/grandsire.list


DATE=`date +%Y-%m-%d`

Rscript checkResultsFromSummary.r results.2018-10-05 results.$DATE