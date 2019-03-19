rm outputs/*

tinypeel=TinyPeel
##AlphaPeel

$tinypeel -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6

$tinypeel -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.estErrorsAndMaf \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6 \
                        -esterrors \
                        -estmaf

#Multi locus with sequence
$tinypeel -seqfile baseData/sequence.2 \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.multi.seq \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6 

#Single locus
$tinypeel -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.single \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -estmaf


##AlphaPeel with sequence

Rscript extractSeg.r

$tinypeel -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.hybrid \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -mapfile baseData/map.txt\
                        -segmapfile baseData/segmap.txt\
                        -segfile outputs/seg.subset.txt

$tinypeel -seqfile baseData/sequence.2 \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.hybrid.seq \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -mapfile baseData/map.txt\
                        -segmapfile baseData/segmap.txt\
                        -segfile outputs/seg.subset.txt

Rscript checkResults.r