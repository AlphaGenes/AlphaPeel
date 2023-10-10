#module load roslin/plink/1.90p

tinypeel=AlphaPeel


# Test 1: Can we read in unrelated individuals from multiple file formats and output the values to a normal dosage file.
rm -rf test1/outputs
mkdir   test1/outputs
$tinypeel -genotypes test1/genotypes.txt \
                          -phasefile test1/phasefile.txt \
                          -penetrance test1/penetrance.txt \
                          -seqfile test1/seqfile.txt \
                          -pedigree test1/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 \
                          -esterrors \
                          -out test1/outputs/output
Rscript checkResults.r 1

# Test 2a: Can we read in a subset of values as in Test 1 output them and make sure it's the same chunk?
rm -rf test2/outputs
mkdir   test2/outputs


$tinypeel -genotypes test2/genotypes.txt \
                          -phasefile test2/phasefile.txt \
                          -penetrance test2/penetrance.txt \
                          -seqfile test2/seqfile.txt \
                          -pedigree test2/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 \
                          -startsnp 2 \
                          -stopsnp 4 \
                          -out test2/outputs/output
Rscript checkResults.r 2


# Test 3: Can we read in values, call the values and output them as binary?

rm -rf test3/outputs
mkdir test3/outputs
$tinypeel -genotypes test3/genotypes.txt \
                          -phasefile test3/phasefile.txt \
                          -penetrance test3/penetrance.txt \
                          -seqfile test3/seqfile.txt \
                          -pedigree test3/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 0.99\
                          -startsnp 2 \
                          -stopsnp 4 \
                          -binary_call_files \
                          -out test3/outputs/output

plink --bfile test3/outputs/output.called.0.1 --real-ref-alleles --recode A --out test3/outputs/output.called.0.1
plink --bfile test3/outputs/output.called.0.99 --real-ref-alleles --recode A --out test3/outputs/output.called.0.99

rm -rf test3b/outputs
mkdir test3b/outputs

# Test 3b: Can we read in a binary file, run the algorithm call the values and check that the output is the same.

for nind in 1 2 3 4; do

  $tinypeel -genotypes test3b/genotypes-$nind.txt \
                            -runType multi \
                            -calling_threshold .1 \
                            -binary_call_files \
                            -out test3b/outputs/output.$nind
  plink --bfile test3b/outputs/output.$nind.called.0.1 --real-ref-alleles --recode A --out test3b/outputs/output.$nind.called.0.1

  $tinypeel -bfile test3b/outputs/output.$nind.called.0.1 \
                            -runType multi \
                            -calling_threshold .1 \
                            -binary_call_files \
                            -out test3b/outputs/round2.$nind

  plink --bfile test3b/outputs/round2.$nind.called.0.1 --real-ref-alleles --recode A --out test3b/outputs/round2.$nind.called.0.1

done

Rscript checkResults.r 3b

# Test 3c: Will the pedigree file be correctly read from the bed file?

#Create the binary file. Re-run.

rm -rf test3c/outputs
mkdir test3c/outputs

$tinypeel -genotypes test3c/genotypes.txt \
                          -runType multi \
                          -calling_threshold .9 \
                          -binary_call_files \
                          -out test3c/outputs/output
plink --bfile test3c/outputs/output.called.0.9 --real-ref-alleles --recode A --out test3c/outputs/output.called.0.9

$tinypeel -bfile test3c/outputs/output.called.0.9 \
                          -runType multi \
                          -calling_threshold .9 \
                          -out test3c/outputs/noFamNoPedigree

$tinypeel -bfile test3c/outputs/output.called.0.9 \
                          -pedigree test3c/pedigree.txt \
                          -runType multi \
                          -calling_threshold .9 \
                          -out test3c/outputs/noFamPedigree

cp test3c/fake.fam test3c/outputs/output.called.0.9.fam
$tinypeel -bfile test3c/outputs/output.called.0.9 \
                          -runType multi \
                          -calling_threshold .9 \
                          -out test3c/outputs/famNoPedigree

Rscript checkResults.r 3c

# Test 4a-d: Can we read in values and return them in the correct order. Check id, pedigree, genotypes, sequence, segregation. Also check onlykeyed.
rm -rf test4/outputs
mkdir   test4/outputs

for method in id pedigree genotypes sequence; do
    $tinypeel -genotypes test4/genotypes.txt \
                              -phasefile test4/phasefile.txt \
                              -penetrance test4/penetrance.txt \
                              -seqfile test4/seqfile.txt \
                              -pedigree test4/pedigree.txt \
                              -runType multi \
                              -calling_threshold .1 \
                              -out test4/outputs/output.$method \
                              -writekey $method
done
$tinypeel -genotypes test4/genotypes.txt \
                          -phasefile test4/phasefile.txt \
                          -penetrance test4/penetrance.txt \
                          -seqfile test4/seqfile.txt \
                          -pedigree test4/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 \
                          -out test4/outputs/output.only \
                          -writekey sequence \
                          -onlykeyed

Rscript checkResults.r 4

# Test 5a: Read in an error rate and output genotypes. Should use some silly values, and some not-so-silly values.
rm -rf test5/outputs
mkdir   test5/outputs

$tinypeel -genotypes test5/genotypes.txt \
                          -phasefile test5/phasefile.txt \
                          -penetrance test5/penetrance.txt \
                          -seqfile test5/seqfile.txt \
                          -pedigree test5/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 \
                          -seqerror 0.5 \
                          -error 1.0 \
                          -nophasefounders \
                          -haps \
                          -out test5/outputs/chance
$tinypeel -genotypes test5/genotypes.txt \
                          -phasefile test5/phasefile.txt \
                          -penetrance test5/penetrance.txt \
                          -seqfile test5/seqfile.txt \
                          -pedigree test5/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 \
                          -seqerror 0.001 \
                          -error 0.01 \
                          -haps \
                          -out test5/outputs/normal
Rscript checkResults.r 5




# Test 6: Sex Chromosome

rm -rf test6/outputs
mkdir   test6/outputs

$tinypeel -genotypes test6/genotypes.txt \
                          -seqfile test6/seqfile.txt \
                          -pedigree test6/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 \
                          -sexchrom \
                          -out test6/outputs/output
Rscript checkResults.r 6


# Test 7: Check -esterrors just to make sure it runs. 


rm -rf test7/outputs
mkdir   test7/outputs

$tinypeel -genotypes test7/genotypes.txt \
                          -phasefile test7/phasefile.txt \
                          -penetrance test7/penetrance.txt \
                          -seqfile test7/seqfile.txt \
                          -pedigree test7/pedigree.txt \
                          -runType multi \
                          -esterrors \
                          -calling_threshold .1 \
                          -out test7/outputs/output
Rscript checkResults.r 7

# Test 7b: Check -estmaf just to make sure it runs. 

rm -rf test7b/outputs
mkdir   test7b/outputs

$tinypeel -genotypes test7b/genotypes.txt \
                          -phasefile test7b/phasefile.txt \
                          -penetrance test7b/penetrance.txt \
                          -seqfile test7b/seqfile.txt \
                          -pedigree test7b/pedigree.txt \
                          -runType multi \
                          -estmaf \
                          -calling_threshold .1 \
                          -out test7b/outputs/output
Rscript checkResults.r 7b

# Test 7c: Check -estmaf just to make sure it runs. 

rm -rf test7c/outputs
mkdir   test7c/outputs

$tinypeel -genotypes test7c/genotypes.txt \
                          -phasefile test7c/phasefile.txt \
                          -penetrance test7c/penetrance.txt \
                          -seqfile test7c/seqfile.txt \
                          -pedigree test7c/pedigree.txt \
                          -runType multi \
                          -length 1.0 \
                          -calling_threshold .1 \
                          -out test7c/outputs/output
Rscript checkResults.r 7c

# Test 8: Check to make sure the no_dosage, no_seg, no_params flags work, and the haps file works.

rm -rf test8/outputs
mkdir   test8/outputs

$tinypeel -genotypes test8/genotypes.txt \
                          -phasefile test8/phasefile.txt \
                          -penetrance test8/penetrance.txt \
                          -seqfile test8/seqfile.txt \
                          -pedigree test8/pedigree.txt \
                          -runType multi \
                          -no_dosage \
                          -out test8/outputs/no_dosage
$tinypeel -genotypes test8/genotypes.txt \
                          -phasefile test8/phasefile.txt \
                          -penetrance test8/penetrance.txt \
                          -seqfile test8/seqfile.txt \
                          -pedigree test8/pedigree.txt \
                          -runType multi \
                          -no_seg \
                          -out test8/outputs/no_seg
$tinypeel -genotypes test8/genotypes.txt \
                          -phasefile test8/phasefile.txt \
                          -penetrance test8/penetrance.txt \
                          -seqfile test8/seqfile.txt \
                          -pedigree test8/pedigree.txt \
                          -runType multi \
                          -no_params \
                          -out test8/outputs/no_params
$tinypeel -genotypes test8/genotypes.txt \
                          -phasefile test8/phasefile.txt \
                          -penetrance test8/penetrance.txt \
                          -seqfile test8/seqfile.txt \
                          -pedigree test8/pedigree.txt \
                          -runType multi \
                          -haps \
                          -out test8/outputs/haps


Rscript checkResults.r 8





Rscript checkResults.r 1
Rscript checkResults.r 2
Rscript checkResults.r 3
Rscript checkResults.r 3b
Rscript checkResults.r 3c

Rscript checkResults.r 4
Rscript checkResults.r 5

Rscript checkResults.r 7
Rscript checkResults.r 7b
Rscript checkResults.r 7c

Rscript checkResults.r 8
