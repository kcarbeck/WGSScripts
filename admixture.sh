## ADMIXTURE is a clustering software similar to STRUCTURE with the aim to infer populations and individual ancestries.
## april 5, 2021

## First we have to prune for linkage. ADMIXTURE requires unlinked SNPs in plink format.
# adadapted from https://speciationgenomics.github.io
  # --double-id: tells plink to duplicate the id of our samples
  # --allow-extra-chr: allows additional chromosomes beyond the human chromosome set that plink expects by default
  # --set-missing-var-ids: sets variant ID for SNPs. Human and model organisms often have annotated SNP names and so plink will look for these. I did not have them so I set default to chromosome:position allele 1, allele2 using the option @:#\$1:\$2
  # --indep-pairwise: for linkage pruning. The first argument, 50 denotes we have set a window of 50 Kb. The second argument, 10 is our window step size - meaning we move 10 bp each time we calculate linkage. Finally, we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate. Here we prune any variables that show an r2 of greater than 0.5.
  # output files: prune.in = sites above threshold (i.e., those we should retain); prune.out = sites below threshold
/programs/plink-1.9-x86_64-beta3.30/plink --vcf FILENAME.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:#\$1:\$2 \
--indep-pairwise 50 10 0.5 --out FILENAME

## Prune sites and make plink file
  # --extract: extract only these positions from our VCF
  # --make-bed: make plink file
/programs/plink-1.9-x86_64-beta3.30/plink --vcf FILENAME.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#\$1:\$2 \
--extract FILENAME.prune.in --make-bed --out FILENAME_LDpruned_plink


### ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1=0;print $0}' FILENAME_LDpruned_plink.bim > FILENAME_LDpruned_plink.bim.tmp
mv FILENAME_LDpruned_plink.bim.tmp FILENAME_LDpruned_plink.bim


## run ADMIXTURE in loop
# default is 5-fold cross validation
# loop for K=8 (we sampled from 8 different locations)
for i in {1..8}
do
 /programs/admixture/admixture --cv=5 FILENAME_LDpruned_plink.bed $i > log${i}.out
done


## identify the best value of k clusters (lowest cross-validation error)
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > admixture.cv.error
