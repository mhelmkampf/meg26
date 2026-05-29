### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 06. Population genomics and genetic diversity                            ###
### ======================================================================== ###


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Establish SSH connection to cluster using your course account (please adjust user name)
ssh -Y user1234@rosa.hpc.uni-oldenburg.de   # -Y enables graphical windows


### Update git repository
cd meg26
git pull


### ALTERNATIVELY: If you receive an error, delete and re-download the repository
# cd
# rm -rf meg26
# git clone https://github.com/mhelmkampf/meg26.git


### Create and navigate to today's working directory (e.g. "gdiv")
#>



### ============================================================================
### Exercise 1: Heterozygosity and inbreeding coefficient at variant sites

module load BCFtools
module load VCFtools


### Create link to today's SNP dataset
ln -s /nfs/data/haex1482/shared/course/hamlets_snps_lg12.vcf.gz \
  hamlets_snps_lg12.vcf.gz


### Filter SNPs without dropping rare alleles (--maf)
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --minQ 40 \
  --minDP 10 \
  --max-missing 0.8 \
  --max-alleles 2 \
  --recode \
  --stdout | bgzip > hamlets_bi_lg12.vcf.gz


### Calculate heterozygosity and Fis per individual
vcftools \
  --gzvcf hamlets_bi_lg12.vcf.gz \
  --het \
  --stdout > Het_hamlets.tsv


### What is unusual about these Fis values, and why?
#>


### Copy population id files to working directory
cp ../../data/meta/*.txt .


### Re-calculate heterozygosity for puebel sub-population only
vcftools \
  --gzvcf hamlets_bi_lg12.vcf.gz \
  --keep hamlet_puebel.txt \
  --het \
  --stdout > Het_puebel.tsv


### Re-calculate heterozygosity for atlahua sub-population (hamlet_atlhua.txt) only
#>



### ============================================================================
### Exercise 2: Assess nucleotide diversity per population

### To calculate nucleotide diversity, we need variant AND invariant sites
ln -s /nfs/data/haex1482/shared/course/hamlets_all_1M.vcf.gz hamlets_all_1M.vcf.gz


### Measure nucleotide diversity per population using shell loop
for i in hamlet_*.txt; do
vcftools \
  --gzvcf hamlets_all_1M.vcf.gz \
  --keep "$i" \
  --window-pi 50000 \
  --stdout > Pi_"$i".tsv
done


### Average pi per population
for i in Pi_*.tsv; do
  echo $i
  awk '$5 != "-nan" { sum += $5 ; next } END { print sum / NR }' $i
  echo
done



### ============================================================================
### Exercise 3: Filter by linkage disequilibrium

### Referring to the linkage decay plot, how many basepairs (bp) apart should the 
### sites be from each other to avoid linkage disequilibrium?


### Remove rare alleles for LD analyses
vcftools --gzvcf hamlets_bi_lg12.vcf.gz \
  --maf 0.05 \
  --max-missing 0.9 \
  --recode \
  --stdout | bgzip > hamlets_maf_lg12.vcf.gz


### Calculate R^2 without thinning
vcftools --gzvcf hamlets_maf_lg12.vcf.gz \
  --geno-r2 \
  --ld-window-bp 10000 \
  --out ld_unthinned


### Thin to 1 SNP per 5000 bp
vcftools --gzvcf hamlets_maf_lg12.vcf.gz \
  --thin 5000 \
  --recode \
  --stdout | bgzip > hamlets_thin_lg12.vcf.gz


### Calculate R^2 after thinning
vcftools --gzvcf hamlets_thin_lg12.vcf.gz \
  --geno-r2 \
  --ld-window-bp 10000 \
  --out ld_thinned


### Calculate mean R^2
# Unthinned
awk 'NR>1 && $5!="-nan" {sum+=$5; n++} END {print "mean r2 =", sum/n}' \
  ld_thinned.geno.ld
  
awk 'NR>1 && $5!="-nan" {sum+=$5; n++} END {print "mean r2 =", sum/n}' \
  ld_unthinned.geno.ld



### ============================================================================
### Solutions:

### Re-calculate heterozygosity for atlahua sub-population only
vcftools \
  --gzvcf hamlets_bi_lg12.vcf.gz \
  --keep hamlet_atlhua.txt \
  --het \
  --stdout > Het_atlhua.tsv


### ============================================================================
### Alternative code and notes

### Summarize in ASCII
for file in ld_thinned.geno.ld ld_unthinned.geno.ld; do
  echo
  echo "$file:"
  awk '$5 != "-nan" { sum += $5; n++; vals[n] = $5 }
  END {
    if (n == 0) exit
    asort(vals)
    median = (n % 2) ? vals[int((n+1)/2)] : (vals[n/2] + vals[n/2 + 1]) / 2
    print "  Mean R²: " sum/n
    print "  Median R²: " median
    print "  Count: " n
  }' "$file"
done
