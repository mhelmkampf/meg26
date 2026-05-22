### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 05. Variant calling and SNPs                                             ###
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


### Create and navigate to today's working directory (e.g. "snp")
#>



### ============================================================================
### Exercise 1: Understanding VCF files

ml VCFtools   # https://vcftools.github.io/man_latest.html
ml BCFtools   # https://samtools.github.io/bcftools/bcftools.html


### Create link to today's SNP dataset
ln -s /nfs/data/haex1482/shared/course/hamlets_snps_lg12.vcf.gz hamlets_snps_lg12.vcf.gz


### Inspect the header (## / #)
zcat hamlets_snps_lg12.vcf.gz | grep "^#" | less


### Inspect the data
zcat hamlets_snps_lg12.vcf.gz | grep  -v "^#" | less
# -v inverts grep's search pattern


# Point out heterozygous and homozygous genotypes


### How many variants are there?
#> 


### Read VCF file with VCFtools
vcftools --gzvcf hamlets_snps_lg12.vcf.gz



### ============================================================================
### Exercise 2: Quality summaries

### Calculate missing data per individual
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --missing-indv \
  --out indiv_missing

less indiv_missing.imiss | column -t


### Calculate depth per individual with --depth
#>

less depth.idepth | column -t


# Are missingness and depth correlated?
ml R

Rscript -e '
  miss  <- read.table("indiv_missing.imiss", header=TRUE)
  depth <- read.table("depth.idepth", header=TRUE)
  d <- merge(miss, depth, by="INDV")

  png("indv_depth_vs_missing.png", width=600, height=500)
  plot(d$MEAN_DEPTH, d$F_MISS,
       xlab="Mean depth", ylab="Fraction missing",
       pch=16, col="red")
  dev.off()
'

display indv_depth_vs_missing.png


### Calculate missing data per site
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --missing-site \
  --out site_missing

head site_missing.lmiss | column -t



### ============================================================================
### Exercise 3: SNP filtering pipeline

### Step 1: remove low-quality sites (Q: probability of site being non-variant)
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --minQ 40 \
  --out 1_qfilt


### Step 2: add minimum depth and missing data filter
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --minQ 40 \
  --minDP 10 \
  --max-missing 0.8 \
  --out 2_mfilt

# Note: genotypes that fail --minDP filter are set to missing,
# sites are dropped if > 20 % missing with --max-missing


### Step 3: keep only biallelic sites with --max-alleles
#>


### Step 4: remove rare variants (MAC filter)
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --minQ 40 \
  --minDP 10 \
  --max-missing 0.8 \
  --max-alleles 2 \
  --mac 2 \
  --out 4_final

# Why might we want to remove rare alleles?


### Write final VCF
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --minQ 40 \
  --minDP 10 \
  --max-missing 0.8 \
  --max-alleles 2 \
  --mac 2 \
  --recode \
  --stdout | bgzip > hamlets_filt_lg12.vcf.gz


### Summary statistics of final VCF
bcftools stats hamlets_filt_lg12.vcf.gz | grep "^SN"


### Calculate allele counts before filtering
vcftools --gzvcf hamlets_snps_lg12.vcf.gz \
  --counts \
  --out raw_snps

less raw_snps.frq.count


### Calculate allele counts after filtering
vcftools --gzvcf hamlets_filt_lg12.vcf.gz \
  --counts \
  --out filt_snps

less filt_snps.frq.count


### Extract minor allele frequency with awk
awk 'NR > 1 {
  split($5, a, ":")
  split($6, b, ":")
  mac = (a[2] < b[2]) ? a[2] : b[2]
  maf = mac / $4
  print maf
}' A_freq.frq > maf.txt


### Plot in R
Rscript -e '
  maf <- scan("maf.txt")
  
  png("maf_distribution.png", width=800, height=500)
  hist(maf,
       breaks=50,
       col="steelblue",
       border="white",
       main="Minor allele frequency distribution",
       xlab="Minor allele frequency (MAF)",
       ylab="Number of sites")
  abline(v=0.05, col="red", lty=2, lwd=2)
  legend("topright", legend="MAF = 0.05 threshold", col="red", lty=2, lwd=2)
  dev.off()
  
  cat("Total sites:", length(maf), "\n")
  cat("Sites with MAF < 0.05:", sum(maf < 0.05), "\n")
  cat("Sites with MAF == 0:", sum(maf == 0), "\n")
'

display maf_distribution.png



### ============================================================================
### Solutions:

### How many variants are there?
zcat hamlets_snps_lg12.vcf.gz | grep  -v "^#" | wc -l


### Calculate depth per individual with --depth
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --depth \
  --out depth


### Step 3: keep only biallelic sites with --max-alleles
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --minQ 40 \
  --minDP 10 \
  --max-missing 0.8 \
  --max-alleles 2 \
  --out 3_bfilt
