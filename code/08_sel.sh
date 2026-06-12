### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 08. Selection                                                            ###
### ======================================================================== ###


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Establish SSH connection to cluster using your course account (please adjust user name)
ssh user1234@rosa.hpc.uni-oldenburg.de


### Update git repository
cd meg26
git pull


### ALTERNATIVELY: If you receive an error, delete and re-download the repository
# cd
# rm -rf meg26
# git clone https://github.com/mhelmkampf/meg26.git


### Start interactive shell session
srun --pty --partition all_cpu.p --ntasks=1 --time=02:00:00 --mem=8G bash


### Load RStudio module
module load RStudio-Server


### Execute script to start RStudio
rstudio-start-on-rosa.sh


### Copy the SSH command provided at step 1) and execute in new terminal tab or window
#> ssh -N -L 8000: ...
#> re-enter your password, and note nothing will change in the terminal on success


### Go to http://localhost:8000 in your browser (ideally, Firefox or Safari)
#> RStudio will launch automatically



### ============================================================================
### Exercise 1: Perform genome-wide Fst scan

### Switch to the bash Terminal built into RStudio


### Create and navigate to today's working directory (e.g. "sel")
#>


### Copy population id files from ~/meg26/data/meta to working directory
#>


### Create link to today's SNP dataset
ln -s /nfs/data/haex1482/shared/course/hamlets_snps_lg12.vcf.gz \
  hamlets_snps_lg12.vcf.gz


### Load required software
module load VCFtools


### Filter SNP in preparation for Fst scan
# key parameters: --maf, --remove, no thinning
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --minQ 40 \
  --minDP 10 \
  --max-missing 0.8 \
  --max-alleles 2 \
  --maf 0.05 \
  --remove pop_atlhua.txt \
  --recode \
  --stdout | bgzip > hamlets_sel_lg12.vcf.gz


### How many sites are retained in the new SNP dataset?
#>


### Calculate joint Fst in sliding windows of 50 kb
# Note: represents overall among-population variance relative to total variance
vcftools \
  --gzvcf hamlets_sel_lg12.vcf.gz \
  --weir-fst-pop pop_gumhon.txt \
  --weir-fst-pop pop_indbel.txt \
  --weir-fst-pop pop_nigbel.txt \
  --weir-fst-pop pop_puebel.txt \
  --weir-fst-pop pop_unibel.txt \
  --fst-window-step 5000 \
  --fst-window-size 50000 \
  --stdout 1> Fst_lg12_50k.tsv 2> Fst_lg12_50k.log


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Switch to the R Console in RStudio


### Set working directory
setwd("~/meg26/work/sel")


### Load packages
library(tidyverse)


### Plot joint Fst, along 50 kb windows (Mantattan plot)
fst_50k <- read_tsv("Fst_lg12_50k.tsv")


### Create Manhattan plot by plotting BIN_START vs. WEIGHTED_FST using ggplot's geom_point()
#> f <- ggplot(..., aes(x = ..., y = ...) + geom_point()


### Assign chromosome-wide, weighted Fst estimate from Fst_lg12_50k.log to new variable (using bash)
#> chrwide_fst <- ...


### Add line representing chromosome-wide Fst to plot
f + geom_hline(yintercept = chrwide_fst, color = "blue")


### Identify 99 % quantile (= 1 % most differentiated sites)
threshold <- quantile(fst_50k$WEIGHTED_FST, 0.99)


### Add line representing 99 % quantile to plot
#> f + geom_hline(yintercept = chrwide_fst, color = "blue") + geom_hline(...)


### Retrieve positions of windows containing Fst peaks
outliers <- fst_50k %>% 
  filter(WEIGHTED_FST >= threshold)


### What are the boundaries of the genomic region around the highest Fst peak?
#> 



### ============================================================================
### Exercise 2: Investigate candidate genes under selection


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Create link to the genome annotation file in GFF format
ln -s /nfs/data/haex1482/shared/course/HypPue1_annotation.gff \
  HypPue1_annotation.gff


### Extract region of interest from GFF file
awk '$1 == "LG12" && $5 >= 20135001 && $4 <= 20295000' HypPue1_annotation.gff \
  > highfst_region.gff


### What genes are found in this region?
awk '$3 == "gene"' highfst_region.gff


### Find out more by searching for "casz1" at
#- Zebrafish Information Network (closest model organism): https://zfin.org
#- NCBI Gene: https://www.ncbi.nlm.nih.gov/gene
#- Uniprot: https://www.uniprot.org



### ============================================================================
### Optional: Extended Haplotype Homozygosity-based test using rehh


### Create link to phased SNP dataset (minor allele count 2, thinned to 2 kb, without H. atlahua)
ln -s /nfs/data/haex1482/shared/course/uni_phased2k_lg12.vcf.gz \
  uni_phased2k_lg12.vcf.gz


### Check for phasing
zless uni_phased2k_lg12.vcf.gz
# Note: phased genotypes use "|" instead of "/"


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Install / Load packages
install.packages("rehh")

library(rehh)
library(vcfR)


### Read in VCF with phased haplotypes (creates haplohh object)
hap_uni <- data2haplohh(hap_file = "uni_phased2k_lg12.vcf.gz",
                        haplotype.in.columns = TRUE,
                        polarize_vcf = TRUE,
                        vcf_reader = "vcfR")


### Scan for EHH (Extended Haplotype Homology) along chromosome
scan_uni <- scan_hh(hap_uni)


### Calculate integrated haplotype score(iHS, standardized across allele frequencies)
ihs_uni <- ihh2ihs(scan_uni)


### Manhattan plot
manhattanplot(ihs_uni)


# Make plot prettier: Remove NAs and filter to iHS results
ihs_df <- ihs_uni$ihs %>%
  filter(!is.na(IHS))

ggplot(ihs_df, aes(x = POSITION, y = IHS)) +
  geom_point(aes(color = abs(IHS) > 2), size = 1.2) +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") +
  labs(title = "iHS along Chromosome (Unicolor)",
       x = "Genomic Position (bp)",
       y = "iHS") +
  theme_minimal()



### ============================================================================
### Solutions

### Create and navigate to today's working directory (e.g. "sel")
cd work
mkdir sel
cd sel


### Copy population id files from ~/meg26/data/meta to working directory
cp ~/meg26/data/meta/*.txt .


### How many sites are retained in the new SNP dataset?
vcftools --gzvcf hamlets_sel_lg12.vcf.gz


### Create Manhattan plot by plotting BIN_START vs. WEIGHTED_FST using ggplot's geom_point()
f <- ggplot(fst_50k, aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_point(size = 0.25, alpha = 0.5) +
  labs(x = "Position", y = "Fst") +
  theme_minimal() +
  theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
        )


### Assign chromosome-wide, weighted Fst estimate from Fst_lg12_50k.log to new variable (using bash)
cat Fst_lg12_50k.log
chrwide_fst <- 0.058942


### Add line representing 99 % quantile to plot
f + geom_hline(yintercept = chrwide_fst, color = "blue") 
  + geom_hline(yintercept = threshold, color = "red")


### What are the boundaries of the genomic region around the highest Fst peak?
#> 20135001 - 20295000 (160 kb)
