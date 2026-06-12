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
### Exercise 1: Genome-wide Fst scan

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


### Plot BIN_START vs. WEIGHTED_FST using ggplot's geom_point()
#>

f <- ggplot(fst_50k, aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_point(size = 0.25, alpha = 0.5) +
  labs(x = "Position", y = "Fst") +
  theme_minimal() +
  theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
        )


### Find the chromosome-wide, weighted Fst estimate in the VCFtools log file (using bash)
#>


### Assign value to a variable called "chrwide_fst"
#> 


### Add chromosome-wide Fst to plot
f + geom_hline(yintercept = chrwide_fst, color = "blue")


### Identify 99.5% quantile
quantile(fst_50k$WEIGHTED_FST, probs = 0.995)


### Add variable to data table indicating outlier status
out <- fst_50k %>%
  mutate(OUTLIER = case_when(
    WEIGHTED_FST > quantile(WEIGHTED_FST, probs = 0.995) ~ "yes",
    TRUE ~ "no")
  )


### Plot with 99.5% quantile highlighted
o <- ggplot(data = out, aes(x = BIN_START, y = WEIGHTED_FST, color = OUTLIER)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = chrwide_fst, color = "blue") +
    labs(x = "Window", y = "Fst") +
    scale_color_manual(values = c("gray20", "red")) +
    guides(color = "none") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()
          )


### Retrieve positions of windows containing Fst peaks
print(out %>% filter(OUTLIER == "yes"), n = nrow(out))


### Genomic region of interest around highest Fst peak
#> 20140001 - 20290000 (150 kb)


### ----------------------------------------------------------------------------
### The following inactive code was used to plot a PCA of the region of interest

### Extract genomic region of interest from VCF file (bash)
# vcftools --gzvcf hamlets_sel_lg12.vcf.gz \
# --chr LG12 \
# --from-bp 20140001 \
# --to-bp 20290000 \
# --recode \
# --stdout | gzip > fstpeak_region.vcf.gz


### Load packages
# library(vcfR)
# library(adegenet)


### Read VCF file into R
# vcf <- read.vcfR("fstpeak_region.vcf.gz")


### Convert from vcfR to genlight object
# data <- vcfR2genlight(vcf)


### Principal Component Analysis (PCA)
# pca <- glPca(data, nf = 2)


### Convert to tibble, add species information
# scores <- as.data.frame(pca$scores) %>%
#   rownames_to_column("Sample") %>%
#   as_tibble() %>%
#   mutate(Species = str_sub(Sample, -6, -4))


### Plot PCA
# (p <- ggplot(data = scores, aes(x = PC1, y = PC2, color = Species)) +
#   geom_point(size = 5, alpha = 0.75) +
#   scale_color_manual(values = c("orange", "royalblue", "gray30", "coral", "gray70")) +
#   theme_light() +
#   theme(
#     text = element_text(color = "gray20"),
#     panel.grid = element_blank(),
#     axis.title.x = element_text(vjust = -1.5)
#   )
# )


### How do the PCAs of the high Fst region and the whole chromosome differ?
### What trait seems to be affected, pointing to genes related to that trait under selection in the region?




### ============================================================================
### Exercise 2: Identify candidate genes under selection


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Create link to the genome annotation file in GFF format
ln -s /nfs/data/haex1482/shared/course/HypPue1_annotation.gff \
  HypPue1_annotation.gff


###


### Extract region of interest from GFF file
awk '$1 == "LG12" && $5 >= 20140001 && $4 <= 20290000' HypPue1_annotation.gff \
  > highfst_region.gff


### What genes are found in this region?


### Find out more by searching for "casz1" at
### NCBI Gene: https://www.ncbi.nlm.nih.gov/gene
### Uniprot: https://www.uniprot.org



### ============================================================================
### Solutions

### Create and navigate to today's working directory (e.g. "sel")
cd work
mkdir sel
cd sel


### Copy population id files from ~/meg26/data/meta to working directory
cp ~/meg26/data/meta/*.txt .