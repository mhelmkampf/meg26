### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 07. Genetic clustering                                                   ###
### ======================================================================== ###


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Establish SSH connection to cluster using your course account (please adjust user name)
ssh user1234@rosa.hpc.uni-oldenburg.de


### Either: Download course materials to course account using Git (if this is the first time)
git clone https://github.com/mhelmkampf/meg26.git


### Or: Update git repository (if it already exists)
cd meg26
git pull


### Start interactive shell session (new step!)
srun --pty --partition all_cpu.p --ntasks=1 --time=02:00:00 --mem=8G bash


### Load RStudio module
module load RStudio-Server


### Execute script to start RStudio
rstudio-start-on-rosa.sh


### Copy the SSH command provided at step 1) and execute in new terminal tab or window
#> ssh -N -L 8000: ...
#> re-enter your password, and note nothing will change in the terminal on success


### Go to http://localhost:8000 in your browser
#> RStudio will launch automatically



### ============================================================================
### Exercise 1: PCA based on SNP data

### Switch to the bash Terminal built into RStudio


### Create and navigate to today's working directory (e.g. "clust")
#>


module load VCFtools


### Create link to today's SNP dataset
ln -s /nfs/data/haex1482/shared/course/hamlets_snps_lg12.vcf.gz \
  hamlets_snps_lg12.vcf.gz


### Filter SNP in preparation for clustering by removing rare and physically linked sites
# key parameters: --maf, --thin
vcftools \
  --gzvcf hamlets_snps_lg12.vcf.gz \
  --minQ 40 \
  --minDP 10 \
  --max-missing 0.8 \
  --max-alleles 2 \
  --maf 0.05 \
  --thin 5000 \
  --recode \
  --stdout | bgzip > hamlets_filt_lg12.vcf.gz


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Switch to the R Console in RStudio

setwd("~/meg26/work/clust")

### Install / load packages
install.packages("vcfR")

library(vcfR)
library(adegenet)
library(tidyverse)


### Read VCF file into R
vcf <- read.vcfR("hamlets_filt_lg12.vcf.gz")


### Convert from vcfR to genlight object
gl <- vcfR2genlight(vcf)


### Principal Component Analysis (PCA)
pca <- glPca(gl, nf = 2)   # retain 2 principal components


### Print principal components to screen
#>


### Convert to tibble, add species information
scores <- pca$scores %>%
  as_tibble(rownames = "Sample") %>%
  mutate(Species = str_sub(Sample, -6, -4))


### Basic PCA plot
#> p <- ggplot(..., aes(x = ..., y = ..., color = ...)) + geom_point()


### Improve plot visually
q <- p + scale_color_manual(values = c("mediumseagreen",
                                       "orange",
                                       "royalblue",
                                       "gray30",
                                       "coral",
                                       "gray70")) +
  theme_light() +
  theme(
    text = element_text(color = "gray20"),
    panel.grid = element_blank(),
    axis.title.x = element_text(vjust = -1.5)
  )


### Eigenvalues: amount of genetic variance along each PC
pca$eig


### Proportion of variance explained, and its distribution
var <- pca$eig / sum(pca$eig)
barplot(var, main = "Proportion of variance explained", las = 2)


### Add proportion of variance explained to PCA axes
q + labs(
    x = paste0("PC1 (", round(pca$eig[1]/sum(pca$eig)*100, 1), "%)"),
    y = paste0("PC2 (", round(pca$eig[2]/sum(pca$eig)*100, 1), "%)")
  )



### ============================================================================
### Exercise 2: Estimate and plot ancestry proportions (ADMIXTURE)

### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Use the bash Terminal built into RStudio (bottom left)

module load BCFtools
module load PLINK
module load ADMIXTURE


### Note: For admixture analysis, sites are assumed to be independent.
# As an approximation for filtering by LD, we continue to use the dataset
# above filtered by physical distance (--thin 5000)


### Convert to BED format with plink (only accepts numbers as contig names)
gzip -cd hamlets_filt_lg12.vcf.gz | \
  sed 's/LG//g' | \
  plink --vcf /dev/stdin \
  --make-bed \
  --out hamlets_filt_lg12


### Run ADMIXTURE with cross-validation (CV)
for k in {1..6}
do
    admixture \
    --cv -j8 \
    hamlets_filt_lg12.bed $k > hamlets_filt_lg12_k${k}.out
done

# Note on CV:
# - for each k, temporarily hides a subset of genotype data
# - fits the model using the remaining data
# - predicts the hidden genotypes
# - calculates a CV error that reflects prediction accuracy (lower = better)


### Print CV error to find best k
for k in {1..6}
do
    grep 'CV' hamlets_filt_lg12_k${k}.out \
done


### What is the best k (i.e., that best explains the data)?
#>


### Retrieve and sort sample ids from plink output
awk '{ print $2 }' hamlets_filt_lg12.fam > sample_ids.txt


### Check whether sample_ids.txt contains all 36 samples
#>


### Add sample ids to ancestry proportions
for k in {1..6}
do
    paste -d " " sample_ids.txt hamlets_filt_lg12.${k}.Q \
    > anprop_k${k}.csv
done


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Read in data for k = 2
admix2 <- read_delim("anprop_k2.csv", delim = " ", col_names = c("Sample", "P1", "P2"))


### Pivot to long format and reformat sample name
long2 <- admix2 %>%
  pivot_longer(cols = starts_with("P"), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Name = paste0(str_sub(Sample, -6), "_", str_sub(Sample, 1, 5)))


### Basic admixture plot
a2 <- ggplot(long2, aes(x = Name, y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity")


### Improve plot visually
a2 + scale_fill_brewer(palette = "Dark2") +   # color palette
  labs(x = NULL, y = NULL, tag = "k = 2") +   # add tag label
  theme_minimal() +                           # basic style package
  theme(                                      # specific changes to style (fonts, grid lines etc.)
    text = element_text(color = "grey20"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, color = "gray20", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "gray20"),
    plot.tag = element_text(angle = -90),
    plot.tag.position = c(1.02, 0.6),
    legend.position = "none",
    plot.margin = unit(c(1, 10, 1, 1), "mm")
    )


### Read in and plot ancestry proportions for k = 3, 4, 5 or 6
#>


### How do the PCA and admixture plots compare?



### ============================================================================
### Optional

### Add labels of interesting samples to PCA plot
labels <- c("54761atlliz", "18267unibel")

(q + geom_text(aes(label = ifelse(Sample %in% labels, Sample, "")),
            size = 3, color = "gray20", vjust = -1)
)

### Zoom into main cluster
(q + coord_cartesian(xlim = c(-6, -3), ylim = c(-3, 3)))



### ============================================================================
### Solutions

### Create and navigate to today's working directory (e.g. "clust")
mkdir clust
cd clust


### Print principal components to screen
pca$scores


### Basic PCA plot
p <- ggplot(data = scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 4, alpha = 0.75)


### What is the best k?
# >>> CV error (K=2): 0.52125 <<<


### Read in and plot ancestry proportions for k = 3, 4, 5 or 6
admix4 <- read_delim("anprop_k4.csv", delim = " ", col_names = c("Sample", "P1", "P2", "P3", "P4"))


### Pivot to long format and reformat sample name
long4 <- admix4 %>%
  pivot_longer(cols = starts_with("P"), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Name = paste0(str_sub(Sample, -6), "_", str_sub(Sample, 1, 5)))


### Basic admixture plot
a4 <- ggplot(long4, aes(x = Name, y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity")


### Improve plot visually
a4 + scale_fill_brewer(palette = "Dark2") +   # color palette
labs(x = NULL, y = NULL, tag = "k = 4") +     # add tag label
  theme_minimal() +                           # basic style package
  theme(                                      # specific changes to style (fonts, grid lines etc.)
    text = element_text(color = "grey20"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, color = "gray20", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "gray20"),
    plot.tag = element_text(angle = -90),
    plot.tag.position = c(1.02, 0.6),
    legend.position = "none",
    plot.margin = unit(c(1, 10, 1, 1), "mm")
    )
