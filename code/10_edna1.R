### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 10. Metabarcoding / eDNA I                                                        ###
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



### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Note: the following steps have already been executed to prepare the 
### data tables for this exercise:

# - Removing MiFish primer sequences from Fastq raw reads using cutadapt
# - Removing low-quality tails and reads that are too short or too long

# - ASV identification using DADA2 (simplified excerpt)

# library(dada2)

# Estimate error rates
# err_fwd <- learnErrors(files$fwd_filt,
#                        errorEstimationFunction = loessErrfun)

# Dereplicate reads
# derep_fwd <- derepFastq(files$fwd_filt)

# Denoise reads
# dada_fwd <- dada(derep_fwd, 
#                  err = err_fwd)

# - Removing chimeric sequences using removeBimeraDenovo()
# - Removing PCR duplicates with LULU



### ============================================================================
### Exercise 1: Inspect sample and ASV count tables


### Set working directory
setwd("~/meg26")


### Load packages
library(tidyverse)


### Read in sample metadata table
sample_table <- read_tsv("data/edna/sample_metadata.tsv")


### Summarize number of samples per variable (Site, Habitat, Season)
sample_table %>%
  count(Site, Habitat, Season)


### Read in and inspect ASV count table
asv_table <- read_tsv("data/edna/asv_table.tsv")
view(asv_table)


### How many ASVs are there?
#>


### Read in key linking ASV IDs to sequences
asv_key <- read_tsv("data/edna/asv_key.tsv")



### ============================================================================
### Exercise 2: Taxonomic assignment using BLAST

### Pick an ASV of your choice from asv_table and BLAST its sequence found in asv_key
# https://blast.ncbi.nlm.nih.gov/Blast.cgi


### According to FishBase, is this result ecologically and geographically plausible?
# https://www.fishbase.se/search.php


### Read in and inspect BLAST taxonomy table
# Note: this was produced by BLAST against Midori2 database
# https://www.reference-midori.info
blast_table <- read_tsv("data/edna/taxonomy_table.tsv")
view(blast_table)



### ============================================================================
### Exercise 3: Construct phyloseq object and filter

install.packages("phyloseq")
library(phyloseq)


### Convert ASV count table to suitable format
asv_mat <- asv_table %>%
  column_to_rownames("sample") %>%
  as.matrix()

asv_ps <- otu_table(asv_mat, taxa_are_rows = FALSE)


### Convert BLAST taxonomy table to suitable format
tax_ps <- blast_table %>%
  filter(qseqid %in% colnames(asv_table)) %>%
  select(qseqid, kingdom, phylum, class, order, family, genus, species) %>%
  column_to_rownames("qseqid") %>%
  as.matrix() %>%
  tax_table()


### Convert sample metadata table to suitable format
sample_ps <- sample_table %>%
  column_to_rownames("Sample") %>%
  sample_data()


### Combine both with sample metadata table to create phyloseq object
ps <- phyloseq(
  OTU = asv_ps,
  TaxonData = tax_ps,
  SampleData = sample_ps
)


### Count total number of ASVs per sample (post-processing sequencing depth)
sample_reads <- tibble(
  Sample = sample_names(ps),
  total_reads = sample_sums(ps)
) %>%
  arrange(total_reads)


### Plot sequencing depth per sample
ggplot(sample_reads, aes(x = fct_reorder(Sample, total_reads), y = total_reads)) +
  geom_col() +
  geom_hline(yintercept = 1000, color = "red") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Total reads"
  ) +
  theme_minimal()


### Are there outlier samples?


### Filter unreliable samples / ASVs
# - Remove samples with < 1000 reads
# - Remove ASVs with < 20 total counts
# - Remove non-fish ASVs (Actinopteri, ray-finned fishes))
ps_filt <- ps %>%
  prune_samples(sample_sums(.) >= 1000, .) %>%
  filter_taxa(function(x) sum(x) >= 20, prune = TRUE) %>%
  subset_taxa(class == "Actinopteri")


### For the sake of this exercise, we will ignore additional important filters:
# - Remove ASVs that are likely contaminants (e.g., from negative controls)
# - Remove ASVs that are likely mis-identified (e.g., from BLAST results)
# - Remove ASVs that are present in only one replicate per sample (i.e., singletons)


### Save filtered phyloseq object for use next time
saveRDS(ps_filt, "work/ps_filt.rds")



### ============================================================================
### Exercise 4: Plot ASV abundance distribution

### Calculate total abundance per ASV, summed across all samples
asv_totals <- tibble(
  ASV = taxa_names(ps_filt),
  total_abundance = taxa_sums(ps_filt)
)


### Plot distribution
ggplot(asv_totals, aes(x = total_abundance)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  labs(
    x = "Total abundance (log scale)",
    y = "Number of ASVs",
    title = "Distribution of ASV abundance across all samples"
  ) +
  theme_minimal()


### Determine the 20 most abundant ASVs
top20 <- asv_totals %>% slice_head(n = 20)


### Plot the 20 most abundant ASVs
#>


### What is the most abundant species in the dataset?



### ============================================================================
### Solutions

ggplot(top20, aes(x = fct_reorder(ASV, total_abundance), y = total_abundance)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "ASV",
    y = "Total abundance",
    title = "Top 20 most abundant ASVs"
  )
