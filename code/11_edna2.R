### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 11. Metabarcoding / eDNA II                                              ###
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

### Set working directory
setwd("~/meg26")


### Load packages
library(tidyverse)
library(phyloseq)


### ============================================================================
### Exercise 1: Alpha diversity


### Load filtered phyloseq object from last time
ps <- readRDS("data/edna/ps.rds")


### How many samples and ASVs (taxa) are in the dataset?
#>


### View the metadata table
sample_data(ps)


### Filter, for now without removing rare ASVs (important for alpha diversity)
ps_filt <- ps %>%
  subset_samples(Habitat != "Negative") %>%
  prune_samples(sample_sums(.) >= 1000, .) %>%
  subset_taxa(class == "Actinopteri")


### Create grouping variable from every metadata column except Replicate
sample_group <- data.frame(sample_data(ps_filt)) %>%
  mutate(sample_group = paste(Site, Habitat, Season, sep = "_")) %>%
  pull(sample_group)


### Merge replicates by summing ASV counts per group
ps_merged <- merge_samples(ps_filt, group = sample_group)
otu_table(ps_merged) <- otu_table(otu_table(ps_merged), 
                                  taxa_are_rows = taxa_are_rows(ps_filt))


### Rebuild metadata table for merged samples
meta_df <- data.frame(sample_data(ps_filt)) %>%
  remove_rownames() %>%
  mutate(sample_group = sample_group) %>%
  distinct(sample_group, .keep_all = TRUE) %>%
  select(-Replicate) %>%
  column_to_rownames("sample_group")
sample_data(ps_merged) <- meta_df[sample_names(ps_merged), ]


### Rarefy to an even depth (default = minimum observed depth across merged samples)
min(sample_sums(ps_merged))

set.seed(1)   # random seed for reproducibility
ps_rare <- rarefy_even_depth(ps_merged, rngseed = 1)


### Aggregate ASV counts to species level
ps_species <- tax_glom(ps_rare, taxrank = "species")


### Calculate observed species richness and Shannon diversity
alpha_div <- estimate_richness(ps_species, measures = c("Observed", "Shannon"))


### Attach sample metadata for plotting
alpha_df <- alpha_div %>%
  bind_cols(data.frame(sample_data(ps_species)))


### Plot alpha diversity (observed species richness) by habitat
# Note: use geom_boxplot() and optionally geom_jitter() to show individual points
#>



### Save work / results
saveRDS(ps_species, "work/ps_species.rds")
write_tsv(alpha_df, "results/edna_alpha.tsv")
ggsave("results/edna_alpha.png", width = 5, height = 5, dpi = 300)



### ============================================================================
### Exercise 2: Beta diversity


### Filter species: remove those present in only 1 sample, or with <20 reads total
ps_beta <- ps_species %>%
  filter_taxa(function(x) sum(x > 0) >= 2, TRUE) %>%
  filter_taxa(function(x) sum(x) >= 20, TRUE)


### How many ASVs were removed?


### Calculate Bray-Curtis dissimilarity matrix
bray_dist <- phyloseq::distance(ps_beta, method = "bray")


### Ordinate using PCoA on Bray-Curtis dissimilarity
bray_ord <- ordinate(ps_beta, method = "PCoA", distance = bray_dist)


#### Extract PCoA axis scores and attach sample metadata
beta_df <- bray_ord$vectors %>%
  data.frame() %>%
  rownames_to_column("sample_group") %>%
  left_join(
    data.frame(sample_data(ps_beta)) %>% rownames_to_column("sample_group"),
    by = "sample_group"
  )


### Percent variance explained, for axis labels
eig <- bray_ord$values$Relative_eig
pct1 <- round(eig[1] * 100, 1)
pct2 <- round(eig[2] * 100, 1)


### Plot beta diversity PCoA
# Note: use geom_point()
#>


### Use color in aes() to color points by Habitat, Site, and Season



### Save work / results
write_tsv(beta_df, "results/edna_beta.tsv")
saveRDS(bray_dist, "results/edna_bc-matrix.rds")
ggsave("results/edna_beta.png", width = 5, height = 4, dpi = 300)


### PERMANOVA: tests for significant variance in community composition
# (i.e. whether group centroids are further apart than expected by chance)

library(vegan)

group_meta <- data.frame(sample_data(ps_beta))

adonis2(bray_dist ~ Season, 
        data = group_meta, permutations = 10000, by = "margin")


### Test the other variables as well
#>


### Test for homogeneity of dispersion (variance) among groups
# Note: if significant, may indicate PERMANOVA result is due to differences in variance 
# rather than community composition
betadisper(bray_dist, group_meta$Season) |> permutest()
betadisper(bray_dist, group_meta$Habitat) |> permutest()
betadisper(bray_dist, group_meta$Site) |> permutest()


###



### ============================================================================
### Exercise 3 (optional): Community composition

### Aggregate to family level (rarefied ASV table, same as species-level step)
tax_df <- data.frame(tax_table(ps_rare)) %>%
  mutate(family = replace_na(family, "Unassigned"))

tax_table(ps_rare) <- as.matrix(tax_df)

ps_family <- tax_glom(ps_rare, taxrank = "family", NArm = FALSE)

view(tax_table(ps_family))


### Convert to long format for ggplot (one row per sample x family; similar to pivot_longer())
family_df <- psmelt(ps_family)


### Identify top 10 families by total abundance across the whole dataset
top10_families <- family_df %>%
  group_by(family) %>%
  summarise(total_abund = sum(Abundance)) %>%
  slice_max(total_abund, n = 10) %>%
  pull(family)


### Collapse everything else into "Other"
family_df <- family_df %>%
  mutate(family = if_else(family %in% top10_families, family, "Other"))


### Plot community composition by family
ggplot(family_df, aes(x = Sample, y = Abundance, fill = family)) +
  geom_col() +
  facet_wrap(~ Habitat, scales = "free_x") +
  labs(x = NULL, y = "Read count", fill = "Family") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


### Save plot
ggsave("results/edna_comp.png", width = 10, height = 6, dpi = 300)



### ============================================================================
### Solutions

### How many samples and ASVs (taxa) are in the dataset?
ps


### Optional: Define colors for plotting
habitat_colours <- c(
  Seagrass = "mediumseagreen",
  Sand     = "#c9a84c",
  Channel  = "coral"
)


### Plot alpha diversity (observed species richness) by habitat
ggplot(alpha_df, aes(x = Habitat, y = Observed, fill = Habitat)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  scale_fill_manual(values = habitat_colours) +
  labs(y = "Species richness (Observed)", x = "Habitat") +
  theme_bw()


### How many ASVs were removed?
394 - 247


### Plot beta diversity PCoA (replace "Habitat" with "Site" or "Season" to color by those variables)
ggplot(beta_df, aes(x = Axis.1, y = Axis.2, color = Habitat)) +
  geom_point(size = 3) +
  # scale_color_manual(values = habitat_colours) +
  labs(x = paste0("PCoA1 (", pct1, "%)"),
       y = paste0("PCoA2 (", pct2, "%)")) +
  theme_bw()


### Test the other variables as well
adonis2(bray_dist ~ Habitat, 
        data = group_meta, permutations = 10000, by = "margin")

adonis2(bray_dist ~ Site, 
        data = group_meta, permutations = 10000, by = "margin")


# as_tibble(data.frame(sample_data(ps_filt)), rownames = "SampleID")
