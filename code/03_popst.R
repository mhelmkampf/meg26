### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 03. Populatioin structure and F-statistics                               ###
### ======================================================================== ###


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Establish SSH connection to cluster using your course account (please adjust user name)
ssh user1234@rosa.hpc.uni-oldenburg.de


### Either: Download course materials to course account using Git (if this is the first time)
git clone https://github.com/mhelmkampf/meg26.git


### Or: Update git repository (if it already exists)
cd meg26
git pull


### Load RStudio module
module load RStudio-Server


### Start interactive shell session
srun


### Execute script to start RStudio
rstudio-start-on-rosa.sh


### Copy the SSH command provided at 1) and execute in new terminal tab or window
#> ssh -N -L 8000: ...
#> re-enter your password, and note nothing will change in the terminal on success


### Go to http://localhost:8000 in your browser
#> RStudio will launch automatically



### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


### ============================================================================
### Exercise 1: Calculate F-statistics and test for population structure

### Load packages (installed last time, otherwise use install.packages())
library(tidyverse)
library(adegenet)
library(genepop)
library(pegas)


### Read in H. puella dataset (Caribbean populations)
caribbean <- read.genepop("data/msats/puella_caribbean.gen", ncode = 3)


### View data; how many populations are included?
# 


### Calculate F-statistics for each locus, across all populations
pegas::Fst(as.loci(caribbean))


### Calculate global F-statistics (averaged across loci)
genepop::Fst("data/msats/puella_caribbean.gen",
             outputFile = "work/Fst_caribbean.txt")


### According to Fst, is there population structure?


### For comparison, read in hamlet species dataset
hamletsp <- read.genepop("data/msats/hamlets_caribbean.gen", ncode = 3)


### Calculate F-statistics for the hamlet species dataset
#


### Test for genetic differentiation (G-test)
test_diff("data/msats/puella_caribbean.gen", outputFile = "work/Diff_caribbean.gen")
test_diff("data/msats/hamlets_caribbean.gen", outputFile = "work/Diff_hamlets.gen")


### Bonus question: Does self-fertilization occur regularly in hamlets?
# Hint: look at Fis in the H. puella dataset (caribbean)



### ============================================================================
### Exercise 2: Calculate population-specific Fst (beta)

### Load packages
install.packages("hierfstat")
library(hierfstat)


### Population-specific Fst (beta)
# How much does each population deviate from the overall mean in terms of allele frequencies?
betas(hamletsp)


### Convert output to tidyverse data frame (tibble) -- execute step by step to follow the pipeline
b <- betas(hamletsp)$betaiovl %>%                    # extract Fsts from betas object
  bind_rows() %>%                                    # convert to tibble
  pivot_longer(cols = everything(),                  # transform data from wide to long format
               names_to = "Population",
               values_to = "Fst") %>%
  mutate(Species = str_sub(Population, 1, 3)) %>%    # add colum with species id
  arrange(desc(Fst))                                 # sort data by descending Fst


### Basic plotting with ggplot
ggplot(data = b, aes(x = Population, y = Fst)) +     # determine basic data structure (mapping)
  geom_bar(stat = "identity")                        # determine type of plot (geom)


### Make plot prettier
ggplot(data = b, 
       aes(x = reorder(Population, -Fst),           # reorder x-axis by descending Fst values
           y = Fst,
           fill = Species)) +                       # add fill mapping (color bars by species)
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("abe" = "#DFDF8D",   # define fill colors for species
                               "chl" = "#8B4513", 
                               "gum" = "#F99729",
                               "ind" = "#22198E", 
                               "nig" = "#333333", 
                               "pue" = "#E48175", 
                               "uni" = "#B3B3B3")) +
  labs(x = "Population") +                          # set x-axis label
  theme_minimal(base_size = 16) +                   # use theme "minimal" (e.g. no frame, white background)
  theme(
    panel.grid.minor = element_blank(),             # adjust theme: remove minor grid lines
    panel.grid.major.x = element_blank(),           # adjust theme: remove major grid lines intersecting x-axis
    axis.text.x = element_text(angle = 45,          # adjust theme: change angle and position of x-axis labels
                               hjust = 1, 
                               vjust = 1.25)    
  )



### ============================================================================
### Exercise 3: Pairwise Fst / Visualize population structure using PCoA

### Calculate matrix of pairwise Fst
d <- genet.dist(hamlets, method = "Nei87")   # pairwise Fst following Nei (1978)
d


### Basic plotting
p <- pcoa(as.matrix(d))


### Create tidyverse data frame (tibble) with first two axes
t <- tibble(pco1 = p$vecp[, 1], 
            pco2 = p$vecp[, 2], 
            population = colnames(as.matrix(d))) %>%
  mutate(Species = str_sub(population, 1, 3))


### Plot (and store plot in variable)
g <-
  ggplot(data = t,
         aes(x = pco1, 
             y = pco2,
             color = Species)) +
  geom_point(size = 3) +
  scale_fill_manual(values = c("abe" = "#DFDF8D",
                               "chl" = "#8B4513", 
                               "gum" = "#F99729",
                               "ind" = "#22198E", 
                               "nig" = "#333333", 
                               "pue" = "#E48175", 
                               "uni" = "#B3B3B3")) +
  labs(x = "PCoA 1",
       y = "PCoA 2") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "transparent")
  )

g


### How does the Fst-based PCoA compare to the Fst barplot from exercise 2?



### ============================================================================
### Exercise 1 solution

### View data, check population slot
caribbean
caribbean@pop
levels(caribbean@pop)


### Calculate F-statistics for the hamlet species dataset
pegas::Fst(as.loci(hamletsp))

genepop::Fst("data/msats/hamlets_caribbean.gen",
             outputFile = "work/Fst_hamlets.txt")


### Bonus question: Does self-fertilization occur regularly in hamlets?
# Hint: look at Fis in the H. puella dataset (caribbean)
pegas::Fst(as.loci(caribbean))

# Self-fertilization can be considered an extreme form of inbreeding, 
# reducing the number of heterozygotes in each generation. The inbreeding coefficient 
# Fis is close to 0 in H. puella, so self-fertilization does not seem to occur regularly
