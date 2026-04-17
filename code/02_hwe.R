### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 01. Introduction and setup                                               ###
### ======================================================================== ###


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Establish SSH connection to cluster using your course account (please adjust user name)
ssh user1234@rosa.hpc.uni-oldenburg.de


### Either: Download course materials to course account using Git (if this is the first time)
git clone https://github.com/mhelmkampf/meg26.git


### Or: Update git repository (if it already exists)
cd meg25
git pull


### Load RStudio module
module load RStudio-Server


### Execute script to start RStudio
rstudio-start-on-rosa.sh


### Copy the SSH command provided at 1) and execute in new terminal tab or window
#> ssh -N -L 8000: ...
#> re-enter your password, and note nothing will change in the terminal on success


### Go to http://localhost:8000 in your browser
#> RStudio will lanuch automatidcally



### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### ============================================================================
### Exercise 1: Manual test for Hardy-Weinberg equilibrium


### What are the observed genotype frequencies?

# Assign value to variable (here, size of the population)
N <- 12

# Recall variable
N

# Calculate the observed frequency of the AA genotype
AA_obscount <- 6
AA_obsfreq <- AA_obscount / N
AA_obsfreq

# Calculate the frequency of the remaining genotypes
# BB_obscount <- 
# BB_obsfreq <- 

# AB_obscount <- 
# AB_obsfreq <- 
  
  
### What are the allele frequencies (p for allele A, q for allele B)?
p <- (AA_obscount * 2 + AB_obscount) / (N * 2)
# q <- 


### Is the population in Hardy-Weinberg equilibrium?

# Calculate expected genotype frequencies under HWE
AA_expfreq <- p ^ 2
# AB_expfreq <- 
# BB_expfreq <- 

# Load tidyverse package required for next step
library(tidyverse)

# Enter data into a table ("tibble") and add expected genotype counts
genotype_data <- tibble(
  genotype = c("AA", "AB", "BB"),
  observed_count = c(AA_obscount, AB_obscount, BB_obscount),
  expected_freq = c(AA_expfreq, AB_expfreq, BB_expfreq),
  expected_count = expected_freq * N
)

# Perform Pearson's chi-squared test of goodness of fit
chi_sq <- sum(
  (genotype_data$observed_count - genotype_data$expected_count)^2 / genotype_data$expected_count
  )

p_value <- 1 - pchisq(chi_sq, df = 1)

chi_sq
p_value

# Note: The chi-square test is not actually recommended for such low sample sizes; 
# we use it here because it is appropriate for most real world data

# How would you interpret the results of the test?


### The degree of deviation from HWE can be expressed as the fixation index
# Fis = (He - Ho) / He
        # He: expected heterozygosity
        # Ho: observed heterozygosity
# Fis > 0: heterozygote deficiency, Fis < 0 heterozygote excess

(Fis <- (AB_expfreq - AB_obsfreq) / AB_expfreq)



### ============================================================================
### Exercise 2: Using the Genepop format to test for HWE


### Install and load required R packages
install.packages("adegenet")
install.packages("pegas")
library(adegenet)
library(pegas)


### Read in data from Genepop format file
help(read.genepop)
starfish <- read.genepop("data/msats/starfish.gen", ncode = 1)


### Compare with input text file on GitHub


### Access data in the new genind object
starfish
starfish@tab
starfish@loc.n.all


### Quick test for HWE (pegas package)
hw.test(starfish)



### ============================================================================
### Exercise 3: Real-world microsatellite data of Caribbean reef fish populations


### Review Genepop text file on GitHub / text editor


### Read in data from Genepop format file
barbados <- read.genepop("data/msats/puella_barbados.gen", ncode = 3)
barbados


### What is the most / least diverse locus in terms of number of alleles?
# Hint: Use the @loc.n.all slot of the genind object
#


### Locus summary using poppr (no. alleles, Simpson's index, heterozygosity, evenness)
install.packages("poppr")
install.packages("genepop")
library(poppr)
library(genepop)

locus_table(barbados)

# Note: Simpson's index and Evenness measures allelic diversity by considering both
# the number of alleles and their relative abundances, with values ranging from 0 to 1
# E.g. Simpson's index: probability that two randomly selected alleles are different


### Is this population in HWE? What does this tell us?


### More detailed test function for HWE, with results over all alleles (genepop package)
test_HW("data/msats/puella_barbados.gen", outputFile = "work/HW_barbados.txt")

# Note: The exact HWE test compares the observed genotype distribution to all possible distributions 
# with the same allele counts, and uses Markov chain sampling to estimate how often equally or more 
# extreme configurations occur under equilibrium.
# How many genotype configurations are as unlikely or more unlikely than the observed one?


### Compare these results to the Cayo de Media Luna population
test_HW("data/msats/puella_medialuna.gen", outputFile = "work/HW_medialuna.txt")


### Test all Caribbean populations. What patterns can you identify regarding loci and populations?
# Use data file: data/msats/puella_caribbean.gen
# test_HW()


### What patterns can you identify regarding individual loci and populations,
# and what could be the underlying reasons?



### ============================================================================
### Solutions — Exercise 1

# Calculate the frequency of the remaining genotypes
BB_obscount <- 4
BB_obsfreq <- BB_obscount / N

AB_obscount <- 2
AB_obsfreq <- AB_obscount / N


# What are the allele frequencies (p for allele A, q for allele B)?
q <- (BB_obscount * 2 + AB_obscount) / (N * 2)   # or q <- 1 - p


# Calculate expected genotype frequencies under HWE
AB_expfreq <- 2 * p * q
BB_expfreq <- q ^ 2


# How would you interpret the results of the test?
# - With a p-value of 0.02, the null hypothesis of Hardy–Weinberg equilibrium can be rejected
# - The fixation index of Fis = 0.66 indicates a strong deficiency of heterozygotes
# - Reasons could be e.g. inbreeding, selection against heterozygotes, or population structure



### ============================================================================
### Solutions — Exercise 3

### What is the most / least diverse locus in terms of number of alleles?
barbados@loc.n.all   # most diverse: pam013 (29 alleles), least diverse: hyp015 / e2 (6 alleles)


### Test all Caribbean populations. What patterns can you identify regarding loci and populations?
test_HW("data/msats/puella_caribbean.gen", outputFile = "work/HW_caribbean.txt")


### What patterns can you identify regarding individual loci and populations?
# - Some loci show significant deviations from HWE in multiple populations, 
# which could indicate selection or technical issues with the locus (e.g. genotyping errors)
# - Some populations show significant deviations from HWE across multiple loci,
# which could indicate population structure, inbreeding, or demographic history



### ============================================================================
### Links and resources

# R cheatsheet: https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# R for Beginners: YaRrr! The Pirate's Guide to R, https://bookdown.org/ndphillips/YaRrr/
# Advanced R: R for Data Science, https://r4ds.had.co.nz/index.html
