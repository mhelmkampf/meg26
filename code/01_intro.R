### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 01. Introduction                                                         ###
### ======================================================================== ###


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Establish SSH connection to cluster using your course account (please adjust)
ssh user1234@rosa.hpc.uni-oldenburg.de


### Download course materials to course account using Git
git clone https://github.com/mhelmkampf/meg26.git


### Load RStudio module
module load RStudio-Server


### Execute script to start RStudio
rstudio-start-on-rosa.sh


### Copy the SSH command provided at 1) and execute in new terminal tab / window
#> ssh -N -L 8000: ...
#> re-enter your password, and note nothing will change in the terminal on success


### Go to http://localhost:8000 in a browser window



### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### ============================================================================
### Exercise 1: Manual test for HWE


### What are the observed genotype frequencies?

# Assign value to variable
N <- 12

# Recall variable
N

# Arithmetic calculation
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

# Enter data into a table ("tibble") and add expected genotype counts
library(tidyverse)

dat <- tibble(
  genotype = c("AA", "AB", "BB"),
  observed_count = c(AA_obscount, AB_obscount, BB_obscount),
  expected_freq = c(AA_expfreq, AB_expfreq, BB_expfreq),
  expected_count = expected_freq * N
)

# Perform Pearson's chi-squared test of goodness of fit
chi_sq <- sum(
  (dat$observed_count - dat$expected_count)^2 / dat$expected_count
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

(Fis <- (AB_expfreq - AB_obsfreq) / AB_expfreq)



### ============================================================================
### Solutions

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
# - Reasons could be e.g. genetic drift, inbreeding, or selection against heterozygotes



### ============================================================================
### Links and resources

# R cheatsheet: https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# R for Beginners: YaRrr! The Pirate's Guide to R, https://bookdown.org/ndphillips/YaRrr/
# Advanced R: R for Data Science, https://r4ds.had.co.nz/index.html
