### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 09. DNA barcoding                                                        ###
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


### Create and navigate to today's working directory
cd work
mkdir barcode
cd barcode



### ============================================================================
### Exercise 1: Extract barcode sequene from Sanger reads

### Background – PCR / sequencing primer cocktail used for COI:

# VF2_t1    tgtaaaacgacggccagtCAACCAACCACAAAGACATTGGCAC
# FishF2_t1 tgtaaaacgacggccagtCGACTAATCATAAAGATATCGGCAC

# FishR2_t1 caggaaacagctatgacACTTCAGGGTGACCGAAGAATCAGAA
# FR1d_t1   caggaaacagctatgacACCTCAGGGTGTCCGAARAAYCARAA

# Amplifies approx. 650 bp fragment in 5' region of COI
# Primers are M13-tailed to allow DNA sequencing
# (Ivanova et al. 2007, Molecluar Ecology Notes, lower case = tail)


### You will be working with one of the following 11 samples:
### 10, 11, 14, 15, 19, 20, 22, 26, 27, 28, 35


### From meg25/data/barcode, copy your pair of ABI trace files to the working directory
cp ~/meg25/data/barcode/seq_10F.ab1 .   # adjust file name here and below
cp ~/meg25/data/barcode/seq_10R.ab1 .   # adjust file name here and below


### View sequence (download raw file from https://github.com/mhelmkampf and upload here:)
#> https://www.gear-genomics.com/teal/


### Extract sequence from trace files
module load EMBOSS

seqret -sequence seq_10F.ab1 -outseq seq_10F.fas
seqret -sequence seq_10R.ab1 -outseq seq_10R.fas


### Trim sequences (adjust limits according to your case)
module load SeqKit

seqkit subseq -r 65:580 seq_10F.fas > trim_10F.fas
seqkit subseq -r 50:470 seq_10R.fas > trim_10R.fas


### Reverse complement reverse sequence
seqkit seq -r -p trim_10R.fas > trim_10R-rc.fas


### Merge foward and reverse sequences to consensus

# Usage:
# merger \
#   -asequence <forward.fas> \
#   -bsequence <reverse.fas> \
#   -outfile <consensus.txt> \
#   -outseq <consensus.fas>

#>



### ============================================================================
### Exercise 2: Match sequence to Barcode of Life Data Systems (BOLD)


### Visit https://v4.boldsystems.org/index.php
### Note: https://boldsystems.org (version 5) does not have full functionality yet


### Select "Identification v4" | Animal Identification (COI) | Species Level Barcode Records (default)


### Paste your consensus sequence and submit


### Explore the results:
### - What species has been identified?
### - Are there multiple species in the top 20? What is the similarity distribution?
### – Is it a plausible match considering the samples were obtained from seafood?


### Optional: Learn more about your species at https://fishbase.de


### Evaluate the reliability of your result by selecting its "Bin page"
### - What is the maximum p-distance within the species?
### – How does the distance distribution compare to the nearest neighbor? Is there a "barcode gap"?
### - How confident are you overall in your id?


### Optional: Identify barcode with BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi



### ============================================================================
### Solutions:

### Create consensus
merger \
  -asequence trim_10F.fas \
  -bsequence trim_10R-rc.fas \
  -outfile consensus_10.txt \
  -outseq consensus_10.fas


### All consensus sequences can be found in data/barcode/co1.fas


### Species IDs, BIN, barcode gap

# 10: Oncorhynchus keta (Chum salmon / Keta-Lachs), BOLD:AAA3872, large barcode gap
# 11: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), BOLD:AAA3069, small barcode gap
# 14: Platichthys flesus (European flounder / Flunder) or Pleuronectes sp., cannot be identified to species
# 15: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), BOLD:AAA3069, small barcode gap
# 19: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), BOLD:AAA3069, small barcode gap
# 20: Platichthys flesus (European flounder / Flunder) or Pleuronectes sp., cannot be identified to species
# 22: Gadus morhua (Atlantic cod / Kabeljau) or Gadus sp., small barcode gap / cannot be identified to species
# 26: Gadus morhua (Atlantic cod / Kabeljau) or Gadus sp. small barcode gap /cannot be identified to species
# 27: Pollachius virens (Pollock / Köhler), BOLD:AAB2980, large barcode gap
# 28: Oncorhynchus keta (Chum salmon / Keta-Lachs), BOLD:AAA3872, large barcode gap
# 35: Coregonus sp. (Cisco), cannot be identified to species