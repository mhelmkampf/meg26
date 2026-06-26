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


### You will be assigned one of the following 11 samples:
### 10, 11, 14, 15, 19, 20, 22, 26, 27, 28, 35


### Assign your sample to a variable for future use
sample=
echo ${sample}   # re-call variable


### Copy your pair of ABI trace files to the working directory
cp ~/meg26/data/barcode/read_${sample}F.ab1 .
cp ~/meg26/data/barcode/read_${sample}R.ab1 .


### View trace file
# Download raw file from https://github.com/mhelmkampf/meg26 and upload here:
#> https://www.gear-genomics.com/teal/


### Extract sequence from trace files
module load EMBOSS

seqret -sequence read_${sample}F.ab1 -outseq read_${sample}F.fas
seqret -sequence read_${sample}R.ab1 -outseq read_${sample}R.fas


### Verify output files are in Fasta format
#>


### Trim low-quality bases from start of read (adjust boundaries according to your case)
module load SeqKit

fstart=
fend=

seqkit subseq -r ${fstart}:${fend} read_${sample}F.fas > trim_${sample}F.fas


### Trim low-quality bases at end of read (adjust boundaries according to your case)
rstart=
rend=

seqkit subseq -r ${rstart}:${rend} read_${sample}R.fas > trim_${sample}R.fas


### Reverse complement reverse sequence
seqkit seq -r -p trim_${sample}R.fas > trim_${sample}R-rc.fas


### Merge foward and reverse sequences to consensus
merger \
  -asequence trim_${sample}F.fas \
  -bsequence trim_${sample}R-rc.fas \
  -outfile con_${sample}.txt \
  -outseq con_${sample}.fas



### ============================================================================
### Exercise 2: Match sequence to Barcode of Life Data Systems (BOLD)


### Visit https://boldsystems.org


### Select "Barcode ID"
# Use default database (Animal Library (Public)) and operating mode (Rapid Species Search)


### Paste your consensus sequence and submit


### Explore the results
### - What species has been identified?
### - Are there multiple species in the top 20? What is the similarity distribution?
### – Is it a plausible match considering the samples were obtained from seafood?


### Optional: Learn more about your species at https://fishbase.de


### Evaluate the reliability of your result on its "BIN page" (via PID > BIN ID)
### - What is the maximum p-distance within the species?
### – How does the distance distribution compare to the nearest neighbor? Is there a "barcode gap"?
### - How confident are you overall in your id?


### Optional: Identify barcode with BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi



### ============================================================================
### Solutions:

### Assign your sample to a variable for future use
sample=10


### Trim low-quality bases from start of read (adjust boundaries according to your case)
fstart=65
fend=580


### Trim low-quality bases from end of read (adjust boundaries according to your case)
rstart=50
rend=470


### All consensus sequences can be found in results/barcodes_co1.fas


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