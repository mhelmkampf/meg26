### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 04. Genome sequencing and assembly                                       ###
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



### ============================================================================
### Exercise 1: Inspect raw sequecing reads on the command line

### Reorganize your work directory
cd work        # change into work/ directory
ls -l          # list directory contents
mkdir msats    # make new subdirectory called "msats"
mv *.* msats   # mv all files into msats


### Create new subdirectory called "asm" in work/ and navigate there
#> 


### Copy raw sequencing data to working directory, and decompress files
cp ~/meg26/data/asm/*.fastq.gz .
gzip -d *.fastq.gz


### Print whole file
cat HypPue1_illumina_raw_F.fastq


### Print only first four lines
head -n 4 HypPue1_illumina_raw_F.fastq
head -n 4 HypPue2_pacbio_hifi.fastq


### How do the Illumina and PacBio reads differ?


### How could we count the number of sequences?
# Hint: grep -c '<pattern>' HypPue1_illumina_raw_F.fastq (replace <pattern>)
#> 



### ============================================================================
### Exercise 2: Compare genome assemblies and calculate assembly metrics

### Create link to hamlet genome assemblies
ln -s /nfs/data/haex1482/shared/course/HypPue1_assembly_ref.fas HypPue1_assembly_ref.fas
ln -s /nfs/data/haex1482/shared/course/HypPue2_assembly_pacbio.fas HypPue2_assembly_pacbio.fas


### Print assembly / first 10 lines to screen
cat HypPue1_assembly_ref.fas
head HypPue1_assembly_ref.fas


### Print only first line (header of first sequence)
#> 


### Install assembly_stats package for Python
module load Python
pip install assembly_stats


### Calculate basic assembly metrics
assembly_stats HypPue1_assembly_ref.fas
assembly_stats HypPue2_assembly_pacbio.fas


### How do the two assemblies differ?


### Create link to a third genome assembly (Blackchin guitarfish)
ln -s /nfs/data/haex1482/shared/course/GlaCem1_assembly_pacbio.fas GlaCem1_assembly_pacbio.fas


### How does this assembly compare to the two hamlet assemblies?



### ============================================================================
### Exercise 3: Assemble genome using hifiasm

mkdir hifiasm
cd hifiasm

module load hifiasm

hifiasm -h


### Create link to subset of PacBio HiFi reads (30000 reads, about 1 GB)
ln -s /nfs/data/haex1482/shared/course/HypPue2_pacbio_30k.fastq.gz HypPue2_pacbio_30k.fastq.gz


### Assemble with hifiasm, a fast and accurate assembler for Hifi reads
hifiasm \
  -o test30k \
  --primary \
  -t 8 \
  HypPue2_pacbio_30k.fastq.gz


### Convert assembly graph to fasta format and calculate metrics
awk '/^S/{ print ">"$2; print $3 }' test30k.p_ctg.gfa > test30k.p_ctg.fas

#>



### ============================================================================
### Optional: Search for microsatellites

cd ~/work/asm

### Highlight AC and GT repeats with at least 10 units in the first 100000 bp
head -n 2 HypPue1_assembly_ref.fas | cut -c -100000 | grep -E '(AC|GT){10,}'



### ============================================================================
### Solutions: Exercise 1

### Create new directory and navigate there
mkdir asm
cd asm


### How do the Illumina and PacBio reads differ?
# - Illumina reads are short, and of equal length
# - PacBio reads are much longer, and of variable length


### Count number of reads in Fastq file
grep -c '^+$' HypPue1_illumina_raw_F.fastq
grep -c '^+$' HypPue2_pacbio_hifi.fastq


### ----------------------------------------------------------------------------
### Solutions: Exercise 2


### Print only first line (header of first sequence)
head -n 1 HypPue1_assembly_ref.fas
head -n 1 HypPue2_assembly_pacbio.fas


### How do the two assemblies differ?
# - the reference genome (Illumina / PacBio hybrid) contains a much larger number of contigs / scaffolds
# - contiguity (N50) is very high for both assemblies
# - the reference genome achieves this by scaffolding, the PacBio assembly by using long reads


### How does the guitarfish genome assembly compare to the two hamlet assemblies?
assembly_stats GlaCem1_assembly_pacbio.fas


### ----------------------------------------------------------------------------
### Solutions: Exercise 3

### Calculate assembly metrics
assembly_stats test30k.p_ctg.fas
