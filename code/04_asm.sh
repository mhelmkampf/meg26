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
### Exercise 1: Handle genome assemblies on the command line

### Reorganize your work directory
cd work        # change into work/ directory
ls -l          # list directory contents
mkdir msats    # make new subdirectory called "msats"
mv *.* msats   # mv all files into msats


### Make new subdirectory called "asm" in work/ and navigate there
#>


### Create link to hamlet genome assemblies
ln -s /fs/dss/home/haex1482/share/HypPue1_assembly_ref.fas HypPue1_assembly_ref.fas
ln -s /fs/dss/home/haex1482/share/HypPue2_assembly_pacbio.fas HypPue2_assembly_pacbio.fas


### Print assembly / first 10 lines to screen
cat HypPue1_assembly_ref.fas
head HypPue1_assembly_ref.fas


### Print only first line (header of first sequence)
### Review usage of "head" with: man head
#>


### Look at a slice of underlying read data
cp ~/meg26/data/asm/*.fastq.gz .
gzip -d *.fastq.gz

cat HypPue1_illumina_raw_F.fastq

head -n 4 HypPue1_illumina_raw_F.fastq
head -n 4 HypPue2_pacbio_hifi.fastq


### How do the underlying Illumina and PacBio reads differ?


### How could we count the number of sequences?
### Use: grep -c '<pattern>' HypPue1_illumina_raw_F.fastq
#>



### ============================================================================
### Exercise 2: Assess read quality and trim Illumina reads (optional)

module load FastQC
module load cutadapt

mkdir fastqc


### Run FastQC, a quality control tool for high-throughput sequence reads
### https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
fastqc -o fastqc HypPue1_illumina_raw_F.fastq HypPue1_illumina_raw_R.fastq
fastqc -o fastqc HypPue2_pacbio_hifi.fastq


### Open a new terminal window on your local computer and download the HTML results
scp user1234@rosa.hpc.uni-oldenburg.de:/user/user1234/work/fastq/*.html .


### Open HTML results with browser


### Alternatively, print results as text only
# lynx -dump fastqc/HypPue1_illumina_raw_F_fastqc.html


### Back on the cluster, remove adapters and low quality bases from Illumina reads
cutadapt -h

cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -q 20 \
  -m 120 \
  -o HypPue1_illumina_trimmed_F.fastq \
  HypPue1_illumina_raw_F.fastq


### Re-run FastQC on trimmed reads and compare before and after
#>



### ============================================================================
### Exercise 3: Calculate and compare assembly metrics

### Install assembly_stats package for Python
pip install assembly_stats

assembly_stats HypPue1_assembly_ref.fas
assembly_stats HypPue2_assembly_pacbio.fas


### How do the two assemblies differ?



### ============================================================================
### Exercise 4: Assemble genome using hifiasm

mkdir hifiasm

module load hifiasm

hifiasm -h


### Create link to subset of PacBio HiFi reads (30000 reads, about 1 GB)
ln -s /fs/dss/home/haex1482/share/HypPue2_pacbio_30k.fastq.gz HypPue2_pacbio_30k.fastq.gz


### Assemble with hifiasm, a fast and accurate assembler for Hifi reads
hifiasm \
  -o hifiasm/test30k \
  --primary \
  -t 10 \
  HypPue2_pacbio_30k.fastq.gz


### Convert assembly graph to fasta format and calculate metrics
cd hifiasm

awk '/^S/{print ">"$2; print $3}' test30k.p_ctg.gfa > test30k.p_ctg.fas

#>



### ============================================================================
### Optional: Search for microsatellites

cd ~/work/asm

### Highlight AC and GT repeats with at least 10 units in the first 100000 bp
head -n 2 HypPue1_assembly_ref.fas | cut -c -100000 | grep -E '(AC|GT){10,}'



### ============================================================================
### Solutions: Exercise 1

### New directory
mkdir asm
cd asm


### Print first line only
head -n 1 HypPue1_assembly_ref.fas


### Counting reads in Fastq file
grep -c '^+$' HypPue1_illumina_raw_F.fastq
grep -c '^+$' HypPue2_pacbio_hifi.fastq



### ----------------------------------------------------------------------------
### Solutions: Exercise 2

### Re-run FastQC on trimmed reads
fastqc -o fastqc HypPue1_illumina_trimmed_F.fastq



### ----------------------------------------------------------------------------
### Solutions: Exercise 4

### Calculate assembly metrics
assembly_stats test30k.p_ctg.fas
