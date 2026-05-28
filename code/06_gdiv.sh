### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2026                             ###
### 06. Population genomics and genetic diversity                            ###
### ======================================================================== ###


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Establish SSH connection to cluster using your course account (please adjust user name)
ssh -Y user1234@rosa.hpc.uni-oldenburg.de   # -Y enables graphical windows


### Update git repository
cd meg26
git pull


### ALTERNATIVELY: If you receive an error, delete and re-download the repository
# cd
# rm -rf meg26
# git clone https://github.com/mhelmkampf/meg26.git


### Prepare today's working directory
cd work/snps
# move all log files into "logs" subdirectory
#>



### ============================================================================
### Exercise 1:

module load BCFtools
module load VCFtools


### Create link to today's SNP dataset (only if file does not exist already)
# ln -s /nfs/data/haex1482/shared/course/hamlets_filt_lg12.vcf.gz hamlets_filt_lg12.vcf.gz


### Copy population id files to working directory
cp ../../meg25/data/meta/*.txt .


### Calculate heterozygosity and Fis per individual
vcftools \
  --gzvcf hamlets_filt_lg12.vcf.gz \
  --keep hamlet_puebel.txt \
  --het \
  --stdout > Het_puebel.tsv

vcftools \
  --gzvcf hamlets_filt_lg12.vcf.gz \
  --keep hamlet_atlhua.txt \
  --het \
  --stdout > Het_atlhua.tsv

column -t Het_puebel.tsv

### From here, we could average and plot Fis per population


ln -s /nfs/data/haex1482/shared/course/guitar_filt_10ctg.vcf.gz guitar_filt_10ctg.vcf.gz

### Calculate heterozygosity and Fis per individual

vcftools \
  --gzvcf guitar_filt_10ctg.vcf.gz \
  --keep guitar_spain.txt \
  --het \
  --stdout > Het_spain.tsv

vcftools \
  --gzvcf guitar_filt_10ctg.vcf.gz \
  --keep guitar_cabov.txt \
  --het \
  --stdout > Het_cabov.tsv


### ============================================================================
### Exercise 2: Assess genetic diversity (nucleotide diversity)

### To calculate nucleotide diversity, we need all sites, not just SNPs:
ln -s /nfs/data/haex1482/shared/course/hamlets_all_1M.vcf.gz hamlets_all_1M.vcf.gz


### Copy population id files to working directory
cp ../../meg25/data/gdiv/*.txt .


### Measure nucleotide divergence per species / population
for i in pop_*.txt; do
vcftools \
  --gzvcf hamlets_LG12-1M_all.vcf.gz \
  --keep "$i" \
  --site-pi \
  --stdout > Pi_"$i".tsv
done


### Average pi per population
for i in Pi_*.tsv; do
  echo $i
  awk '$3 != "-nan" { sum += $3 ; next } END { print sum / NR }' $i
  echo
done



### ============================================================================
### Exercise 3: Filter by linkage disequilibrium

### Referring to the linkage decay plot, how many basepairs (bp) apart should the 
### sites be from each other to avoid linkage disequilibrium? Use vcftools --thin
### to create such a dataset


### Calculate the correlation coefficient of linkage disequilibrium in the first
### 100000 bp of both datasets (before and after thinning)
vcftools --gzvcf hamlets_LG12_snp.vcf.gz \
  --hap-r2 \
  --stdout > R2_before_thinning.tsv

vcftools --gzvcf hamlets_LG12_snp_2kb.vcf.gz \
  --hap-r2 \
  --stdout > R2_after_thinning.tsv


### Plot R2 distribution before and after thinning, in ASCII
for file in R2_before_thinning.tsv R2_after_thinning.tsv; do
  echo
  echo "==> $file <=="
  awk '
  BEGIN { bin_size=0.05 }
  /^CHR/ { next }  # skip header
  {
    if ($5 != "-nan") {
      r2 = $5 + 0
      bin = int(r2 / bin_size)
      counts[bin]++
      total++
    }
  }
  END {
    max_percent = 0
    for (i in counts) {
      percent = 100 * counts[i] / total
      if (percent > max_percent) max_percent = percent
    }
    for (i = 0; i <= int(1 / bin_size); i++) {
      percent = 100 * counts[i] / total
      bar_len = int((percent / max_percent) * 50)
      bin_label = sprintf("%.2f", i * bin_size)
      printf "%5s | %s (%.1f%%)\n", bin_label, substr("##################################################", 1, bar_len), percent
    }
  }' "$file"
  echo
done



### ============================================================================
### Solutions:

### VCF stats
bcftools stats hamlets_LG12_snp.vcf.gz | grep 'SN'
#> 1032804 sites

zcat hamlets_LG12_snp.vcf.gz | grep -c -v '#'


### Include only sites with minor allele count of at least 2
vcftools \
    --gzvcf hamlets_LG12_snp.vcf.gz \
    --mac 2 \
    --recode \
    --stdout | bgzip > hamlets_LG12_snp_mac2.vcf.gz


### Summarize
bcftools stats hamlets_LG12_snp_mac2.vcf.gz | grep 'SN'
#> 319960 sites


### Spot check
zcat hamlets_LG12_snp_mac2.vcf.gz | tail


### Thin the filtered data so that no two sites are within 2000 bp from each other
vcftools \
    --gzvcf hamlets_LG12_snp.vcf.gz \
    --thin 2000 \
    --recode \
    --stdout | bgzip > hamlets_LG12_snp_2kb.vcf.gz



### ============================================================================
### Alternative code and notes

### Summarize in ASCII
for file in R2_before_thinning.tsv R2_after_thinning.tsv; do
  echo
  echo "$file:"
  awk '$5 != "-nan" { sum += $5; n++; vals[n] = $5 }
  END {
    if (n == 0) exit
    asort(vals)
    median = (n % 2) ? vals[int((n+1)/2)] : (vals[n/2] + vals[n/2 + 1]) / 2
    print "  Mean R²: " sum/n
    print "  Median R²: " median
    print "  Count: " n
  }' "$file"
done
