Downloaded raw data from SRA. 

Files are: 
SRS3346514
SRS3346515
SRS3346516
SRS3346517
SRS3346518
SRS3346519
SRS3346520
SRS3346521
SRS3346522

Create BWA index of Xmac genome: 
GCF_002775205.1_X_maculatus-5.0-male_genomic.fna

Headers of downloaded SRA data are malformatted for bwa, fix them.
```
zcat SRS3346514_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346514_1.fix.fastq 
zcat SRS3346514_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346514_2.fix.fastq 
zcat SRS3346515_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346515_1.fix.fastq 
zcat SRS3346515_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346515_2.fix.fastq 
zcat SRS3346516_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346516_1.fix.fastq 
zcat SRS3346516_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346516_2.fix.fastq 
zcat SRS3346517_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346517_1.fix.fastq 
zcat SRS3346517_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346517_2.fix.fastq 
zcat SRS3346518_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346518_1.fix.fastq 
zcat SRS3346518_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346518_2.fix.fastq 
zcat SRS3346519_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346519_1.fix.fastq 
zcat SRS3346519_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346519_2.fix.fastq 
zcat SRS3346520_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346520_1.fix.fastq 
zcat SRS3346520_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346520_2.fix.fastq 
zcat SRS3346521_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346521_1.fix.fastq 
zcat SRS3346521_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346521_2.fix.fastq 
zcat SRS3346522_1.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346522_1.fix.fastq 
zcat SRS3346522_2.fastq.gz | perl -ne 's/\.([12]) /\/$1 /; print $_'  > SRS3346522_2.fix.fastq 
```
```
gzip *.fix.fastq
```

BWA MEM
```
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346514_1.fix.fastq.gz SRS3346514_2.fix.fastq.gz -o SRS3346514.sam
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346515_1.fix.fastq.gz SRS3346515_2.fix.fastq.gz -o SRS3346515.sam
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346516_1.fix.fastq.gz SRS3346516_2.fix.fastq.gz -o SRS3346516.sam
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346517_1.fix.fastq.gz SRS3346517_2.fix.fastq.gz -o SRS3346517.sam
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346518_1.fix.fastq.gz SRS3346518_2.fix.fastq.gz -o SRS3346518.sam
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346519_1.fix.fastq.gz SRS3346519_2.fix.fastq.gz -o SRS3346519.sam
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346520_1.fix.fastq.gz SRS3346520_2.fix.fastq.gz -o SRS3346520.sam
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346521_1.fix.fastq.gz SRS3346521_2.fix.fastq.gz -o SRS3346521.sam
bwa mem GCF_002775205.1_X_maculatus-5.0-male_genomic.fna SRS3346522_1.fix.fastq.gz SRS3346522_2.fix.fastq.gz -o SRS3346522.sam
```

Convert sam to bam and sort.
```
samtools view -bSh --threads 20 raw_data/SRS3346513.sam | samtools sort --threads 20 -o raw_data/SRS3346513.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346514.sam | samtools sort --threads 20 -o raw_data/SRS3346514.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346515.sam | samtools sort --threads 20 -o raw_data/SRS3346515.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346516.sam | samtools sort --threads 20 -o raw_data/SRS3346516.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346517.sam | samtools sort --threads 20 -o raw_data/SRS3346517.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346518.sam | samtools sort --threads 20 -o raw_data/SRS3346518.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346519.sam | samtools sort --threads 20 -o raw_data/SRS3346519.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346520.sam | samtools sort --threads 20 -o raw_data/SRS3346520.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346521.sam | samtools sort --threads 20 -o raw_data/SRS3346521.sorted.bam - 
samtools view -bSh --threads 20 raw_data/SRS3346522.sam | samtools sort --threads 20 -o raw_data/SRS3346522.sorted.bam - 
```

Check mapping stats, if they look reasonable remove the sam files and the original (not fixed) fastq files
```
module load picard/2.21.2

mkdir stats

picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346513.sorted.bam O=stats/SRS3346513.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346514.sorted.bam O=stats/SRS3346514.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346515.sorted.bam O=stats/SRS3346515.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346516.sorted.bam O=stats/SRS3346516.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346517.sorted.bam O=stats/SRS3346517.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346518.sorted.bam O=stats/SRS3346518.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346519.sorted.bam O=stats/SRS3346519.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346520.sorted.bam O=stats/SRS3346520.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346521.sorted.bam O=stats/SRS3346521.summary
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=raw_data/SRS3346522.sorted.bam O=stats/SRS3346522.summary
```

Create file with bam locations
vi bam_locations.txt

```
raw_data/SRS3346513.sorted.bam SRS3346513
raw_data/SRS3346514.sorted.bam SRS3346514
raw_data/SRS3346515.sorted.bam SRS3346515
...
```

Loop over files to index bam files, collect bamtools stats and add read groups
```
for i in {1..10}

do

bamfile=$(cat bam_locations.txt | sed -n ''$i'p' | awk '{print $1}')
samplename=$(cat bam_locations.txt | sed -n ''$i'p' | awk '{print $2}')

# Index the bam files
samtools index ${bamfile}

# Check alignment statistics using bamtools
bamtools stats -in ${bamfile} > ${samplename}.STATS

# AddOrReplaceReadGroups: Add readgroups to the bam files for SNP calling
module load picard
picard AddOrReplaceReadGroups INPUT=${bamfile} OUTPUT=${samplename}.readgroup.bam SORT_ORDER=coordinate RGPL=Illumina RGLB=${samplename} RGSM=${samplename} RGPU="Illumina"

samtools index ${samplename}.readgroup.bam

done
```


Call SNPs with Unified Genotyper for all files (only one shown here)
```
module load gatk/3.8.0
gatk -T UnifiedGenotyper \
-R GCF_002775205.1_X_maculatus-5.0-male_genomic.fna \
-o SRS3346520.vcf.gz \
--output_mode EMIT_ALL_SITES \
-I SRS3346520.readgroup.bam
```

Merge and index merged vcf 
```
module load htslib/1.8
vcf-merge SRS3346513.vcf.gz SRS3346514.vcf.gz SRS3346515.vcf.gz SRS3346516.vcf.gz SRS3346517.vcf.gz SRS3346518.vcf.gz SRS3346519.vcf.gz SRS3346520.vcf.gz SRS3346521.vcf.gz SRS3346522.vcf.gz | bgzip -c > merged.vcf.gz
tabix -p vcf merged.vcf.gz
```

Select biallelic sites
```
module load gatk/3.8.0
gatk -T SelectVariants \
    -R GCF_002775205.1_X_maculatus-5.0-male_genomic.fna \
    -V merged.vcf.gz \
    -selectType SNP \
    -restrictAllelesTo BIALLELIC \
    -o merged_just_biallelic_snps.vcf
```

Apply GATK filters
```
module load gatk/3.8.0
gatk -T VariantFiltration \
    -R GCF_002775205.1_X_maculatus-5.0-male_genomic.fna \
    -V merged_just_biallelic_snps.vcf \
    -filter "QD<2.0||FS>60.0||MQ<40.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" \
    --filterName "AllFilter" \
    -o new_filter_applied_biallelic_snps.vcf
```

Use SelectVariants to exclude sites that didn't pass the hard filters
```
module load gatk/3.8.0
gatk -T SelectVariants \
    -R GCF_002775205.1_X_maculatus-5.0-male_genomic.fna \
    -V new_filter_applied_biallelic_snps.vcf \
    --excludeFiltered \
    -o filter_excluded_biallelic_snps.vcf
```


Additional filtering with vcftools
```
vcftools --gzvcf filter_excluded_biallelic_snps.vcf.gz --min-alleles 2 --max-alleles 2 --minQ 30 --minDP 10 --max-missing 0.9 --recode --recode-INFO-all --out merged_snps_minq30_minDP10_maxmiss0.9
```

Gzip file
```
gzip merged_snps_minq30_minDP10_maxmiss0.9.recode.vcf
```

Remove fixed nonreference sites, first attempt did not seem to work, tried again with a different flag
```
vcftools --gzvcf merged_snps_minq30_minDP10_maxmiss0.9.recode.vcf.gz --max-non-ref-ac 19 --recode --recode-INFO-all --out merged_snps_minq30_minDP10_maxmiss0.9.nofixednonref
gzip merged_snps_minq30_minDP10_maxmiss0.9.nofixednonref.recode.vcf
# After looking at Fst, there are sites with -nan for Fst (all fixed non-ref) so the filtering did not work the way I anticipated. I only want variable sites. 
rm merged_snps_minq30_minDP10_maxmiss0.9.nofixednonref.recode.vcf.gz

vcftools --gzvcf merged_snps_minq30_minDP10_maxmiss0.9.recode.vcf.gz --max-non-ref-af 0.999 --recode --recode-INFO-all --out merged_snps_minq30_minDP10_maxmiss0.9.maxnonrefaf0.999
gzip merged_snps_minq30_minDP10_maxmiss0.9.maxnonrefaf0.999.recode.vcf

```

Looks like many sites are 0/1 in all individuals. Do we need more stringent filtering in the GATK step?
Calculate Fst - generate two files to define the populations
```
sulfidic.txt
SRS3346513
SRS3346516
SRS3346517
SRS3346518
SRS3346519

nonsulfidic.txt
SRS3346514
SRS3346515
SRS3346520
SRS3346521
SRS3346522

vcftools --gzvcf merged_snps_minq30_minDP10_maxmiss0.9.maxnonrefaf0.999.recode.vcf.gz --weir-fst-pop sulfidic.txt --weir-fst-pop nonsulfidic.txt --out sulfvsnonsulf.fst
gzip sulfvsnonsulf.fst.weir.out

vcftools --gzvcf merged_snps_minq30_minDP10_maxmiss0.9.maxnonrefaf0.999.recode.vcf.gz --weir-fst-pop sulfidic.txt --weir-fst-pop nonsulfidic.txt --out sulfvsnonsulf.5kb.1kbSTEP.fst --fst-window-size 5000 --fst-window-step 1000
```

### Fst sliding window no overlap
```
vcftools --gzvcf merged_snps_minq30_minDP10_maxmiss0.9.maxnonrefaf0.999.recode.vcf.gz --weir-fst-pop sulfidic.txt --weir-fst-pop nonsulfidic.txt --out sulfvsnonsulf.25kb.nooverlap.fst --fst-window-size 25000 --fst-window-step 25000
cp sulfvsnonsulf.25kb.5kbstep.fst.windowed.weir.fst sulfvsnonsulf.25kb.5kbstep.fst.windowed.weir.bed
# remove first line (header) from sulfvsnonsulf.25kb.5kbstep.fst.windowed.weir.bed to convert to bed file
gzip sulfvsnonsulf.25kb.5kbstep.fst.windowed.weir.fst
```

### Fst sliding window 5kb step
```
vcftools --gzvcf merged_snps_minq30_minDP10_maxmiss0.9.maxnonrefaf0.999.recode.vcf.gz --weir-fst-pop sulfidic.txt --weir-fst-pop nonsulfidic.txt --out sulfvsnonsulf.25kb.5kbstep.fst --fst-window-size 25000 --fst-window-step 5000
```

Create bed file of top Fst outliers (either 1% or 0.1%)
```
bedtools intersect -a top1FST.bed -b GCF_002775205.1_X_maculatus-5.0-male_genomic.gff
bedtools intersect -a GCF_002775205.1_X_maculatus-5.0-male_genomic.gff -b top1FST.bed


grep -P "\tgene\t" GCF_002775205.1_X_maculatus-5.0-male_genomic.gff > GCF_002775205.1_X_maculatus-5.0-male_genomic.onlyGenes.gff

bedtools intersect -a GCF_002775205.1_X_maculatus-5.0-male_genomic.onlyGenes.gff -b top1FST.bed

bedtools intersect -a GCF_002775205.1_X_maculatus-5.0-male_genomic.onlyGenes.gff -b top1FST.bed > top1FST.genes.gff

bedtools intersect -a GCF_002775205.1_X_maculatus-5.0-male_genomic.onlyGenes.gff -b top0.01FST.bed > top0.01FST.genes.gff
```

Create the overlap between Fst and genes for all genes, first need to modify the Fst outfile from vcftools to remove the first line to turn it into a bed file

sulfvsnonsulf.25kb.5kbstep.fst.windowed.weir.bed

```
bedtools intersect -a GCF_002775205.1_X_maculatus-5.0-male_genomic.onlyGenes.gff -b sulfvsnonsulf.25kb.5kbstep.fst.windowed.weir.bed -wb -wa > fst_for_all_genes.txt
gzip fst_for_all_genes.txt
```



# Add Guppy SNPs to the analysis (ultimately unnecessary as we did not include guppy in the saguaro analysis)
```
fastq-dump --outdir sra_guppy_data --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR1171024

bwa mem -R '@RG\tID:guppy\tSM:guppy\tPU:lane\tPL:illumina' GCF_002775205.1_X_maculatus-5.0-male_genomic.fna sra_guppy_data/SRR1171024_pass_1.fastq.gz  sra_guppy_data/SRR1171024_pass_2.fastq.gz -o guppy.sam
samtools view -bSh --threads 20 guppy.sam | samtools sort --threads 20 -o guppy.sorted.bam - 


module load picard/2.21.2
picard CollectAlignmentSummaryMetrics R=GCF_002775205.1_X_maculatus-5.0-male_genomic.fna I=guppy.sorted.bam O=guppy.summary


samtools index guppy.sorted.bam

# Check alignment statistics using bamtools
bamtools stats -in guppy.sorted.bam > guppy.STATS

module load gatk/3.8.0
gatk -T UnifiedGenotyper \
-R GCF_002775205.1_X_maculatus-5.0-male_genomic.fna \
-o /scratch/joanna.l.kelley_626782/guppy.vcf.gz \
--output_mode EMIT_ALL_SITES \
-I guppy.sorted.bam

# Merge and index merged vcf

module load htslib/1.8
vcf-merge merged.vcf.gz guppy.vcf.gz | bgzip -c > combined_merged_pmexANDguppy.vcf.gz
tabix -p vcf combined_merged_pmexANDguppy.vcf.gz
```


# PhyloSNP
Need to create a file for the chromosome names because PhyloSNP expects numeric chromosome names
```
zcat combined_merged_pmexANDguppy.vcf.gz | sed 's/\|/ /'|awk '{print $1}' | uniq > header.txt
vi header.txt #(remove lines with #) 
nl -ba header.txt > chromosome_rename.txt
awk ' { t = $1; $1 = $2; $2 = t; print; } ' chromosome_rename.txt > chrom_rename_correct_order.txt

module load bcftools/1.6
bcftools annotate --rename-chrs chrom_rename_correct_order.txt combined_merged_pmexANDguppy.vcf.gz > \
merged.v2_pmexANDguppy_renamed_for_snphylo.vcf
/SNPhylo/snphylo.sh -v merged.v2_pmexANDguppy_renamed_for_snphylo.vcf -d snphylo.output.gds -o guppy -a 102 -b -B 1000
```

# Saguaro 

Ultimately ran this for the saguaro input
```
/saguaro/saguarogw-code/VCF2HMMFeature -i merged_just_biallelic_snps.vcf -o merged.saguaro 
```

Run Saguaro, -f indicates the input HMM file, -o indicates the output directory, -iter indicates the number of iterations to go through (29 results in 30 cacti)
```
/saguaro/saguarogw-code/Saguaro -f merged.saguaro -o saguaro_output -iter 29
```

convert the saguaro results into Phylip format (this will output in whatever your current directory is) however, the individual names are too long so I need to reformat first. 

```
cp saguaro.cactus saguaro.renamed.cactus
vi saguaro.renamed.cactus
# Replace SRS with S in vi using %s/SRS/S/g
```

```
/saguaro/saguarogw-code/Saguaro2Phylip -i saguaro_output_2019_12_17/saguaro.renamed.cactus -outgroup 0
```


Need to download Phylip (I use the GUI, http://evolution.genetics.washington.edu/phylip.html) to use the Neighbor module to go from Phylip files to neighbor joining trees
The final output should be trees in newick format

Need to identify which regions are best described to which cacti
```
grep score saguaro_output_2019_12_17/LocalTrees.out \
> cactus_locations.txt
```

Create the trees using /phylip-3.696/src/neighbor

 
 Separate the cactus_locations.txt files into columns to create a bed file with regions (cactus_locations.bed).
 Create bedfile of candidate regions from Michi (candidates.bed)

```
edtools intersect -a cactus_locations.bed -b candidates.bed -wb -wa > cactus_intersect_candidates.txt
 ```

## Remove guppy and run saguaro again
```
vcftools --gzvcf merged_just_biallelic_snps.vcf.gz --remove guppy --recode --out merged_just_biallelic_snps.noguppy.vcf
/data/kelley/projects/programs/saguaro/saguarogw-code/VCF2HMMFeature -i merged_just_biallelic_snps.noguppy.vcf.recode.vcf -o merged.noguppy.saguaro 

/data/kelley/projects/programs/saguaro/saguarogw-code/Saguaro -f merged.noguppy.saguaro -o saguaro_output_noguppy -iter 29

cp saguaro.cactus saguaro.renamed.cactus
vi saguaro.renamed.cactus
# Replace SRS with S in vi using %s/SRS/S/g

/data/kelley/projects/programs/saguaro/saguarogw-code/Saguaro2Phylip -i /scratch/joanna.l.kelley_267592/saguaro_output_noguppy/saguaro.renamed.cactus -outgroup 1

grep score /scratch/joanna.l.kelley_267592/saguaro_output_noguppy/LocalTrees.out \
> cactus_locations_noguppy.txt
```

Move cactus_locations_noguppy.txt to my laptop and manipulate in excel & R to get the sums: 
```
tt = read.csv("cactus_locations_noguppy.csv", header = T)
aggregate(tt$length, by=list(Category=tt$cactus), FUN=sum)
``` 

Create phylip trees for each cactus using phylip neighbor. An example script to run is included in this folder maketree.sh that was run for each cactus [0-29]

Separate the cactus_locations_noguppy.txt files into columns to create a bed file with regions (cactus_locations_noguppy.bed).
 Create bedfile of candidate regions from Michi (candidates.bed)
 ```
 bedtools intersect -a cactus_locations_noguppy.bed -b candidates.bed -wb -wa > cactus_intersect_candidates_noguppy.txt
```

Move cactus_intersect_candidates_noguppy.txt to my computer and manipulate in excel. Remove all redundant rows to figure out the cacti that span the region. Make a bed file for each region and import into IGV with the gene region GFF. 




