# code for proc b paper

trim

```bash

cd ~/hudsonica_transcriptomics/data/trimmed

cp /data/programs/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa .

for i in $(ls /data/copepods/hudsonica_multigen_transcriptomics/usftp21.novogene.com/raw_data | grep -v 'Rawdata_Readme.pdf' )

do {

    read1=$(ls /data/copepods/hudsonica_multigen_transcriptomics/usftp21.novogene.com/raw_data/${i}/*.fq.gz | grep '_1.fq.gz')
    read2=$(ls /data/copepods/hudsonica_multigen_transcriptomics/usftp21.novogene.com/raw_data/${i}/*.fq.gz| grep '_2.fq.gz')
    base1=$(ls /data/copepods/hudsonica_multigen_transcriptomics/usftp21.novogene.com/raw_data/${i}/*.fq.gz| grep '_1.fq.gz' | cut -f 8- -d "/" | cut -f 1 -d ".")
    base2=$(ls /data/copepods/hudsonica_multigen_transcriptomics/usftp21.novogene.com/raw_data/${i}/*.fq.gz | grep '_2.fq.gz' | cut -f 8- -d "/" | cut -f 1 -d ".")


  java -jar /data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
   -threads 8 \
    ${read1} ${read2} \
    ${base1}.qc.fq.gz s1_se \
    ${base2}.qc.fq.gz s2_se \
    HEADCROP:12 \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:31

}

done

```
## assemble transcriptome

```bash
cd /data/project_data/RNAseq/assembly

screen

/data/popgen/trinityrnaseq-v2.13.2/Trinity --seqType fq \
--samples_file /data/project_data/RNAseq/assembly/ahud_samples.txt \
--CPU 16 --max_memory 90G --bflyCalculateCPU --min_kmer_cov 2 > run11-21.log 2>&1 &
```

### Assess assembly
```bash
/data/popgen/trinityrnaseq-v2.13.2/util/TrinityStats.pl  ahud_Trinity.fasta

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	130580
Total trinity transcripts:	349516
Percent GC: 35.57

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 4647
	Contig N20: 3149
	Contig N30: 2356
	Contig N40: 1791
	Contig N50: 1356

	Median contig length: 430
	Average contig: 801.91
	Total assembled bases: 280279107


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 4679
	Contig N20: 3115
	Contig N30: 2247
	Contig N40: 1613
	Contig N50: 1057

	Median contig length: 320
	Average contig: 626.33
	Total assembled bases: 81786351
```

Install busco
```bash
conda install -c conda-forge -c bioconda busco=5.2.2

sudo ln -s /data/popgen/busco-5.2.2/bin/busco /usr/local/bin/
```
Test that it works:

```bash
busco -h
```
Yes!

```bash
busco -m transcriptome -i ahud_Trinity.fasta -o BUSCOarthropoda -l arthropoda_odb10
```

```bash
busco -m transcriptome -i ahud_Trinity.fasta -o BUSCOeuk --auto-lineage-euk
```

```bash
cat short_summary.specific.eukaryota_odb10.BUSCOeuk.txt 
# BUSCO version is: 5.2.2 
# The lineage dataset is: eukaryota_odb10 (Creation date: 2020-09-10, number of genomes: 70, number of BUSCOs: 255)
# Summarized benchmarking in BUSCO notation for file /data/project_data/RNAseq/assembly/ahud_Trinity.fasta
# BUSCO was run in mode: transcriptome

	***** Results: *****

	C:97.7%[S:7.5%,D:90.2%],F:2.4%,M:0.1%,n:255	   
	249	Complete BUSCOs (C)			   
	19	Complete and single-copy BUSCOs (S)	   
	230	Complete and duplicated BUSCOs (D)	   
	6	Fragmented BUSCOs (F)			   
	0	Missing BUSCOs (M)			   
	255	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.1
	metaeuk: 5.34c21f2
	```
	
```bash
cat short_summary.specific.arthropoda_odb10.BUSCOarthropoda.txt 
# BUSCO version is: 5.2.2 
# The lineage dataset is: arthropoda_odb10 (Creation date: 2020-09-10, number of genomes: 90, number of BUSCOs: 1013)
# Summarized benchmarking in BUSCO notation for file /data/project_data/RNAseq/assembly/ahud_Trinity.fasta
# BUSCO was run in mode: transcriptome

	***** Results: *****

	C:96.9%[S:7.1%,D:89.8%],F:1.1%,M:2.0%,n:1013	   
	982	Complete BUSCOs (C)			   
	72	Complete and single-copy BUSCOs (S)	   
	910	Complete and duplicated BUSCOs (D)	   
	11	Fragmented BUSCOs (F)			   
	20	Missing BUSCOs (M)			   
	1013	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.1
	metaeuk: 5.34c21f2
```
	
#### BLAST to sprot for annotation

Also examine the number of transcripts that were assembled that appear to be full-length or nearly full-length by comparing to sprot

```bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot
```
 Perform the blast search, reporting only the top alignment:
```bash
blastx -query /data/project_data/RNAseq/assembly/ahud_Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 -evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6
        
```

## map
Make the reference, prep the assembly for mapping, then map using salmon
```
/data/popgen/trinityrnaseq-v2.13.2/util/align_and_estimate_abundance.pl --transcripts /data/project_data/RNAseq/assembly/ahud_Trinity.fasta \
  --est_method salmon \
  --trinity_mode \
  --prep_reference

/data/popgen/trinityrnaseq-v2.13.2/util/align_and_estimate_abundance.pl --transcripts /data/project_data/RNAseq/assembly/ahud_Trinity.fasta \
  --seqType fq --samples_file /data/project_data/RNAseq/assembly/ahud_samples.txt \
  --est_method salmon \
  --output_dir salmon-full-trinity \
  --thread_count 16 \
  --trinity_mode
```

Make file list to point to just the `quant.sf` files
```
/data/popgen/trinityrnaseq-v2.13.2/util/abundance_estimates_to_matrix.pl --est_method salmon \
  --gene_trans_map /data/project_data/RNAseq/assembly/ahud_Trinity.fasta.gene_trans_map \
  --quant_files /data/project_data/RNAseq/mapping/salmon_results_filelist2.txt \
  --name_sample_by_basedir
```
Output combined matrix.
```
/data/popgen/trinityrnaseq-v2.13.2/util/support_scripts/run_TMM_scale_matrix.pl --matrix salmon.isoform.TPM.not_cross_norm > salmon.isoform.TMM.EXPR.matrixCMD: R --no-save --no-restore --no-site-file --no-init-file -q < salmon.isoform.TPM.not_cross_norm.runTMM.R 1>&2 

```
#### make supertranscripts for the reference. 

parse down the transcriptome to only decent genes:

first need to quantify expression:

```bash

#!/bin/bash
######
#
# quantify each sample with salmon
#
#######

# -l A tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.)
# -p 8 says uses 8 threads
# -o indicates the directory and name of output
# -l A means automatically infer library type
# seqbias corrects for random hexamer priming
# gcbias corrects for gcbias, but only when present.

conda activate salmon

for i in $(ls ~/hudsonica_transcriptomics/data/trimmed | grep '.fq.gz' | cut -f 1-3 -d "_"| uniq);
do

    echo "starting sample ${i}"
    #starting with only name of rep. need to pull out files

    read1=$(ls ~/hudsonica_transcriptomics/data/trimmed | grep ${i} | grep '_1.qc.fq.gz')
    read2=$(ls ~/hudsonica_transcriptomics/data/trimmed | grep ${i} | grep '_2.qc.fq.gz')

    salmon quant -i /data/copepods/RNAseq/ahud_index \
        -l A \
         -1 ~/hudsonica_transcriptomics/data/trimmed/${read1} \
         -2 ~/hudsonica_transcriptomics/data/trimmed/${read2} \
         -p 8  \
         --softclip \
         --seqBias \
         --gcBias \
         -o ~/hudsonica_transcriptomics/analysis/salmon/${i}

    echo "sample ${i} done"

done

```

Then filter gene counts

```r

library(tidyr)
library(DESeq2)

dat <- read.table("/data/copepods/hudsonica_multigen_transcriptomics/salmon/salmon.gene.counts.matrix")

dat <- (round(dat, 0))

samples <- colnames(dat)

samples <- gsub("F0", "F00",samples)
samples <- gsub("F2", "F02",samples)
samples <- gsub("F4", "F04",samples)
samples <- substr(samples, 1,11)
colnames(dat) <- samples
# making labels for sample table
id <- separate(data=as.data.frame(samples),col="samples", sep="_", into = c("Population", "Generation", "Replicate"))

sampleTable <- data.frame(
        sampleName = samples,
        Treatment = paste(id$Population, id$Generation, sep="_"))

# assign row names to sample table
rownames(sampleTable) <- samples

# convert to DEseq
dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = sampleTable,
                              design = ~ Treatment)



# remove rows where count < 10 in more than 90% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > (round(ncol(dat)*0.9)), TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)
#[1] 36636

write.table(file="~/hudsonica_genomics/analysis/RNA/expressed.txt", rownames(dds), 
    sep="\t",
    row.names=F, quote=FALSE, col.names=F) 

```

create super transcripts from filtered gene set:

```bash
cd /data/copepods/hudsonica_multigen_transcriptomics

# collapse transcripts to genes
~/bin/trinityrnaseq-Trinity-v2.8.0/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
       --trinity_fasta /data/copepods/hudsonica_multigen_transcriptomics/ahud_Trinity.fasta \
       --out_prefix ahud_super_transcript

cat ahud_super_transcript.fasta | grep '^>' | wc -l
# 130580
# parse to only expressed genes:

grep -Fw -A 1 -f ~/hudsonica_genomics/analysis/RNA/expressed.txt ahud_super_transcript.fasta > ahud_super_transcript.ExpressedOnly.fasta

cat ahud_super_transcript.ExpressedOnly.fasta | grep '^>' | wc -l

~/bin/bwa/bwa index ahud_super_transcript.ExpressedOnly.fasta

```

### map
```bash

my_bwa=~/bin/bwa/bwa
my_samtools=~/bin/samtools-1.6/samtools
bwagenind=/data/copepods/hudsonica_multigen_transcriptomics/ahud_super_transcript.ExpressedOnly.fasta
my_samblstr=~/bin/samblaster/samblaster


cd ~/hudsonica_transcriptomics/data/trimmed

for sample in `ls ~/hudsonica_transcriptomics/data/trimmed | grep '.fq.gz' | cut -f 1-3 -d "_"| uniq`

do

    echo "starting sample ${sample}"
    #starting with only name of rep. need to pull out files

    lane1_1=$(ls ~/hudsonica_transcriptomics/data/trimmed | grep ${sample} | grep '_1')
    lane1_2=$(ls ~/hudsonica_transcriptomics/data/trimmed | grep ${sample} | grep '_2')

    # align lane 1
    echo "starting lane 1 for ${sample}"

    lib1=$(echo $lane1_1 | cut -f 1-3,6 -d "_")

    rg=$(echo \@RG\\tID:$sample\\tPL:Illumina\\tPU:x\\tLB:$lib1\\tSM:$sample)

    $my_bwa mem -t 4 -R $rg $bwagenind $lane1_1 $lane1_2 | \
    $my_samblstr |\
    $my_samtools view -h -u - | \
    $my_samtools sort - -O bam -o ~/hudsonica_genomics/analysis/RNA/aligned/${lib1}.bam


done

```

## call snps

```bash

samtools mpileup  -B --max-depth 50000 --skip-indels -f /data/copepods/hudsonica_multigen_transcriptomics/ahud_super_transcript.ExpressedOnly.fasta AA_F0_Rep1.bam AA_F0_Rep2.bam AA_F0_Rep3.bam AA_F2_Rep2.bam AA_F2_Rep3.bam AA_F4_Rep1.bam AA_F4_Rep2.bam AA_F4_Rep3.bam AA_F11_Rep1.bam AA_F11_Rep2.bam AA_F11_Rep3.bam AH_F0_Rep1.bam AH_F0_Rep2.bam AH_F0_Rep3.bam AH_F2_Rep1.bam AH_F2_Rep2.bam AH_F2_Rep3.bam AH_F4_Rep1.bam AH_F4_Rep2.bam AH_F4_Rep3.bam HA_F0_Rep1.bam HA_F0_Rep2.bam HA_F0_Rep3.bam HA_F2_Rep1.bam HA_F2_Rep2.bam HA_F2_Rep3.bam HA_F4_Rep1.bam HA_F4_Rep2.bam HA_F4_Rep3.bam HH_F0_Rep1.bam HH_F0_Rep2.bam HH_F0_Rep3.bam HH_F4_Rep1.bam  HH_F4_Rep2.bam  HH_F4_Rep3.bam HH_F11_Rep1.bam  HH_F11_Rep2.bam HH_F11_Rep3.bam |\
    java -jar ~/bin/VarScan.v2.3.9.jar mpileup2snp --mpileup 1 --min-coverage 20 --min-reads 2 --min-avg-qual 15 --min-var-freq 0.001 --variants --p-value 0.1 > ~/hudsonica_genomics/analysis/RNA/snp_all_out

```

## filter snps

```bash
cat ~/hudsonica_genomics/analysis/RNA/snp_all_out | wc -l
# 5,843,400

# below here, NR==1 keeps the header, column 10 is SamplesNC

# how many have no samples missing
awk '(NR==1) || ($10 < 1 ) '  ~/hudsonica_genomics/analysis/RNA/snp_all_out | wc -l
# 1,164,229

awk '(NR==1) || ($10 < 1 ) ' ~/hudsonica_genomics/analysis/RNA/snp_all_out > ~/hudsonica_genomics/analysis/RNA/snp_out.nomiss


```


```r


## filtering raw variants

library(stringr)
library(data.table)

dat <- read.table("~/hudsonica_genomics/analysis/RNA/snp_out.nomiss", stringsAsFactors=FALSE, skip=1)
datnm <- read.table("~/hudsonica_genomics/analysis/RNA/snp_out.nomiss", stringsAsFactors=FALSE, nrows=1)

pops <- c(
            "AA_F00_Rep1","AA_F00_Rep2","AA_F00_Rep3", 
            "AA_F02_Rep2","AA_F02_Rep3", 
            "AA_F04_Rep1","AA_F04_Rep2","AA_F04_Rep3", 
            "AA_F11_Rep1","AA_F11_Rep2","AA_F11_Rep3", 
            "AH_F00_Rep1","AH_F00_Rep2","AH_F00_Rep3",  
            "AH_F02_Rep1","AH_F02_Rep2","AH_F02_Rep3", 
            "AH_F04_Rep1","AH_F04_Rep2","AH_F04_Rep3", 
            "HA_F00_Rep1","HA_F00_Rep2","HA_F00_Rep3", 
            "HA_F02_Rep1","HA_F02_Rep2","HA_F02_Rep3",  
            "HA_F04_Rep1","HA_F04_Rep2","HA_F04_Rep3", 
            "HH_F00_Rep1","HH_F00_Rep2","HH_F00_Rep3", 
            "HH_F04_Rep1","HH_F04_Rep2","HH_F04_Rep3", 
            "HH_F11_Rep1","HH_F11_Rep2","HH_F11_Rep3")



colnames(dat) <- c(datnm[1,1:10], pops)

dat2 <- dat

# filter by coverage:

depthout <- as.data.frame(matrix(nrow=nrow(dat2), ncol=length(pops)))
colnames(depthout) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat2[,grep(i_pop, colnames(dat2))]
        depth <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,2])
        # sum up reads
        depthout[,grep(i_pop, colnames(depthout))] <- depth

    }

head(depthout)
colMeans(depthout)
hist(as.matrix(depthout), xlim=c(0, 5000), breaks=50)
upp_depth <- quantile(as.matrix(depthout), c(0.95, 0.975, 0.99))
#   797  1419  2618

hi_cv <- (apply(depthout, 1, function(x) {ifelse((length(which(x > upp_depth[1])) > 0), FALSE, TRUE)}))
sum(hi_cv)
#[1] 1,056,758
dat3 <- dat2[hi_cv,]

#many of these are skewed by indels. only keep reads where depth of actual bialleleic snps > 40
# from the manual: Also, VarScan reports variants on a biallelic basis.
    #That is, for a given SNP call, the "reads1" column is the number of
    #reference-supporting reads (RD), and the "reads2" column is the number of
    #variant-supporting reads (AD).
    #There may be additional reads at that position showing other bases (SNP or indel variants).
    #If these other variants meet the calling criteria, they will be reported in
    #their own line. If not, it may look like you have "missing" reads.
# columns for each call are:
    #consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
keep <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(keep) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])
        # sum up reads
        keep[,grep(i_pop, colnames(keep))] <- (maj+ minor)

    }

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 40)) > 0), FALSE, TRUE)}))
sum(low_cv)
#[1] 478,643
dat4 <- dat3[low_cv,]
nrow(dat4)
dat3 <- dat4
nrow(dat3)

# to get rid of multiallele sites or deletions/insertions
dat4 <- dat3[(which(as.numeric(lapply(strsplit(dat3[,4],","), length)) == 1)),]
nrow(dat4)
# [1] 447,940
dat3 <- dat4

# here calculate allele freqs
# columns for each call are: consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
af <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(af) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])

        # calculate af
        af[,grep(i_pop, colnames(af))] <- maj/(maj+ minor)

    }

# get rid of invariant sites
dat4 <- dat3[(which(rowSums(af) > 0)),]
dat3 <- dat4
af <- af[(which(rowSums(af) > 0)),]
nrow(dat3)
# 447933

af.out <- (cbind(paste(dat3$Chrom, dat3$Position, sep=":"),af))

afct.maf <- (sapply(af,function(x)
          ifelse(x > 0.5, (1-x), x)))

# low maf cut off of < 0.025 in at least 4 groups.
## note that this corresponds to 2 reads at 40x, which seems reasonable.
low_maf <- (apply(afct.maf, 1, function(x) {ifelse((length(which(x > 0.025)) < 4), FALSE, TRUE)}))
sum(low_maf)
# 286139

dat4 <- dat3[low_maf,]
nrow(dat4)
#[1] 286139

af_f <- af[low_maf,]

af.out <- (cbind(paste(dat4$Chrom, dat4$Position, sep=":"),af_f))

colnames(af.out) <- c("SNP", colnames(af_f))

###########
######
# save filtered genotypes
#####
##########

write.table(file="~/hudsonica_genomics/analysis/RNA/filtered_variants.txt", dat4, sep="\t",
              row.names=FALSE, quote=FALSE)

write.table(file="~/hudsonica_genomics/analysis/RNA/filtered_allele_freqs.txt", af.out, sep="\t",
              row.names=FALSE, quote=FALSE)

```


## pca

```r

library(ggplot2)
library(tidyr)

##### snps

af <- read.table("~/hudsonica_genomics/analysis/RNA/filtered_allele_freqs.txt", header=TRUE)

pops <- colnames(af)[2:ncol(af)]

freqs <- t(af[,2:ncol(af)])
colnames(freqs) <-  af$SNP

####
##
## plot pca
##
####

# all gens alone

pcaResult <- prcomp(freqs, scale=TRUE)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(id=pops, Line=substr(pops, 1,2),
        gen=substr(pops, 4,6),
        PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

data <- data %>% separate(id, c("line", "generation", "rep"), remove=F)


#data$Treatment <- gsub("AA", "AM", data$Treatment)

d <- ggplot(data, aes(PC1, PC2, fill=line, shape=line)) +
        geom_point(size=3, color="black") +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
        scale_fill_manual(values=c('#6699CC',"#F2AD00",
                            "#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "OWA"))+
        scale_shape_manual(values=c( 21,22, 23,24), guide = "none")+
        guides(fill=guide_legend(override.aes=list(
            shape=c(21,22, 23,24),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
            fill=c('#6699CC',"#F2AD00",
                            "#00A08A", "#CC3333")),order = 2)) +
        facet_wrap(~gen)


ggsave("~/hudsonica_genomics/figures/RNA/pca.pdf",d, h=5, w=8)

```

