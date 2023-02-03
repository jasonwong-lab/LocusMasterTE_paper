##### Please assign right path to the files accordingly #####

### 1. Processing Long-read ###
## Prepare data
$fa=genome_fasta_file_hg38
$long_fq=long_read_fastq
$te_gtf=TE_individual_gtf

# Minimap2
minimap2 -ax splice $fa $long_fq > processed.sam
samtools view -bS processed.sam > processed.bam

# featureCounts
featureCounts -T 20 -R BAM --verbose -O -L -M -a $te_gtf processed.bam

### 2. Processing Short-read ###
$gene_gtf=Gene_gtf_from_Gencode
$short_fq1=short_read_fastq_pair1
$short_fq2=short_read_fastq_pair2

## star commands(STAR-2.7.9a Version)
# star index
STAR --genomeChrBinNbits 15 --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles $fa --sjdbGTFfile $gene_gtf --limitGenomeGenerateRAM 3762779324690
# run star
STAR --genomeDir STAR_index --runThreadN 15 --outFilterMultimapNmax 500 --readFilesIn $short_fq1 $short_fq2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

### 3. Generating Simulated Short-read ###
$gtf_combined=Combination_TE+Gene_GTF
# Below link is an actual link to TPM file that our study used
$TPM=https://figshare.com/ndownloader/files/39053627
spankisim_transcripts -g $gtf_combined -f $fa -bp 76 -frag 200 -ends 2 -t $TPM