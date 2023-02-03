library(tidyverse)
library(edgeR)
My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

### 1. Preparing the datasets ### 
original_telescope_norm_raw <- read.delim("https://figshare.com/ndownloader/files/39052721")
modified_telescope_norm_raw <- read.delim("https://figshare.com/ndownloader/files/39052739")

### 2. DESEQ ###
library("DESeq2")
TCGA_COAD_paired <- read.delim("https://figshare.com/ndownloader/files/39038870")

rownames(TCGA_COAD_paired) <- TCGA_COAD_paired$submitter_id.samples
TCGA_COAD_paired <- TCGA_COAD_paired[,c(2,3)]
TCGA_COAD_paired <- TCGA_COAD_paired[order(rownames(TCGA_COAD_paired), decreasing=TRUE),]

colnames(original_telescope_norm_raw) <- substr(colnames(original_telescope_norm_raw),1,16)
original_telescope_norm_raw <- original_telescope_norm_raw[rowSums(original_telescope_norm_raw[])>0,]
original_telescope_norm_raw <- original_telescope_norm_raw[,order(colnames(original_telescope_norm_raw), decreasing=TRUE)]

original_telescope_norm_raw <- as.data.frame(lapply(original_telescope_norm_raw, as.integer))
colnames(original_telescope_norm_raw) <- gsub(".","-", colnames(original_telescope_norm_raw), fixed=TRUE)
rownames(original_telescope_norm_raw) <- rownames(original_telescope_norm_raw)
dds <- DESeqDataSetFromMatrix(countData = original_telescope_norm_raw,
                              colData = TCGA_COAD_paired,
                              design = ~ type)
dds <- DESeq(dds)
res <- results(dds)
res <- data.frame(res)
res_org <- na.omit(res)
res_org <- res_org[order(res_org$padj),]
deseq_original <- res_org[which(abs(res_org$log2FoldChange) > 1 & res_org$padj < 0.05),]


colnames(modified_telescope_norm_raw) <- substr(colnames(modified_telescope_norm_raw),1,16)
modified_telescope_norm_raw <- modified_telescope_norm_raw[rowSums(modified_telescope_norm_raw[])>0,]
modified_telescope_norm_raw <- modified_telescope_norm_raw[,order(colnames(modified_telescope_norm_raw), decreasing=TRUE)]

modified_telescope_norm_raw <- as.data.frame(lapply(modified_telescope_norm_raw, as.integer))
colnames(modified_telescope_norm_raw) <- gsub(".","-", colnames(modified_telescope_norm_raw), fixed=TRUE)
rownames(modified_telescope_norm_raw) <- rownames(modified_telescope_norm_raw)
modified_telescope_norm_raw <- modified_telescope_norm_raw[rowSums(modified_telescope_norm_raw[])>0,]

dds <- DESeqDataSetFromMatrix(countData = modified_telescope_norm_raw,
                              colData = TCGA_COAD_paired,
                              design = ~ type)
dds <- DESeq(dds)
res <- results(dds)
res <- data.frame(res)
res_mod <- na.omit(res)
res_mod <- res_mod[order(res_mod$padj),]
deseq_modified <- res_mod[which(abs(res_mod$log2FoldChange) > 1 & res_mod$padj < 0.05),]

### only select DETEs with 80%
original_telescope_norm_raw <- original_telescope_norm_raw[rownames(original_telescope_norm_raw) %in% rownames(deseq_original),]
modified_telescope_norm_raw <- modified_telescope_norm_raw[rownames(modified_telescope_norm_raw) %in% rownames(deseq_modified),]
threshold <- ncol(original_telescope_norm_raw) * 0.2
original_telescope_norm_raw <- original_telescope_norm_raw[rowSums(original_telescope_norm_raw == 0) <= threshold, ]
modified_telescope_norm_raw <- modified_telescope_norm_raw[rowSums(modified_telescope_norm_raw == 0) <= threshold, ]

## Filtered out intergenic only 
hg38_gencode_rmsk_indi_loc_annotate <- read.delim("https://figshare.com/ndownloader/files/39052802", header=FALSE)
original_telescope_norm_raw <- original_telescope_norm_raw[rownames(original_telescope_norm_raw) %in% hg38_gencode_rmsk_indi_loc_annotate$V4[hg38_gencode_rmsk_indi_loc_annotate$V5 =="intergenic"],]
modified_telescope_norm_raw <- modified_telescope_norm_raw[rownames(modified_telescope_norm_raw) %in% hg38_gencode_rmsk_indi_loc_annotate$V4[hg38_gencode_rmsk_indi_loc_annotate$V5 =="intergenic"],]

## 3. Correlating with DNA methylation values
library(GenomicRanges)
HM450.hg38.manifest <- read.delim("https://figshare.com/ndownloader/files/39038831")

HM450.hg38.manifest <- HM450.hg38.manifest[,c(1:5)]
HM450.hg38.manifest <- na.omit(HM450.hg38.manifest)
HM450.hg38.manifest_GR <- GRanges(seqnames=Rle(HM450.hg38.manifest$CpG_chrm), ranges=IRanges(start = HM450.hg38.manifest$CpG_beg, end = HM450.hg38.manifest$CpG_end), names=HM450.hg38.manifest$probeID)

hg38_gencode_rmsk_indi <- read.delim("https://figshare.com/ndownloader/files/39038720", header=FALSE)
hg38_gencode_rmsk_indi$label <- paste(hg38_gencode_rmsk_indi$V1, paste(hg38_gencode_rmsk_indi$V4, hg38_gencode_rmsk_indi$V5, sep="|"), sep="|")
hg38_gencode_rmsk_indi$indi_label <- unlist(sapply(strsplit(hg38_gencode_rmsk_indi$V9, " ", fixed=TRUE), function(x) x[2], simplify=FALSE))
hg38_gencode_rmsk_indi$indi_label_final <- unlist(sapply(strsplit(hg38_gencode_rmsk_indi$indi_label, ";", fixed=TRUE), function(x) x[1], simplify=FALSE))
hg38_gencode_rmsk_indi <- data.frame(hg38_gencode_rmsk_indi)

deseq_modified_up <- deseq_modified[deseq_modified$log2FoldChange < 0,]
deseq_original_up <- deseq_original[deseq_original$log2FoldChange < 0,]

intersected_te <- intersect(rownames(deseq_modified_up), rownames(deseq_original_up))
only_mod <- rownames(deseq_modified_up)[!rownames(deseq_modified_up) %in% rownames(deseq_original_up)]

TCGA_COAD_methylation <- read.delim("https://figshare.com/ndownloader/files/39038867")
rownames(TCGA_COAD_methylation) <- TCGA_COAD_methylation$sample
TCGA_COAD_methylation <- TCGA_COAD_methylation[,colnames(TCGA_COAD_methylation) %in% substr(colnames(original_telescope_norm_raw),1,15)]
TCGA_COAD_methylation <- na.omit(TCGA_COAD_methylation)
TCGA_COAD_methylation <- TCGA_COAD_methylation[,order(colnames(TCGA_COAD_methylation))]

## performing correlation between DNA methylation and DETEs among intersected TEs
all_te <- union(intersected_te, only_mod)
all_te_df <- hg38_gencode_rmsk_indi[match(all_te,hg38_gencode_rmsk_indi$indi_label_final),c(1,4,5,12)]
all_te_df$V4 <- as.numeric(all_te_df$V4) - 500
all_te_df$V5 <- as.numeric(all_te_df$V5) + 500
all_te_df_GR <- GRanges(seqnames = Rle(all_te_df$V1), ranges=IRanges(start=all_te_df$V4, end=all_te_df$V5), names=all_te_df$indi_label_final)

temp_all <- findOverlaps(HM450.hg38.manifest_GR, all_te_df_GR)
temp_all_df <- data.frame(HM450.hg38.manifest_GR[queryHits(temp_all)])
temp_all_df2 <- data.frame(all_te_df_GR[subjectHits(temp_all)])
temp_all_combined <- cbind(temp_all_df, temp_all_df2)

modified_telescope_norm_raw <- modified_telescope_norm_raw[,substr(colnames(modified_telescope_norm_raw),1,15) %in% colnames(TCGA_COAD_methylation)]

## total 559 temp_all_combined
temp_all_combined$group_label <- paste(temp_all_combined[,6], temp_all_combined[,12], sep="_")
mod <- na.omit(temp_all_combined)
mod$label <- unlist(sapply(strsplit(mod[,6],".", fixed = TRUE), function(x) x[1], simplify=FALSE))

## 329 CpG island found to have 
mod <- mod[mod$label %in% rownames(TCGA_COAD_methylation),]
mod$cor <- 0
for(i in 1:nrow(mod)){
  mod$cor[i] <- cor.test(as.numeric(TCGA_COAD_methylation[match(mod$label[i],rownames(TCGA_COAD_methylation)),]),as.numeric(modified_telescope_norm_raw[match(mod[i,12],rownames(modified_telescope_norm_raw)),]),method="spearman")$estimate
}

mod$label <- "common"
mod$label[mod[,12] %in% only_mod] <- "only in lasTEq"
mod_nec <- mod[,c(13:15)]
p_cor <- ggplot(mod_nec, aes(x=label, y=cor, fill=label)) + geom_boxplot(width=0.2, lwd=1) +geom_quasirandom(method = "pseudorandom", color = "#7E6148FF", alpha=0.8, size=0.5)
p_cor <- p_cor + My_Theme + xlab("TEs")+ylab("cor(TE and beta(T-N))") + ggtitle("Cor.value of intergenic TE and Methylation\n(500bp +/- TEs)")
p_cor <- p_cor + scale_fill_manual(values=c("#44AA99","#CC6677"))