My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
### 1. Age Analysis ###
combined_dataset_ERR6937788 <- read.delim("https://figshare.com/ndownloader/files/39359537", comment.char="#")
combined_dataset_ERR6937791 <- read.delim("https://figshare.com/ndownloader/files/39359543", comment.char="#")
combined_dataset_ERR6937794 <- read.delim("https://figshare.com/ndownloader/files/39359549", comment.char="#")
combined_dataset_ERR6937797 <- read.delim("https://figshare.com/ndownloader/files/39359555", comment.char="#")
combined_dataset_ERR6937800 <- read.delim("https://figshare.com/ndownloader/files/39359558", comment.char="#")
combined_dataset_HCT116 <- read.delim("https://figshare.com/ndownloader/files/39359567", comment.char="#")

age_process <- function(combined_age){
  combined_age <- combined_age[which(combined_age$short != 0 | combined_age$lasteq != 0),]
  combined_age$diff_raw <- as.numeric(combined_age$short) - as.numeric(combined_age$lasteq)
  combined_age$corrected_by_modified_telescope <- "not captured"
  combined_age$corrected_by_modified_telescope[which(combined_age$diff_raw != 0)] <- "captured"
  combined_age$originally_correct <- "not captured"
  combined_age$originally_correct[which(combined_age$diff_raw == 0)] <- "captured"
  combined_age <- combined_age[,c(2,7,8)]
  colnames(combined_age)[2] <- "TEs Corrected by lasTEq"
  colnames(combined_age)[3] <- "TEs Originally Correct"
  combined_age <- melt(combined_age,id = "label")
  combined_age <- combined_age[which(combined_age$value != 0 & combined_age$value != "not captured"),]
  return(combined_age)
}


combined_dataset_ERR6937788_age <- age_process(combined_dataset_ERR6937788)
combined_dataset_ERR6937791_age <- age_process(combined_dataset_ERR6937791)
combined_dataset_ERR6937794_age <- age_process(combined_dataset_ERR6937794)
combined_dataset_ERR6937797_age <- age_process(combined_dataset_ERR6937797)
combined_dataset_ERR6937800_age <- age_process(combined_dataset_ERR6937800)
combined_dataset_HCT116_age <- age_process(combined_dataset_HCT116)
combined_dataset_ERR6937788_age <- cbind("group"=rep("SG-NEX A549"), combined_dataset_ERR6937788_age)
combined_dataset_ERR6937791_age <- cbind("group"=rep("SG-NEX HepG2"), combined_dataset_ERR6937791_age)
combined_dataset_ERR6937794_age <- cbind("group"=rep("SG-NEX HCT116"), combined_dataset_ERR6937794_age)
combined_dataset_ERR6937797_age <- cbind("group"=rep("SG-NEX K562"), combined_dataset_ERR6937797_age)
combined_dataset_ERR6937800_age <- cbind("group"=rep("SG-NEX MCF-7"), combined_dataset_ERR6937800_age)
combined_dataset_HCT116_age <- cbind("group"=rep("Own HCT116"), combined_dataset_HCT116_age)

age_combined <- rbind(combined_dataset_ERR6937788_age,combined_dataset_ERR6937791_age)
age_combined <- rbind(age_combined,combined_dataset_ERR6937794_age)
age_combined <- rbind(age_combined,combined_dataset_ERR6937797_age)
age_combined <- rbind(age_combined,combined_dataset_ERR6937800_age)
age_combined <- rbind(age_combined,combined_dataset_HCT116_age)

hg38_repeatmasker <- read.delim("https://figshare.com/ndownloader/files/39038723", header=FALSE, comment.char="#")
hg38_repeatmasker$V7 <- hg38_repeatmasker$V7+1
hg38_repeatmasker$label <- paste(hg38_repeatmasker$V6, paste(hg38_repeatmasker$V7, hg38_repeatmasker$V8, sep="|"), sep="|")
age_combined$millidiv <- hg38_repeatmasker[match(age_combined$label, hg38_repeatmasker$label),3]

age_combined <- na.omit(age_combined)
age_combined[is.na(age_combined)] <- 0
stats <- compare_means(millidiv ~ variable, group.by = "group", data = age_combined, method = "t.test")
give.n <- age_combined %>% group_by(group, variable) %>% summarise(Freq=n())

age_combined_plot <- ggplot(age_combined, aes(x=variable, y=millidiv, fill=variable)) + geom_boxplot(width=0.6, lwd=1)+facet_wrap(~group,nrow = 1)
age_combined_plot <- age_combined_plot + My_Theme + scale_fill_npg() + scale_y_continuous(limits=c(0,300))
age_combined_plot <- age_combined_plot + theme(axis.title.x=element_blank(),
                                               axis.text.x=element_blank(),
                                               axis.ticks.x=element_blank(),strip.background = element_blank(),
                                               strip.text.x = element_blank()) + ylab("Age(milliDiv)") +stat_compare_means(label =  "p.signif", label.x = 1.5)
age_combined_plot <- age_combined_plot + ggtitle("Age Distribution of TEs")

ggsave(file="/storage/jwlab/sandy/Figures/final/age.pdf", plot=age_combined_plot, bg = 'white', width = 12, height = 8, units = 'cm', dpi = 600)

### Correlation Plot ###

library(corrplot)
cell_line_name <- "own_Hct116"
combined_dataset <- combined_dataset_HCT116

cell_line_name <- "A549"
combined_dataset <- combined_dataset_ERR6937788

cell_line_name <- "HepG2"
combined_dataset <- combined_dataset_ERR6937791

cell_line_name <- "Hct116"
combined_dataset <- combined_dataset_ERR6937794

cell_line_name <- "K562"
combined_dataset <- combined_dataset_ERR6937797

cell_line_name <- "MCF7"
combined_dataset <- combined_dataset_ERR6937800

correlation <- function(combined_dataset,cell_line_name){
  combined_dataset_cor <- combined_dataset[,c(1,2,3:5)]
  combined_dataset_cor$str <- unlist(sapply(strsplit(combined_dataset_cor$label,"|", fixed = TRUE), function(x) x[2], simplify=FALSE))
  combined_dataset_cor$end <- unlist(sapply(strsplit(combined_dataset_cor$label,"|", fixed = TRUE), function(x) x[3], simplify=FALSE))
  combined_dataset_cor$length <- as.numeric(combined_dataset_cor$end) - as.numeric(combined_dataset_cor$str)+1
  
  combined_dataset_cor$lasteq <- round(combined_dataset_cor$lasteq)
  
  ### change to TPM ### use sum of coding genes as library size
  long_read_fc_coding_gene_cleanup <- read.delim("https://figshare.com/ndownloader/files/39037997", comment.char="#")
  short_read_fc_coding_gene_cleanup <- read.delim("https://figshare.com/ndownloader/files/39149033", comment.char="#")
  short_read_fc_coding_gene_cleanup$rate <- as.numeric(short_read_fc_coding_gene_cleanup$X.storage2.jwlab.sandy.HCT116.HCT116_14_12.final_star.HCT116_14_12Aligned.sortedByCoord.out.bam) / as.numeric(short_read_fc_coding_gene_cleanup$Length)
  long_read_fc_coding_gene_cleanup$rate <- as.numeric(long_read_fc_coding_gene_cleanup$X.storage.jwlab.sandy.HCT116.minimap2_genome.HCT116_mRNA_long.bam) / as.numeric(long_read_fc_coding_gene_cleanup$Length)
  
  library_size_short = sum(short_read_fc_coding_gene_cleanup$rate)
  library_size_long = sum(long_read_fc_coding_gene_cleanup$rate)
  
  tpm3 <- function(counts,len) {
    x <- as.numeric(counts)/as.numeric(len)
    return(t(t(x)*1e6))
  }
  
  combined_dataset_cor$longTPM <- (tpm3(as.numeric(combined_dataset_cor$long), as.numeric(combined_dataset_cor$length)))/library_size_long
  combined_dataset_cor$shortTPM <- (tpm3(as.numeric(combined_dataset_cor$short), as.numeric(combined_dataset_cor$length)))/library_size_short
  combined_dataset_cor$lasteqTPM <- (tpm3(as.numeric(combined_dataset_cor$lasteq), as.numeric(combined_dataset_cor$length)))/library_size_short
  
  combined_dataset_cor$diff <- as.numeric(combined_dataset_cor$shortTPM) - as.numeric(combined_dataset_cor$lasteqTPM)
  combined_dataset_cor$diff_label <- "Same"
  combined_dataset_cor[which(combined_dataset_cor$diff != 0),"diff_label"] <- "Changed"
  combined_dataset_cor_common <- combined_dataset_cor[which(abs(combined_dataset_cor$diff)> 1e-5),]
  combined_dataset_cor_common <- combined_dataset_cor_common[,c(1,9:11)]
  combined_dataset_cor_common <- na.omit(combined_dataset_cor_common)
  rownames(combined_dataset_cor_common) <- combined_dataset_cor_common[,1]
  combined_dataset_cor_common <- combined_dataset_cor_common[,2:4]
  
  combined_dataset_cor_common <- combined_dataset_cor_common[!(apply(combined_dataset_cor_common, 1, function(y) any(y == 0))),]
  combined_dataset_cor_common <- log(combined_dataset_cor_common)
  p1 <- ggscatter(combined_dataset_cor_common[,c(1,2)], y = "longTPM", x = "shortTPM",
                  add = "reg.line", color = "#7E6148B2",
                  add.params = list(color = "#3C5488B2", fill = "#8491B4B2", size=4),
                  conf.int = TRUE)+My_Theme
  p1 <- p1 + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size = 4)
  p1 <- p1 + ggtitle("Telescope") + xlab("log(TPM by Telescope)") + ylab("log(TPM by Long Read)")

  
  p2 <- ggscatter(combined_dataset_cor_common[,c(1,3)], y = "longTPM", x = "lasteqTPM",
                  add = "reg.line", color = "#7E6148B2",
                  add.params = list(color = "#DC0000FF", fill = "#8491B4B2", size=4),
                  conf.int = TRUE)+My_Theme
  p2 <- p2 + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size = 4)
  p2 <- p2 + ggtitle("lasTEq") + xlab("log(TPM by lasTEq)") + ylab("log(TPM by Long Read)")
}

### 3. Correlation Bubble Plot ###

### heatmap for correlation
cor_df <- data.frame(rep(0, 2))
for(i in 1:5){
  cor_df <- cbind(cor_df, rep(0,2))
}
colnames(cor_df) <- c("inhouse","A549","HCT116","HepG2","K562","MCF-7")
rownames(cor_df) <- c("Telescope","lasTEq")
cor_df[1,] <- c(0.37,0.53, 0.5,0.46,0.4,0.55)
cor_df[2,] <- c(0.39,0.57, 0.58,0.5,0.44,0.59)
cor_df <- cbind("label"=rownames(cor_df),cor_df)
cor_df$label <- factor(cor_df$label, levels=c("Telescope", "lasTEq"))
cor_df <- cor_df[,c(1,3:7)]
cor_df <- melt(cor_df)
cor_df$log_value <- log10(cor_df$value)+10
cor_df_heatmap <- ggplot(cor_df, aes(x =label, y = variable, color=label,size=log_value)) +geom_point()+scale_size(range = c(10,16))
cor_df_heatmap <- cor_df_heatmap + My_Theme
cor_df_heatmap <- cor_df_heatmap + theme(strip.background = element_blank(), strip.text.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="bottom")
cor_df_heatmap <- cor_df_heatmap + scale_color_manual(values=c("#3C5488FF","#DC0000FF"))+geom_text(aes(label=value), size=4, color="white",fontface = "bold")
cor_df_heatmap <- cor_df_heatmap + ggtitle("Correlation with Long Read")+theme(legend.position="none")

### 4. Histone bed file ###
### bed file extraction for histone ###
histone_own_HCT116 <- combined_dataset_HCT116[,c(1,2,4,5)]
histone_own_HCT116$diff <- histone_own_HCT116$lasteq - histone_own_HCT116$short
histone_own_HCT116$label_diff <- "same"
histone_own_HCT116$label_diff[histone_own_HCT116$diff > 1] <- "high_mixed"
histone_own_HCT116$label_diff[histone_own_HCT116$diff < -1] <- "high_short"
histone_own_HCT116$chr <- unlist(sapply(strsplit(histone_own_HCT116$label, "|", fixed=TRUE), function(x) x[1], simplify=FALSE))
histone_own_HCT116$str <- unlist(sapply(strsplit(histone_own_HCT116$label, "|", fixed=TRUE), function(x) x[2], simplify=FALSE))
histone_own_HCT116$end <- unlist(sapply(strsplit(histone_own_HCT116$label, "|", fixed=TRUE), function(x) x[3], simplify=FALSE))

#write.table(histone_own_HCT116[histone_own_HCT116$label_diff == "high_mixed", 7:9], "bigger1_Corrected_in_modified_telescope.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
#write.table(histone_own_HCT116[histone_own_HCT116$label_diff == "high_short", 7:9], "bigger1_Incorrectly_measured_in_original_telescope.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
