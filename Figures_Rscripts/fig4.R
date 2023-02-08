My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
### 1. Age Analysis ###
combined_dataset_ERR6937788 <- read.delim("https://figshare.com/ndownloader/files/39038408", comment.char="#")
combined_dataset_ERR6937791 <- read.delim("https://figshare.com/ndownloader/files/39038402", comment.char="#")
combined_dataset_ERR6937794 <- read.delim("https://figshare.com/ndownloader/files/39038390", comment.char="#")
combined_dataset_ERR6937797 <- read.delim("https://figshare.com/ndownloader/files/39038396", comment.char="#")
combined_dataset_ERR6937800 <- read.delim("https://figshare.com/ndownloader/files/39038411", comment.char="#")
combined_dataset_HCT116 <- read.delim("https://figshare.com/ndownloader/files/39038414", comment.char="#")

age_process <- function(combined_age){
  combined_age$lasteq[combined_age$lasteq < 0.01] <- 0
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

### Correlation Plot ###

library(corrplot)
cell_line_name <- "own_Hct116"
combined_dataset <- combined_dataset_HCT116

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
colnames(cor_df) <- c("ownHCT116","A549","HCT116","HepG2","K562","MCF-7")
rownames(cor_df) <- c("Telescope","lasTEq")
cor_df[1,] <- c(0.32, 0.57, 0.57,0.46,0.44,0.57)
cor_df[2,] <- c(0.37, 0.58, 0.61,0.47,0.45,0.59)
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


### bed file extraction for histone ###
histone_own_HCT116 <- combined_dataset_HCT116[,c(1,2,4,5)]
histone_own_HCT116$diff <- histone_own_HCT116$lasteq - histone_own_HCT116$short
histone_own_HCT116$label_diff <- "same"
histone_own_HCT116$label_diff[histone_own_HCT116$diff > 1] <- "high_mixed"
histone_own_HCT116$label_diff[histone_own_HCT116$diff < -1] <- "high_short"
histone_own_HCT116$chr <- unlist(sapply(strsplit(histone_own_HCT116$label, "|", fixed=TRUE), function(x) x[1], simplify=FALSE))
histone_own_HCT116$str <- unlist(sapply(strsplit(histone_own_HCT116$label, "|", fixed=TRUE), function(x) x[2], simplify=FALSE))
histone_own_HCT116$end <- unlist(sapply(strsplit(histone_own_HCT116$label, "|", fixed=TRUE), function(x) x[3], simplify=FALSE))

## save into bed files 
bigger1_Corrected_in_modified_telescope <- histone_own_HCT116[histone_own_HCT116$label_diff == "high_mixed", 7:9]
bigger1_Incorrectly_measured_in_original_telescope <- histone_own_HCT116[histone_own_HCT116$label_diff == "high_short", 7:9]

### 4. Histone Heatmap Plot ###

## histone heatmap
histone_df <- data.frame(rep(0, 5))
for(i in 1:11){
  histone_df <- cbind(histone_df, rep(0,5))
}
rownames(histone_df) <- c("H3K27ac","H3K27me3","H3K4me3","H3K9ac","H3K9me3")
colnames(histone_df) <- c("correct_ownHCT116","incorrect_ownHCT116","correct_A549","incorrect_A549","correct_HCT116","incorrect_HCT116","correct_HepG2","incorrect_HepG2","correct_K562","incorrect_K562","correct_MCF-7","incorrect_MCF-7")

histone_df[,1] <- c(1.97431, 0.361352, 1.79114, 1.74377, 0.437465)
histone_df[,2] <- c(1.55554, 0.366998, 1.31739, 1.35212, 0.446571)

histone_df[,3] <- c(0.353323, 0.0345828, 0.532635, 0.337217, 0.0954647)
histone_df[,4] <- c(0.38188, 0.0559413, 0.350027, 0.335743, 0.146679)

histone_df[,5] <- c(0.239359, 0.0498001, 0.146767, 0.157627, 0.0597366)
histone_df[,6] <- c(0.251522, 0.0684533, 0.195611, 0.185631, 0.0964365)

histone_df[,7] <- c(0.242406, 0.040903, 0.114888, 0.192506, 0.0673679)
histone_df[,8] <- c(0.148777, 0.0540486, 0.125958, 0.179852, 0.116947)

histone_df[,9] <- c(0.113632, 0.0494671, 0.126207, 0.107562, 0.1461)
histone_df[,10] <- c(0.18671, 0.0871798, 0.177329, 0.147815, 0.437465)

histone_df[,11] <- c(0.196137, 0.0499715, 0.238204, 0.255244, 0.0783445)
histone_df[,12] <- c(0.14549, 0.0930822, 0.200189, 0.164347, 0.141475)

histone_df <- histone_df[,c(3:12)]

histone_df <- cbind("label"=rownames(histone_df),histone_df)
histone_df <- histone_df[,c(1,3,5,7,9,11,2,4,6,8,10)]
colnames(histone_df) <- factor(colnames(histone_df), levels=colnames(histone_df))

histone_df <- melt(histone_df)
histone_df$variable <- as.character(histone_df$variable)
histone_df$cell <- unlist(sapply(strsplit(histone_df$variable,"_", fixed = TRUE), function(x) x[2], simplify=FALSE))
histone_df$group <- unlist(sapply(strsplit(histone_df$variable,"_", fixed = TRUE), function(x) x[1], simplify=FALSE))
histone_df$group[histone_df$group=="correct"] <- "Originally_Correctly\nMeasured"
histone_df$group[histone_df$group=="incorrect"] <- "Incorrectly\nMeasured"


#H3K27ac
histone_df_nec <- histone_df[which(histone_df$label =="H3K27ac"),]
histone_heatmap <- ggplot(histone_df_nec, aes(x = group, y = cell, fill=value)) + geom_tile(color = "white", lwd = 1.5,linetype = 1)+scale_x_discrete(labels=c("Incorrectly\nMeasured", "Correctly\nMeasured")) 
histone_heatmap <- histone_heatmap +My_Theme + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
histone_heatmap <- histone_heatmap + labs(fill = "Mean\n(Coverage)") + scale_fill_gradient2(limits=c(0, 0.4),breaks=c(0, 0.2,0.4), labels=c(0,0.2, 0.4) ,low = "white",
                                                                                            high = "#DC0000FF",
                                                                                            guide = "colorbar") 
histone_heatmap <- histone_heatmap + ggtitle("H3K27ac")

#H3K27me3
histone_df_nec <- histone_df[which(histone_df$label =="H3K27me3"),]
range(histone_df_nec$value)
histone_heatmap <- ggplot(histone_df_nec, aes(x = group, y = cell, fill=value)) + geom_tile(color = "white", lwd = 1.5,linetype = 1)+scale_x_discrete(labels=c("Incorrectly\nMeasured", "Correctly\nMeasured")) 
histone_heatmap <- histone_heatmap +My_Theme + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
histone_heatmap <- histone_heatmap + labs(fill = "Mean\n(Coverage)") + scale_fill_gradient(limits=c(0.02, 0.1),breaks=c(0.02, 0.06,0.1), labels=c(0.02,0.06, 0.1) ,low = "#3C5488FF",
                                                                                           high = "white",
                                                                                           guide = "colorbar") 
histone_heatmap <- histone_heatmap + ggtitle("H3K27me3")

#H3K4me3
histone_df_nec <- histone_df[which(histone_df$label =="H3K4me3"),]
range(histone_df_nec$value)
histone_heatmap <- ggplot(histone_df_nec, aes(x = group, y = cell, fill=value)) + geom_tile(color = "white", lwd = 1.5,linetype = 1)+scale_x_discrete(labels=c("Incorrectly\nMeasured", "Correctly\nMeasured")) 
histone_heatmap <- histone_heatmap +My_Theme + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
histone_heatmap <- histone_heatmap + labs(fill = "Mean\n(Coverage)") + scale_fill_gradient(limits=c(0, 0.6),breaks=c(0, 0.3,0.6), labels=c(0,0.3, 0.6) ,low = "white",
                                                                                           high = "#DC0000FF",
                                                                                           guide = "colorbar") 
histone_heatmap <- histone_heatmap + ggtitle("H3K4me3")

#H3K9ac
histone_df_nec <- histone_df[which(histone_df$label =="H3K9ac"),]
range(histone_df_nec$value)
histone_heatmap <- ggplot(histone_df_nec, aes(x = group, y = cell, fill=value)) + geom_tile(color = "white", lwd = 1.5,linetype = 1)+scale_x_discrete(labels=c("Incorrectly\nMeasured", "Correctly\nMeasured")) 
histone_heatmap <- histone_heatmap +My_Theme + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
histone_heatmap <- histone_heatmap + labs(fill = "Mean\n(Coverage)") + scale_fill_gradient(limits=c(0, 0.4),breaks=c(0, 0.2,0.4), labels=c(0,0.2, 0.4) ,low = "white",
                                                                                           high = "#DC0000FF",
                                                                                           guide = "colorbar") 
histone_heatmap <- histone_heatmap + ggtitle("H3K9ac")

#H3K9me3
histone_df_nec <- histone_df[which(histone_df$label =="H3K9me3"),]
range(histone_df_nec$value)
histone_heatmap <- ggplot(histone_df_nec, aes(x = group, y = cell, fill=value)) + geom_tile(color = "white", lwd = 1.5,linetype = 1)+scale_x_discrete(labels=c("Incorrectly\nMeasured", "Correctly\nMeasured")) 
histone_heatmap <- histone_heatmap +My_Theme + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
histone_heatmap <- histone_heatmap + labs(fill = "Mean\n(Coverage)") + scale_fill_gradient(limits=c(0, 0.5),breaks=c(0, 0.25,0.5), labels=c(0,0.25, 0.5) ,low = "#3C5488FF",
                                                                                           high = "white",
                                                                                           guide = "colorbar") 
histone_heatmap <- histone_heatmap + ggtitle("H3K9me3")
