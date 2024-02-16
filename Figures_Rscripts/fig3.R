My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

## --- Simulated Reads --- ##
## Figure 3 A-E ##
combined_extract <- read.delim("https://figshare.com/ndownloader/files/44176667")

f1mat <- data.frame(rep(0,3))
rownames(f1mat) <- c("Precision","Recall","F1")
for(i in c(1:6)){
  f1mat <- cbind(f1mat,rep(0,3))
}
colnames(f1mat) <- colnames(combined_extract)[3:ncol(combined_extract)]
## do it for 6 different quantification methods
for(idx in c(3:ncol(combined_extract))){
  color_temp <- "#91D1C2FF"
  scatterdf <- combined_extract[,c(1:2,idx)]
  colnames(scatterdf)[3] <- "Short"
  scatterdf_conf <- scatterdf
  
  ## confusion matrix
  scatterdf_conf$label <- "TP"
  scatterdf_conf$label[scatterdf_conf$Long != 0 & scatterdf_conf$Short ==0] <- "FN"
  scatterdf_conf$label[scatterdf_conf$Long == 0 & scatterdf_conf$Short !=0] <- "FP"
  scatterdf_conf$label[scatterdf_conf$Short == 0 & scatterdf_conf$Long ==0] <- "TN"
  
  confusion_mat_temp <- scatterdf_conf %>% dplyr::group_by(label) %>% dplyr::summarise(Freq=n())
  confusion_mat_temp <- data.frame(confusion_mat_temp)
  confusion_mat <- data.frame(label=c("TP","FN","FP","TN"))
  confusion_mat$Freq <- unlist(confusion_mat_temp[match(confusion_mat$label, confusion_mat_temp$label),2])
  confusion_mat[is.na(confusion_mat)] <- 0
  confusion_mat
  precision <- as.numeric(confusion_mat[confusion_mat$label=="TP",2]) / (as.numeric(confusion_mat[confusion_mat$label=="TP",2])+as.numeric(confusion_mat[confusion_mat$label=="FP",2]))
  recall <- as.numeric(confusion_mat[confusion_mat$label=="TP",2]) / (as.numeric(confusion_mat[confusion_mat$label=="TP",2])+as.numeric(confusion_mat[confusion_mat$label=="FN",2]))
  f1 <- (2*precision*recall)/(precision+recall)
  
  f1mat[,idx-2] <- c(precision, recall, f1)
}

graph_list <- data.frame(c(3:9))
graph_list$color <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
graph_list$title <- c("Telescope","LocusMasterTE", "SQuIRE","Unique Counts","featureCounts","RSEM","SalmonTE")
graph_list$file_title <- c("Telescope","LocusMasterTE","SQuIRE","Unique Counts","featureCounts","RSEM","SalmonTE")

for(j in c(1:nrow(graph_list))){
  idx <- graph_list$c.3.9.[j]
  color_temp <- graph_list$color[j]
  scatterdf <- combined_extract[,c(1:2,idx)]
  colnames(scatterdf)[3] <- "Short"
  scatterdf$label <- "TP"
  scatterdf$label[scatterdf$Long != 0 & scatterdf$Short ==0] <- "FN"
  scatterdf$label[scatterdf$Long == 0 & scatterdf$Short !=0] <- "FP"
  scatterdf$label[scatterdf$Short == 0 & scatterdf$Long ==0] <- "TN"
  print(graph_list$title[j])
  print(length(scatterdf$label[scatterdf$label=="TP"]))
  print(length(scatterdf$label[scatterdf$label=="TN"]))
  print(length(scatterdf$label[scatterdf$label=="FP"]))
  print(length(scatterdf$label[scatterdf$label=="FN"]))
  scatterdf <- scatterdf[scatte1rdf$label !="TN",]
  scatterdf[,2:3] <- log(scatterdf[,2:3]+1e-7)
  scatterdf2 <- scatterdf[scatterdf$label == "TP",]
  combined_scatter <- ggplot(data=scatterdf, aes(x=Long, y=Short)) + geom_pointdensity(color=color_temp, alpha=0.5)
  combined_scatter <- combined_scatter + geom_hline(yintercept=-5, linetype="dashed", color = "black")
  combined_scatter <- combined_scatter + geom_vline(xintercept=-5, linetype="dashed", color = "black")
  combined_scatter <- combined_scatter + geom_smooth(data=scatterdf2, method = "lm", se = FALSE, fullrange=TRUE,linetype="dashed", color="black")
  combined_scatter <- combined_scatter + xlab("log(TPM by Long Read)") + ylab("log(TPM by Short Read)") + My_Theme
  combined_scatter <- combined_scatter + ggtitle(graph_list$title[j])
}


### precision vs recall graph ###
## precision = TP/TP+FP
## recall = TP/TP+FN
## f1 score = 2*(precision * recall)/(precision+recall)
f1mat <- data.frame(t(f1mat))
f1mat_plot <- cbind("label"=rownames(f1mat),f1mat)
f1mat_plot$label <- c("Telescope","LocusMasterTE", "SQuIRE","Unique Counts","RSEM","SalmonTE","featureCounts")
f1mat_plot$label <- factor(f1mat_plot$label, levels=c("Telescope","LocusMasterTE", "SQuIRE","Unique Counts","featureCounts","RSEM","SalmonTE"))

f1mat_precision <- ggplot(aes(x = Recall, y = Precision,label=label), data = f1mat_plot) +geom_point(aes(fill=label), shape=21, size= 6, color="black")
f1mat_precision <- f1mat_precision +scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF"))
f1mat_precision <- f1mat_precision + ggtitle("Precision vs Recall")

f1mat_precision <- f1mat_precision + theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="blank",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                                           axis.title.x = element_text(colour = "black", size = 10),
                                           axis.text.x = element_text(colour = "black",size = 10),
                                           axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
f1mat_precision <- f1mat_precision + geom_label_repel(force = 50, size=3)


## --- Cell Line --- ##

### 1. Age Analysis ###
combined_dataset_ERR6937788 <- read.delim("https://figshare.com/ndownloader/files/44176907", comment.char="#")
combined_dataset_ERR6937791 <- read.delim("https://figshare.com/ndownloader/files/44176916", comment.char="#")
combined_dataset_ERR6937794 <- read.delim("https://figshare.com/ndownloader/files/44176913", comment.char="#")
combined_dataset_ERR6937797 <- read.delim("https://figshare.com/ndownloader/files/44176919", comment.char="#")
combined_dataset_ERR6937800 <- read.delim("https://figshare.com/ndownloader/files/44176922", comment.char="#")
combined_dataset_HCT116 <- read.delim("https://figshare.com/ndownloader/files/44176910", comment.char="#")

age_process <- function(combined_age){
  combined_age <- combined_age[which((combined_age$Telescope != 0 | combined_age$LocusMasterTE != 0) & combined_age$long !=0),]
  combined_age$diff_raw <- as.numeric(combined_age$Telescope) - as.numeric(combined_age$LocusMasterTE)
  combined_age$corrected_by_modified_telescope <- "not captured"
  combined_age$corrected_by_modified_telescope[which(combined_age$diff_raw != 0)] <- "captured"
  combined_age$originally_correct <- "not captured"
  combined_age$originally_correct[which(combined_age$diff_raw == 0)] <- "captured"
  combined_age <- combined_age[,c(2,7,8)]
  colnames(combined_age)[2] <- "TEs Corrected by LocusMasterTE"
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
age_combined$age <- unlist(hg38_repeatmasker[match(age_combined$label, hg38_repeatmasker$label),3])

age_combined <- na.omit(age_combined)
age_combined[is.na(age_combined)] <- 0
stats <- compare_means(millidiv ~ variable, group.by = "group", data = age_combined, method = "t.test")
give.n <- age_combined %>% dplyr::group_by(group, variable) %>% dplyr::summarise(Freq=n())

age_combined_plot <- ggplot(age_combined, aes(x=variable, y=age, fill=variable)) + geom_boxplot(width=0.6, lwd=1)+facet_wrap(~group,nrow = 1)
age_combined_plot <- age_combined_plot + My_Theme + scale_fill_npg() + scale_y_continuous(limits=c(0,350))
age_combined_plot <- age_combined_plot + theme(axis.title.x=element_blank(),
                                               axis.text.x=element_blank(),
                                               axis.ticks.x=element_blank(),strip.background = element_blank(),
                                               strip.text.x = element_blank()) + ylab("Age (milliDiv)") +stat_compare_means(label =  "p.signif", label.x = 1.5)
age_combined_plot <- age_combined_plot + ggtitle("Age Distribution of TEs")
age_combined_plot

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


correlation <- function(combined_dataset,cell_line_name, short, long){
  combined_dataset_cor <- combined_dataset
  combined_dataset_cor$str <- unlist(sapply(strsplit(combined_dataset_cor$label,"|", fixed = TRUE), function(x) x[2], simplify=FALSE))
  combined_dataset_cor$end <- unlist(sapply(strsplit(combined_dataset_cor$label,"|", fixed = TRUE), function(x) x[3], simplify=FALSE))
  combined_dataset_cor$length <- as.numeric(combined_dataset_cor$end) - as.numeric(combined_dataset_cor$str)+1
  
  combined_dataset_cor$LocusMasterTE <- round(combined_dataset_cor$LocusMasterTE)
  
  ### change to TPM ### use sum of coding genes as library size
  long_read_fc_coding_gene_cleanup <- read.delim("https://figshare.com/ndownloader/files/44176547", comment.char="#")
  short_read_fc_coding_gene_cleanup <- read.delim("https://figshare.com/ndownloader/files/44176556", comment.char="#")
  short_read_fc_coding_gene_cleanup$rate <- as.numeric(short_read_fc_coding_gene_cleanup[,7]) / as.numeric(short_read_fc_coding_gene_cleanup$Length)
  long_read_fc_coding_gene_cleanup$rate <- as.numeric(long_read_fc_coding_gene_cleanup[,7]) / as.numeric(long_read_fc_coding_gene_cleanup$Length)
  
  library_size_short = sum(short_read_fc_coding_gene_cleanup$rate)
  library_size_long = sum(long_read_fc_coding_gene_cleanup$rate)
  
  tpm3 <- function(counts,len) {
    x <- as.numeric(counts)/as.numeric(len)
    return(t(t(x)*1e6))
  }
  
  combined_dataset_cor$longTPM <- (tpm3(as.numeric(combined_dataset_cor$long), as.numeric(combined_dataset_cor$length)))/library_size_long
  combined_dataset_cor$shortTPM <- (tpm3(as.numeric(combined_dataset_cor$Telescope), as.numeric(combined_dataset_cor$length)))/library_size_short
  combined_dataset_cor$LocusMasterTETPM <- (tpm3(as.numeric(combined_dataset_cor$LocusMasterTE), as.numeric(combined_dataset_cor$length)))/library_size_short
  
  combined_dataset_cor$diff <- as.numeric(combined_dataset_cor$shortTPM) - as.numeric(combined_dataset_cor$LocusMasterTETPM)
  combined_dataset_cor$diff_label <- "Same"
  combined_dataset_cor[which(combined_dataset_cor$diff != 0),"diff_label"] <- "Changed"
  combined_dataset_cor_common <- combined_dataset_cor[which(abs(combined_dataset_cor$diff)> 1e-5),]
  combined_dataset_cor_common <- combined_dataset_cor_common[,c(1,9:11)]
  combined_dataset_cor_common <- na.omit(combined_dataset_cor_common)
  rownames(combined_dataset_cor_common) <- combined_dataset_cor_common[,1]
  combined_dataset_cor_common <- combined_dataset_cor_common[rowSums(combined_dataset_cor_common[,2:4])>0,]
  combined_dataset_cor_common <- combined_dataset_cor_common[,2:4]
  
  combined_dataset_cor_common <- combined_dataset_cor_common[!(apply(combined_dataset_cor_common, 1, function(y) any(y == 0))),]
  for(i in 1:ncol(combined_dataset_cor_common)){
    combined_dataset_cor_common[,i] <- as.numeric(log(combined_dataset_cor_common[,i]+1e-5, base=10))
  }
  p1 <- ggscatter(combined_dataset_cor_common[,c(1,2)], y = "longTPM", x = "shortTPM",
                  add = "reg.line", color = "#7E6148B2",
                  add.params = list(color = "#3C5488B2", fill = "#8491B4B2", size=4),
                  conf.int = TRUE)+My_Theme
  p1 <- p1 + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size = 4)
  p1 <- p1 + ggtitle("Telescope") + xlab("log10 (TPM by Telescope)") + ylab("log10 (TPM by Long Read)")

  
  p2 <- ggscatter(combined_dataset_cor_common[,c(1,3)], y = "longTPM", x = "LocusMasterTETPM",
                  add = "reg.line", color = "#7E6148B2",
                  add.params = list(color = "#DC0000FF", fill = "#8491B4B2", size=4),
                  conf.int = TRUE)+My_Theme
  p2 <- p2 + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size = 4)
  p2 <- p2 + ggtitle("LocusMasterTE") + xlab("log10 (TPM by LocusMasterTE)") + ylab("log10 (TPM by Long Read)")

}

### 3. Correlation Bubble Plot ###

### heatmap for correlation
cor_df <- data.frame(rep(0, 2))
for(i in 1:5){
  cor_df <- cbind(cor_df, rep(0,2))
}
colnames(cor_df) <- c("inhouse","A549","HCT116","HepG2","K562","MCF-7")
rownames(cor_df) <- c("Telescope","LocusMasterTE")
cor_df[1,] <- c(0.5,0.51, 0.48,0.44,0.39,0.55)
cor_df[2,] <- c(0.54,0.55, 0.56,0.48,0.43,0.58)
cor_df <- cbind("label"=rownames(cor_df),cor_df)
cor_df$label <- factor(cor_df$label, levels=c("Telescope", "LocusMasterTE"))
cor_df <- cor_df[,c(1,3:7)]
cor_df <- melt(cor_df)
cor_df$log_value <- log10(cor_df$value)+10
cor_df_heatmap <- ggplot(cor_df, aes(x =label, y = variable, color=label,size=log_value)) +geom_point()+scale_size(range = c(9,15))
cor_df_heatmap <- cor_df_heatmap + My_Theme
cor_df_heatmap <- cor_df_heatmap + theme(strip.background = element_blank(), strip.text.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="bottom")
cor_df_heatmap <- cor_df_heatmap + scale_color_manual(values=c("#3C5488FF","#DC0000FF"))+geom_text(aes(label=value), size=4, color="white",fontface = "bold")
cor_df_heatmap <- cor_df_heatmap + ggtitle("Correlation with Long Read")+theme(legend.position="none")

### 4. Histone bed file ###
### bed file extraction for histone one example ###

histone_own <- combined_dataset_ERR6937800[,c(1,2,4,5)]
histone_own$diff <- histone_own$LocusMasterTE - histone_own$Telescope
histone_own$label_diff <- "same"
histone_own$label_diff[histone_own$diff > 1] <- "high_mixed"
histone_own$label_diff[histone_own$diff < -1] <- "high_short"
histone_own$chr <- unlist(sapply(strsplit(histone_own$label, "|", fixed=TRUE), function(x) x[1], simplify=FALSE))
histone_own$str <- unlist(sapply(strsplit(histone_own$label, "|", fixed=TRUE), function(x) x[2], simplify=FALSE))
histone_own$end <- unlist(sapply(strsplit(histone_own$label, "|", fixed=TRUE), function(x) x[3], simplify=FALSE))

write.table(histone_own[histone_own$label_diff == "high_mixed", 7:9], "~/bigger1_Corrected_in_modified_telescope.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
write.table(histone_own[histone_own$label_diff == "high_short", 7:9], "~/bigger1_Incorrectly_measured_in_original_telescope.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)


