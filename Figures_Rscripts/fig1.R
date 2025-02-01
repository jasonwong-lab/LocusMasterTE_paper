library(tidyverse)
library(ggpubr)
library(ggpointdensity)

My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

### ---- Figure1-B ---- ###

## Load the Datasets ##
short_read <- read.delim("https://figshare.com/ndownloader/files/52072676", comment.char="#")
long_read <- read.delim("https://figshare.com/ndownloader/files/52072673", comment.char="#")

## change to TPM counts
long_read_fc_coding_gene_cleanup <- read.delim("https://figshare.com/ndownloader/files/52072667", comment.char="#")
short_read_fc_coding_gene_cleanup <- read.delim("https://figshare.com/ndownloader/files/52072670", comment.char="#")
short_read_fc_coding_gene_cleanup$rate <- as.numeric(short_read_fc_coding_gene_cleanup[,7]) / as.numeric(short_read_fc_coding_gene_cleanup$Length)
long_read_fc_coding_gene_cleanup$rate <- as.numeric(long_read_fc_coding_gene_cleanup[,7]) / as.numeric(long_read_fc_coding_gene_cleanup$Length)

library_size_short = sum(short_read_fc_coding_gene_cleanup$rate)
library_size_long = sum(long_read_fc_coding_gene_cleanup$rate)

tpm3 <- function(counts,len) {
  x <- as.numeric(counts)/as.numeric(len)
  return(t(t(x)*1e6))
}

short_read$TPM <- (tpm3(as.numeric(short_read[,7]), as.numeric(short_read$Length)))/library_size_short
short_read <- short_read[!is.infinite(rowSums(short_read$TPM)),]

long_read$TPM <- (tpm3(as.numeric(long_read[,7]), as.numeric(long_read$Length)))/library_size_long
long_read <- long_read[!is.infinite(rowSums(long_read$TPM)),]

## merge two datasets for comparison ##
combined_short_long <- long_read[,c(1,8)]
colnames(combined_short_long) <- c("TE_name","long_read")
combined_short_long$short_read <- short_read[match(combined_short_long$TE_name,short_read$Geneid),8]
combined_short_long[is.na(combined_short_long)] <- 0

## excluding simple repeats
hg38_repeatmasker <- read.delim("https://figshare.com/ndownloader/files/39038723", header=FALSE, comment.char="#")
hg38_repeatmasker$V7 <- hg38_repeatmasker$V7+1
hg38_repeatmasker$label <- paste(hg38_repeatmasker$V6, paste(hg38_repeatmasker$V7, hg38_repeatmasker$V8, sep="|"), sep="|")

hg38_gencode_rmsk_indi <- read.delim("https://figshare.com/ndownloader/files/39038720", header=FALSE)
hg38_gencode_rmsk_indi$label <- paste(hg38_gencode_rmsk_indi$V1, paste(hg38_gencode_rmsk_indi$V4, hg38_gencode_rmsk_indi$V5, sep="|"), sep="|")
hg38_gencode_rmsk_indi$subF <- hg38_repeatmasker[match(hg38_gencode_rmsk_indi$label,hg38_repeatmasker$label),11]
hg38_gencode_rmsk_indi$indi_label <- unlist(sapply(strsplit(hg38_gencode_rmsk_indi$V9, " ", fixed=TRUE), function(x) x[2], simplify=FALSE))
hg38_gencode_rmsk_indi$indi_label_final <- unlist(sapply(strsplit(hg38_gencode_rmsk_indi$indi_label, ";", fixed=TRUE), function(x) x[1], simplify=FALSE))
hg38_gencode_rmsk_indi$Fam <- hg38_repeatmasker[match(hg38_gencode_rmsk_indi$label,hg38_repeatmasker$label),13]

combined_short_long$subF <- hg38_gencode_rmsk_indi[match(combined_short_long$TE_name, hg38_gencode_rmsk_indi$indi_label_final),11]
combined_short_long$Fam <- hg38_gencode_rmsk_indi[match(combined_short_long$TE_name, hg38_gencode_rmsk_indi$indi_label_final),14]

## summarise data for TE individual
combined_short_long <- combined_short_long[combined_short_long$Fam != "Simple_repeat",]

## summarise data for subF
combined_short_long_subF <- data.frame(cbind(combined_short_long %>% dplyr::group_by(subF) %>% dplyr::summarise("long_read_subF" = sum(long_read)), combined_short_long %>% dplyr::group_by(subF) %>% dplyr::summarise("short_read_subF" = sum(short_read))))
combined_short_long_subF <- combined_short_long_subF[,c(1,2,4)]

## data preparation for coding genes - counted by featureCounts
long_read_fc_coding_gene_cleanup$TPM <- (tpm3(as.numeric(long_read_fc_coding_gene_cleanup[,7]), as.numeric(long_read_fc_coding_gene_cleanup$Length)))/library_size_long
short_read_fc_coding_gene_cleanup$TPM <- (tpm3(as.numeric(short_read_fc_coding_gene_cleanup[,7]), as.numeric(short_read_fc_coding_gene_cleanup$Length)))/library_size_short


### 1. Venn Diagram ###
## individual 3.5 4.5 plot
long_short_gg <- list("Long Read" = c(combined_short_long$TE_name[combined_short_long$long_read != 0]), "Short Read" = c(combined_short_long$TE_name[combined_short_long$short_read != 0]))
p <- ggvenn(long_short_gg,stroke_size = 0.5, set_name_size=4, text_size = 4) +scale_fill_npg()

## subF
long_short_gg_subf <- list("Long Read" = c(combined_short_long_subF$subF[combined_short_long_subF$long_read_subF != 0]), "Short Read" = c(combined_short_long_subF$subF[combined_short_long_subF$short_read_subF != 0]))
p_subF <- ggvenn(long_short_gg_subf,stroke_size = 0.5, set_name_size=4, text_size = 4) +scale_fill_npg()

### 2. Correlation Plot ###

## indi --> 5.5 6.5 plot ## - log(only captured one)
combined_short_long <- na.omit(combined_short_long)
combined_short_long_cor <- combined_short_long[,2:3]
rownames(combined_short_long_cor) <- combined_short_long$TE_name
combined_short_long_cor <- combined_short_long_cor[!(apply(combined_short_long_cor, 1, function(y) any(y == 0))),]
combined_short_long_cor$long_read <- log(combined_short_long_cor$long_read + 1e-7, base=10)
combined_short_long_cor$short_read <- log(combined_short_long_cor$short_read + 1e-7, base=10)
combined_short_long_cor_indi_graph_common <- ggplot(combined_short_long_cor, aes(x=long_read,y=short_read)) + geom_pointdensity(color = "#7E6148B2")
combined_short_long_cor_indi_graph_common <- combined_short_long_cor_indi_graph_common + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size = 4)
combined_short_long_cor_indi_graph_common <- combined_short_long_cor_indi_graph_common + geom_smooth(method='lm', formula= y~x, color="#3C5488B2", size=4)
combined_short_long_cor_indi_graph_common <- combined_short_long_cor_indi_graph_common + xlab("Long Read at TE Individual (log10 TPM)") + ylab("Short Read at TE Individual (log10 TPM)")+My_Theme

## subF --> 5.5 6.5 plot ## - log(only captured one)
combined_short_long_subF <- na.omit(combined_short_long_subF)
combined_short_long_cor_subF <- combined_short_long_subF[,2:3]
rownames(combined_short_long_cor_subF) <- combined_short_long_subF[,1]
combined_short_long_cor_subF <- combined_short_long_cor_subF[!(apply(combined_short_long_cor_subF, 1, function(y) any(y == 0))),]
combined_short_long_cor_subF$long_read_subF <- log(combined_short_long_cor_subF$long_read_subF + 1e-7, base=10)
combined_short_long_cor_subF$short_read_subF <- log(combined_short_long_cor_subF$short_read_subF + 1e-7, base=10)
combined_short_long_cor_graph_common <- ggplot(combined_short_long_cor_subF, aes(x=long_read_subF,y=short_read_subF)) + geom_pointdensity(color = "#7E6148B2")
combined_short_long_cor_graph_common <- combined_short_long_cor_graph_common + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size = 4)
combined_short_long_cor_graph_common <- combined_short_long_cor_graph_common + geom_smooth(method='lm', formula= y~x, color="#3C5488B2", size=4)
combined_short_long_cor_graph_common <- combined_short_long_cor_graph_common + xlab("Long Read at Subfamily (log10 TPM)") + ylab("Short Read at Subfamily (log10 TPM)")+My_Theme

### ---- Figure1-E ---- ###
### Age Histogram ###
combined_short_long <- na.omit(combined_short_long)
combined_short_long$label <- hg38_gencode_rmsk_indi[match(combined_short_long$TE_name,hg38_gencode_rmsk_indi$indi_label_final),10]
combined_short_long$only_long <- 0
combined_short_long[which(combined_short_long$long_read != 0 & combined_short_long$short_read == 0),7] <- 1
combined_short_long$only_short <- 0
combined_short_long[which(combined_short_long$long_read == 0 & combined_short_long$short_read != 0),8] <- 1
combined_short_long_age <- combined_short_long[,c(6,2:3,7:8)]
colnames(combined_short_long_age)[2:5] <- c("Long","Short", "Only in\nLong","Only in\nShort")
combined_short_long_age <- melt(combined_short_long_age)
combined_short_long_age <- combined_short_long_age[combined_short_long_age$value !=0,]
combined_short_long_age$milidiv <- hg38_repeatmasker[match(combined_short_long_age$label,hg38_repeatmasker$label),3]
combined_short_long_age$milidiv <- combined_short_long_age$milidiv*0.1
my_comparisons <- list(c("Long", "Short"), c("Only in\nLong", "Only in\nShort"))
p3_box <- ggboxplot(combined_short_long_age, x="variable", y="milidiv", fill="variable", width=0.6, lwd=1) +scale_fill_npg()
p3_box <- p3_box + stat_compare_means(comparisons = my_comparisons,  label.y=c(31, 31), label = "p.signif", method="t.test")
p3_box <- p3_box+ xlab("")+ylab("%Divergence")+My_Theme+ theme(legend.position="none")
p3_box <- p3_box+ggtitle("Mean (Detected TEs divergence)")+scale_y_continuous(limits=c(15, 34), breaks=c(20, 25, 30), labels=c(20, 25, 30))

### ---- Figure1-A ---- ###
### Unique+Multi
map_1_short = 320320-105560
map_m_short = 105560
map_tot_short = 320320
map_1_short = map_1_short*100/map_tot_short
map_m_short = map_m_short*100/map_tot_short

map_1_long = 300656-43826
map_m_long = 43826
map_tot_long = 300656
map_1_long = map_1_long*100/map_tot_long
map_m_long = map_m_long*100/map_tot_long

p4_df <- c(map_1_short,map_m_short)
p4_df <- data.frame(p4_df)
p4_df$label <- c("Unique-mapped","Multi-mapped")
colnames(p4_df)[1] <- "Short"
p4_df <- p4_df[,c(2,1)]
p4_df$Long <- c(map_1_long,map_m_long)
p4_df <- melt(p4_df)
p4 <- ggplot(p4_df, aes(x = variable, y=value, fill=label, group=label)) + geom_col(width=0.7) + scale_fill_manual(values=c("#8491B4FF", "#91D1C2FF"))
p4 <- p4 + xlab("")+ylab("Percentage(%)") + My_Theme
p4 <- p4 + theme(legend.position="right")


### ---- Figure1-D and E ---- ###
### EM can not be the ultimate solution ###

random <- read.delim("https://figshare.com/ndownloader/files/52072664", header=FALSE, comment.char="#")
telescope_final <- read.delim("https://figshare.com/ndownloader/files/52072661")
locusmasterte_final <- read.delim("https://figshare.com/ndownloader/files/52072658", comment.char="#")

EM_not_ultimate <- hg38_gencode_rmsk_indi[,c(13,10)]
EM_not_ultimate$best <- unlist(random[match(EM_not_ultimate$indi_label_final,random$V1),9])
EM_not_ultimate$final_count <- unlist(telescope_final[match(EM_not_ultimate$indi_label_final,telescope_final$transcript),2])
EM_not_ultimate$locusmaster <- unlist(locusmasterte_final[match(EM_not_ultimate$indi_label_final,locusmasterte_final$transcript),2])

EM_not_ultimate <- EM_not_ultimate[,c(1,3,4,5)]
EM_not_ultimate[is.na(EM_not_ultimate)] <- 0
EM_not_ultimate <- EM_not_ultimate[rowSums(EM_not_ultimate[,2:4])>0,]
EM_not_ultimate$em_fail <- as.numeric(EM_not_ultimate$best)-as.numeric(EM_not_ultimate$final_count)
EM_not_ultimate$em_solve <- as.numeric(EM_not_ultimate$locusmaster)-as.numeric(EM_not_ultimate$final_count)
EM_not_ultimate <- EM_not_ultimate[order(EM_not_ultimate$indi_label_final),]
EM_not_ultimate <- EM_not_ultimate[21142:nrow(EM_not_ultimate),]

add_fig <- data.frame(c("Total","telescope_suc","Em_fail","unsolved"))
add_fig <- cbind(add_fig, "num"=c(nrow(EM_not_ultimate), nrow(EM_not_ultimate[EM_not_ultimate$em_fail==0 ,]),nrow(EM_not_ultimate[EM_not_ultimate$em_fail !=0,]),nrow(EM_not_ultimate[EM_not_ultimate$em_fail !=0 & EM_not_ultimate$em_solve ==0,])))
add_fig <- rbind(add_fig, c("locusmaster_suc",(add_fig[1,2]-add_fig[4,2])))
add_fig_pie <- add_fig[c(2:5),]
add_fig_pie$num <- as.numeric(add_fig_pie$num)

## EM effected ones age
EM_not_ultimate$label <- hg38_gencode_rmsk_indi[match(EM_not_ultimate$indi_label_final,hg38_gencode_rmsk_indi$indi_label_final),"label"]
EM_not_ultimate$age <- unlist(hg38_repeatmasker[match(EM_not_ultimate$label,hg38_repeatmasker$label),3])
EM_not_ultimate$age <- EM_not_ultimate$age *0.1
EM_not_ultimate$em <- "EM successfully\nconclude"
EM_not_ultimate$em[EM_not_ultimate$em_fail !=0] <- "EM fails\nto conclude"

em_age <- ggboxplot(EM_not_ultimate, x="em", y="age", fill="em", width=0.3, lwd=1) +scale_fill_manual(values=c("#661100", "#44AA99"))
em_age <- em_age + stat_compare_means( label = "p.signif", method="t.test")
em_age <- em_age+ xlab("")+ylab("%Divergence")+My_Theme+ theme(legend.position="none") 
em_age <- em_age+ggtitle("Mean (Detected TEs divergence)")

## bar chart
add_fig_pie$text <- add_fig_pie$num*100/as.numeric(add_fig[1,2])
add_fig_pie$text <- format(round(as.numeric(add_fig_pie$text),3), nsmall = 2)
add_fig_pie$final_text <- paste(add_fig_pie$num, paste0(round(as.numeric(add_fig_pie$text)),"%"), sep="\n")
add_fig_pie$group <- c("EM successfully\nconclude","EM fails\nto conclude","EM fails\nto conclude","EM successfully\nconclude")
add_fig_pie$tool <- c("Telescope","Telescope","LocusMasterTE","LocusMasterTE")
add_fig_pie$text <- as.numeric(add_fig_pie$text)
p_bar_chart <- ggplot(add_fig_pie, aes(x=tool, y=text, fill=group))+geom_col(width=0.4)
p_bar_chart <- p_bar_chart +My_Theme + scale_fill_manual(values=c("#661100", "#44AA99"))
p_bar_chart <- p_bar_chart+xlab("Captured TEs (N=625591)") + ylab("Percentage (%)")
p_bar_chart <- p_bar_chart +theme(legend.position="bottom") + coord_flip()
