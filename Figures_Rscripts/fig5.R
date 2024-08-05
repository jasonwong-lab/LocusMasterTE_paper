library(reshape)
My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(gridExtra)
library(lsa)

mut_mat <- read.delim("https://figshare.com/ndownloader/files/44481518")
mut_mat$average_common <- rowMeans(mut_mat[,str_detect(colnames(mut_mat), "common")])
mut_mat$average_lasTEq <- rowMeans(mut_mat[,str_detect(colnames(mut_mat), "mod")])
mut_mat$average_coding_gene <- rowMeans(mut_mat[,str_detect(colnames(mut_mat), "coding_gene")])

colnames(mut_mat)[c(1,(ncol(mut_mat)-2):ncol(mut_mat))] <- c("Reference","common\nTEs", "LocusMaster\nTE only", "Coding\nGenes")
plot_96_profile(mut_mat[,c(1,2,3)])/plot_96_profile(mut_mat[,c((ncol(mut_mat)-2):ncol(mut_mat))], ymax = 0.1)

mut_mat <- mut_mat[,1:(ncol(mut_mat)-3)]

mut_mat <- rbind("cosine" = rep(0),mut_mat)
for(i in 1:ncol(mut_mat)){
  mut_mat[1,i] <- cosine(mut_mat$Reference, mut_mat[,i])
}

## box-plot
mut_mat_graph <- mut_mat[1,2:ncol(mut_mat)]
mut_mat_graph <- melt(cbind(rownames(mut_mat_graph),mut_mat_graph))
mut_mat_graph$label <- "coding_genes"
mut_mat_graph$label[str_detect(mut_mat_graph$variable,"common")] <- "common"
mut_mat_graph$label[str_detect(mut_mat_graph$variable,"mod")] <- "LocusMasterTE only"

p3_box <- ggplot(mut_mat_graph, aes(x=label, y=value, fill=label)) + geom_quasirandom(alpha=0.7,color = "#B09C85FF") +scale_fill_npg() + geom_boxplot(width=0.2, lwd=1)
p3_box <- p3_box + xlab("") + ylab("Cosine")+My_Theme
