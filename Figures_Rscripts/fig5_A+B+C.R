library(reshape)
My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

## GSVA Graph ##

### compare correlation between DNA RNA sensing and TE
DNA_RNA_gene_TE_cor_GSVA_mod <- read.delim("https://figshare.com/ndownloader/files/39370034")
DNA_RNA_gene_TE_cor_GSVA_org <- read.delim("https://figshare.com/ndownloader/files/39038912")
DNA_RNA_gene_TE_cor_GSVA_org$label <- "Telescope"
DNA_RNA_gene_TE_cor_GSVA_mod$label <- "lasTEq"

DNA_RNA_gene_TE_cor_GSVA_org$mod <- unlist(DNA_RNA_gene_TE_cor_GSVA_mod[match(DNA_RNA_gene_TE_cor_GSVA_org$name, DNA_RNA_gene_TE_cor_GSVA_mod$name),3])
p_cor <- ggplot(DNA_RNA_gene_TE_cor_GSVA_org, aes(x = abs(R_value), y=abs(mod))) + geom_point(shape=1, size=2, color="#00A087FF")+geom_abline(slope=1, linewidth =1,linetype = 'dashed', color = "#7E6148FF")
p_cor <- p_cor + My_Theme + xlab("abs(cor.value) by Telescope")+ylab("abs(cor.value) by lasTEq") + ggtitle("Cor.value of TE with GSVA score \n for DNA/RNA sensing gene")


## Intergenic TEs with 5kb regions
intergenic_TE_gene_5kb <- read.delim("https://figshare.com/ndownloader/files/39370052", comment.char="#")
UTR_3 <- read.delim("https://figshare.com/ndownloader/files/39044861", header=FALSE)
UTR_3$gene <- unlist(sapply(strsplit(UTR_3$V4, "_", fixed=TRUE), function(x) x[1], simplify=FALSE))
UTR_5 <- read.delim("https://figshare.com/ndownloader/files/39044867", header=FALSE)
UTR_5$gene <- unlist(sapply(strsplit(UTR_5$V4, "_", fixed=TRUE), function(x) x[1], simplify=FALSE))

intergenic_TE_gene_5kb <- intergenic_TE_gene_5kb[,c(1:4,11:16)]
intergenic_TE_gene_5kb <- unique(intergenic_TE_gene_5kb)

UTR_5$ENST <- substr(UTR_5$gene, 1, 15)
UTR_3$ENST <- substr(UTR_3$gene, 1, 15)

intergenic_TE_gene_5kb$UTR5 <- unlist(UTR_5[match(intergenic_TE_gene_5kb$V12,UTR_5$ENST),2])
intergenic_TE_gene_5kb$UTR3 <- unlist(UTR_3[match(intergenic_TE_gene_5kb$V12,UTR_3$ENST),2])

intergenic_TE_gene_5kb <- na.omit(intergenic_TE_gene_5kb)
intergenic_TE_gene_5kb$diff3 <- abs(as.numeric(intergenic_TE_gene_5kb$V2)-as.numeric(intergenic_TE_gene_5kb$UTR3))
intergenic_TE_gene_5kb$diff5 <- abs(as.numeric(intergenic_TE_gene_5kb$V2)-as.numeric(intergenic_TE_gene_5kb$UTR5))
intergenic_TE_gene_5kb$pos <- "upstream"
intergenic_TE_gene_5kb$pos[intergenic_TE_gene_5kb$diff3 < intergenic_TE_gene_5kb$diff5] <- "downstream"

intergenic_TE_gene_5kb$diff <- abs(as.numeric(intergenic_TE_gene_5kb$V15)) - abs(as.numeric(intergenic_TE_gene_5kb$V16))

intergenic_TE_gene_5kb <- intergenic_TE_gene_5kb[order(intergenic_TE_gene_5kb$diff, decreasing=TRUE),]
intergenic_TE_gene_5kb$order_x <- c(1:nrow(intergenic_TE_gene_5kb))

intergenic_TE_gene_5kb$group <- "abs(<0.5)"
intergenic_TE_gene_5kb$group[abs(intergenic_TE_gene_5kb$V15) >= 0.5 | abs(intergenic_TE_gene_5kb$V16) >= 0.5] <- "abs(>=0.5) in either"
intergenic_TE_gene_5kb$high <- "Same"
intergenic_TE_gene_5kb$high[intergenic_TE_gene_5kb$diff > 0.01] <- "high in Our"
intergenic_TE_gene_5kb$high[intergenic_TE_gene_5kb$diff < -0.01] <- "high in Telescope"

intergenic_TE_gene_5kb$color_label <- "1Others"
intergenic_TE_gene_5kb$color_label[intergenic_TE_gene_5kb$group == "abs(>=0.5) in either" & intergenic_TE_gene_5kb$high=="high in Our"] <- "3abs(>=0.5) & \nhigh in lasTEq"
intergenic_TE_gene_5kb$color_label[intergenic_TE_gene_5kb$group == "abs(>=0.5) in either" & intergenic_TE_gene_5kb$high=="high in Telescope"] <- "2abs(>=0.5) & \nhigh in Telescope"

p_cor_5kb <- ggplot(intergenic_TE_gene_5kb[intergenic_TE_gene_5kb$pos=="upstream",], aes(x = abs(V16), y=abs(V15),color=color_label)) +  geom_pointdensity(alpha=0.8) +geom_abline(slope=1, linewidth =1,linetype = 'dashed', color = "#7E6148FF")
p_cor_5kb <- p_cor_5kb + My_Theme + xlab("abs(cor.value) by Telescope")+ylab("abs(cor.value) by lasTEq") + ggtitle("Cor.value of intergenic TE +/- 5kb \nwith upstream of coding genes") + theme(legend.position = "none")
p_cor_5kb <- p_cor_5kb + scale_colour_manual(values = c("grey", "#3C5488FF", "#DC0000FF"))+scale_fill_discrete(labels=c('Others', 'abs(>=0.5) & \nhigh in Telescope', 'abs(>=0.5) & \nhigh in lasTEq'))

intergenic_TE_gene_5kb_df_up <- intergenic_TE_gene_5kb[intergenic_TE_gene_5kb$pos =="upstream",] %>% group_by(group, high) %>% summarise(Freq=n())
intergenic_TE_gene_5kb_df_up <- intergenic_TE_gene_5kb_df_up[intergenic_TE_gene_5kb_df_up$high != "Same",]
intergenic_TE_gene_5kb_df_up <- cast(intergenic_TE_gene_5kb_df_up, group~high) 

chisq.test(intergenic_TE_gene_5kb_df_up)
intergenic_TE_gene_5kb_df_up

p_cor_5kb_down <- ggplot(intergenic_TE_gene_5kb[intergenic_TE_gene_5kb$pos=="downstream",], aes(x = abs(V16), y=abs(V15),color=color_label)) +  geom_pointdensity(alpha=0.8) +geom_abline(slope=1, linewidth =1,linetype = 'dashed', color = "#7E6148FF")
p_cor_5kb_down <- p_cor_5kb_down + My_Theme + xlab("abs(cor.value) by Telescope")+ylab("abs(cor.value) by lasTEq") + ggtitle("Cor.value of intergenic TE +/- 5kb \nwith downstream of coding genes") + theme(legend.position = "none")
p_cor_5kb_down <- p_cor_5kb_down + scale_colour_manual(values = c("grey", "#3C5488FF", "#DC0000FF"))+scale_fill_discrete(labels=c('Others', 'abs(>=0.5) & \nhigh in Telescope', 'abs(>=0.5) & \nhigh in lasTEq'))

intergenic_TE_gene_5kb_df_down <- intergenic_TE_gene_5kb[intergenic_TE_gene_5kb$pos =="downstream",] %>% group_by(group, high) %>% summarise(Freq=n())
intergenic_TE_gene_5kb_df_down <- intergenic_TE_gene_5kb_df_down[intergenic_TE_gene_5kb_df_down$high != "Same",]
intergenic_TE_gene_5kb_df_down <- cast(intergenic_TE_gene_5kb_df_down, group~high) 

chisq.test(intergenic_TE_gene_5kb_df_down)
