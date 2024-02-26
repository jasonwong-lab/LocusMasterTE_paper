library(ggsci)
My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

## L1PA13 TCGA-COAD strata ##
L1PA13_strata <- read.delim("https://figshare.com/ndownloader/files/44622562")
L1PA13_strata_mat <- data.frame(melt(cbind(name=rownames(L1PA13_strata), L1PA13_strata[,1:4])))
L1PA13_strata_mat$strata <- unlist(L1PA13_strata[match(L1PA13_strata_mat$name, rownames(L1PA13_strata)),5])
l1pa13 <- ggplot(L1PA13_strata_mat, aes(x=strata, y=value, fill=strata)) + geom_quasirandom()+ylab("log2(TPM)") + xlab("")+ ggtitle("L1PA13 Expression")+ geom_boxplot(width=0.2) + facet_wrap(~variable,scales="free_y")+My_Theme+scale_fill_jco()+ stat_compare_means(label =  "p.signif", label.x = 1.5)

survival_four_raw <- read.delim("https://figshare.com/ndownloader/files/39038849")
survival_four_raw$sample <- gsub("-",".",survival_four_raw$sample, fixed=TRUE)
survival_four_raw_nec <- survival_four_raw[survival_four_raw$sample %in% substr(L1PA13_strata_mat$name, 1, 15),]
survival_four_raw_nec$group <- unlist(L1PA13_strata_mat[match(survival_four_raw_nec$sample, substr(L1PA13_strata_mat$name, 1, 15)),4])

## Survival analysis ##
library(survminer)
library(survival)
fit <- survfit(Surv(DFI.time, DFI) ~ group, data = survival_four_raw_nec)

ggsurvplot(fit,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           palette = c("#DC0000FF", "#3C5488FF"),
           ggtheme = My_Theme,xlab = "DFI(time)")


## correlation test between L1PA13 and genes #3
cor_test_L1PA13 <- read.delim("https://figshare.com/ndownloader/files/44569085")
cor_test_L1PA13$name <- data.frame(mapIds(org.Hs.eg.db, keys = cor_test_L1PA13$gene, keytype="ENSEMBL", column = "SYMBOL"))
cor_test_L1PA13$name <- unlist(cor_test_L1PA13$name)
cor_test_L1PA13 <- cor_test_L1PA13[order(abs(cor_test_L1PA13$val), decreasing=TRUE),]
cor_test_L1PA13$name[1:2] <- c("FAR2P1","KLF2P1")
cor_test_L1PA13$gene <- factor(cor_test_L1PA13$gene, levels = cor_test_L1PA13$gene)
cor_test_L1PA13 <- na.omit(cor_test_L1PA13)
cor_test_L1PA13 <- data.frame(cor_test_L1PA13)
cor_test_L1PA13 <- cor_test_L1PA13[order(abs(cor_test_L1PA13$val), decreasing=TRUE),]
cor_test_L1PA13_nec <- cor_test_L1PA13[cor_test_L1PA13$val>0.5,]
cor_test_L1PA13_nec$name <- factor(cor_test_L1PA13_nec$name, levels = cor_test_L1PA13_nec$name)

cor_test <- ggplot(cor_test_L1PA13_nec, aes(x=name, y=val)) + geom_bar(stat="identity", fill="#7E6148FF", width=0.5) + My_Theme+ggtitle("Cor.var with L1PA13") + xlab("Gene") + ylab("cor.value")+ theme(axis.text.x = element_text(angle =45, hjust = 1))

## FAR2P1 expression graph
GENE_log2TPM <- read.delim("https://figshare.com/ndownloader/files/44569130")
GENE_log2TPM_far2p1 <- GENE_log2TPM[rownames(GENE_log2TPM)=="ENSG00000180178",]
GENE_log2TPM_far2p1 <- data.frame(melt(cbind("far2p1", GENE_log2TPM_far2p1)))
GENE_log2TPM_far2p1$strata <- unlist(potef_locusmaster_mat[match(GENE_log2TPM_far2p1$variable,potef_locusmaster_mat$Var2),4])
GENE_log2TPM_far2p1 <- na.omit(GENE_log2TPM_far2p1)
far2p1 <- ggplot(GENE_log2TPM_far2p1, aes(x=strata, y=value, fill=strata)) + geom_quasirandom()+ylab("log2 (TPM)") + xlab("") + ggtitle("FAR2P1 Expression")+ geom_boxplot(width=0.1) +My_Theme+scale_fill_jco()+ stat_compare_means(label =  "p.signif", label.x = 1.5)


## Check other factors ##
stage <- ggplot(L1PA13_strata, aes(x=strata, fill=stage)) +  geom_bar(position = "fill", width=0.7) +scale_x_discrete(guide = guide_axis(angle = 45))+xlab("")+ylab("Proportion (%)")+scale_y_continuous(labels = scales::percent) + My_Theme + scale_fill_nejm()+theme(legend.position = "right")
L1PA13_strata$MSI[L1PA13_strata$MSI=="MSS" | L1PA13_strata$MSI=="MSI-L"] <- "MSS"
L1PA13_strata$MSI[L1PA13_strata$MSI=="MSI-H"] <- "MSI"
L1PA13_strata$MSI[L1PA13_strata$MSI=="none" | L1PA13_strata$MSI=="Indeterminate"] <- "Indeterminate"
msi <- ggplot(L1PA13_strata, aes(x=strata, fill=MSI)) +  geom_bar(position = "fill", width=0.7) +scale_x_discrete(guide = guide_axis(angle = 45))+xlab("")+ylab("Proportion (%)")+scale_y_continuous(labels = scales::percent) + My_Theme + scale_fill_ordinal()+theme(legend.position = "right")
