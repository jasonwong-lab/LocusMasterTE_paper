library(reshape)
My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

## ------1. GSVA Graph------ ##

### compare correlation between DNA RNA sensing and TE
DNA_RNA_gene_TE_cor_GSVA_mod <- read.delim("https://figshare.com/ndownloader/files/44481905")
DNA_RNA_gene_TE_cor_GSVA_org <- read.delim("https://figshare.com/ndownloader/files/44481881")
DNA_RNA_gene_TE_cor_GSVA_org$label <- "Telescope"
DNA_RNA_gene_TE_cor_GSVA_mod$label <- "LocusMasterTE"

DNA_RNA_gene_TE_cor_GSVA_org$mod <- unlist(DNA_RNA_gene_TE_cor_GSVA_mod[match(DNA_RNA_gene_TE_cor_GSVA_org$name, DNA_RNA_gene_TE_cor_GSVA_mod$name),3])
p_cor <- ggplot(DNA_RNA_gene_TE_cor_GSVA_org, aes(x = abs(R_value), y=abs(mod))) + geom_point(shape=1, size=2, color="#00A087FF")+geom_abline(slope=1, linewidth =1,linetype = 'dashed', color = "#7E6148FF")
p_cor <- p_cor + My_Theme + xlab("abs(cor.value) by Telescope")+ylab("abs(cor.value) by LocusMasterTE") + ggtitle("Cor.value of TE with GSVA score \n for DNA/RNA sensing gene")

## ------2. Correlation Test------ ##

## Intergenic TEs with 5kb regions
intergenic_TE_gene <- read.delim("https://figshare.com/ndownloader/files/44481866")

intergenic_TE_gene$diff <- abs(as.numeric(intergenic_TE_gene$LocusMasterTE_cor)) - abs(as.numeric(intergenic_TE_gene$Telescope_cor))

intergenic_TE_gene$group <- "abs(<0.5)"
intergenic_TE_gene$group[abs(intergenic_TE_gene$LocusMasterTE_cor) >= 0.5 | abs(intergenic_TE_gene$Telescope_cor) >= 0.5] <- "abs(>=0.5) in either"
intergenic_TE_gene$high <- "Same"
intergenic_TE_gene$high[intergenic_TE_gene$diff > 0.01] <- "high in Our"
intergenic_TE_gene$high[intergenic_TE_gene$diff < -0.01] <- "high in Telescope"

intergenic_TE_gene$color_label <- "1Others"
intergenic_TE_gene$color_label[intergenic_TE_gene$group == "abs(>=0.5) in either" & intergenic_TE_gene$high=="high in Our"] <- "3abs(>=0.5) & \nhigh in LocusMasterTE"
intergenic_TE_gene$color_label[intergenic_TE_gene$group == "abs(>=0.5) in either" & intergenic_TE_gene$high=="high in Telescope"] <- "2abs(>=0.5) & \nhigh in Telescope"

p_cor_5kb_all <- ggplot(intergenic_TE_gene, aes(x = abs(Telescope_cor), y=abs(LocusMasterTE_cor),color=color_label)) +  geom_pointdensity(alpha=0.8) +geom_abline(slope=1, linewidth =1,linetype = 'dashed', color = "#7E6148FF")
p_cor_5kb_all <- p_cor_5kb_all + My_Theme + xlab("abs(cor.value) by Telescope")+ylab("abs(cor.value) by LocusMasterTE") + ggtitle("Cor.value of intergenic TE +/- 5kb \nwith coding genes") + theme(legend.position = "none")
p_cor_5kb_all <- p_cor_5kb_all + scale_colour_manual(values = c("grey", "#3C5488FF", "#DC0000FF"))+scale_fill_discrete(labels=c('Others', 'abs(>=0.5) & \nhigh in Telescope', 'abs(>=0.5) & \nhigh in LocusMasterTE'))

intergenic_TE_gene_all <- intergenic_TE_gene %>% dplyr::group_by(group, high) %>% dplyr::summarise(Freq=n())
intergenic_TE_gene_all <- intergenic_TE_gene_all[intergenic_TE_gene_all$high != "Same",]
intergenic_TE_gene_all <- cast(intergenic_TE_gene_all, group~high) 

chisq.test(intergenic_TE_gene_all)
intergenic_TE_gene_all


## ------3. Survival------ ##
TCGA_COAD_TE_exp <- read.delim("https://figshare.com/ndownloader/files/44569130")
locusmasterte_te <- read.table("https://figshare.com/ndownloader/files/44481584", quote="\"", comment.char="")
diff_TE <- read.delim("https://figshare.com/ndownloader/files/44481530", header=FALSE)
only_mod_imp <- diff_TE[diff_TE$V5 == "intergenic",]

TCGA_COAD_TE_exp_nec_no0 <- TCGA_COAD_TE_exp[rownames(TCGA_COAD_TE_exp) %in% only_mod_imp$V4,]
TCGA_COAD_TE_exp_nec_no0 <- TCGA_COAD_TE_exp_nec_no0[,str_detect(colnames(TCGA_COAD_TE_exp_nec_no0),"01A")]
colnames(TCGA_COAD_TE_exp_nec_no0) <- substr(colnames(TCGA_COAD_TE_exp_nec_no0),1,12)

### run Lasso-cox regression model ###
library(glmnet)
library(survival)
survival_four_raw <- read.delim("https://figshare.com/ndownloader/files/44483729")
survival_four_raw <- survival_four_raw[str_detect(survival_four_raw$sample,"-01"),]
survival_four_raw$X_PATIENT <- gsub("-",".", survival_four_raw$X_PATIENT, fixed=TRUE)
survival_four_raw <- survival_four_raw[survival_four_raw$OS.time != 0,]
TCGA_COAD_TE_exp_nec_no0 <- TCGA_COAD_TE_exp_nec_no0[,colnames(TCGA_COAD_TE_exp_nec_no0) %in% survival_four_raw$X_PATIENT]
TCGA_COAD_TE_exp_nec_no0 <- data.frame(TCGA_COAD_TE_exp_nec_no0)

## only select TE with 80% expressed
threshold <- ncol(TCGA_COAD_TE_exp_nec_no0) * 0.8
TCGA_COAD_TE_exp_nec_no0 <- TCGA_COAD_TE_exp_nec_no0[rowSums(TCGA_COAD_TE_exp_nec_no0 == 0) <= threshold, ]

survival_four_raw <- survival_four_raw[survival_four_raw$X_PATIENT %in% colnames(TCGA_COAD_TE_exp_nec_no0),]
for(i in 1:nrow(TCGA_COAD_TE_exp_nec_no0)){
  survival_four_raw <- cbind(survival_four_raw, name = unlist(TCGA_COAD_TE_exp_nec_no0[i,match(survival_four_raw$X_PATIENT,colnames(TCGA_COAD_TE_exp_nec_no0))]))
}
colnames(survival_four_raw)[12:ncol(survival_four_raw)] <- rownames(TCGA_COAD_TE_exp_nec_no0)

survival_four_raw <- survival_four_raw[survival_four_raw$PFI.time > 0,]

te <- survival_four_raw[,12:ncol(survival_four_raw)]
te <- as.matrix(sapply(te, as.numeric))  

## save memory
TCGA_COAD_TE_exp <- NULL

### Look at PFI ###
## PFI ##
# "AluSx_44769"  "AluSp_29088"  "AluSg_30598"  "AluY_88947"   "AluY_91831"   "L1MC4a_10165"
for(i in c(1:100)){
  cv.fit2 = cv.glmnet(te,
                      Surv(survival_four_raw$PFI.time, survival_four_raw$PFI),
                      alpha = 1,
                      family = "cox")
  small.lambda.index <- which(cv.fit2$lambda == cv.fit2$lambda.min)
  small.lambda.betas <- cv.fit2$glmnet.fit$beta[,small.lambda.index]
  small.lambda.betas <- data.frame(small.lambda.betas)
  small.lambda.betas <- rownames(small.lambda.betas)[small.lambda.betas$small.lambda.betas !=0]
  print(small.lambda.betas)
}


### Selected TE
tot <- c("AluSx_44769","AluSp_29088","AluSg_30598","AluY_88947","AluY_91831","L1MC4a_10165")
### find optimal clusters ###
library(ConsensusClusterPlus)
TCGA_COAD_TE_extract <- TCGA_COAD_TE_exp_nec_no0[rownames(TCGA_COAD_TE_exp_nec_no0) %in% tot,]
TCGA_COAD_TE_extract = sweep(TCGA_COAD_TE_extract,1, apply(TCGA_COAD_TE_extract,1,mean,na.rm=T))
TCGA_COAD_TE_extract <- data.matrix(TCGA_COAD_TE_extract)

results = ConsensusClusterPlus(TCGA_COAD_TE_extract,maxK=10,reps=100,pFeature=1,distance="euclidean",clusterAlg="pam",plot="pdf")

temp <- data.frame(results[[2]]$consensusClass)
survival_four <- survival_four_raw[,c(1:10)]
survival_four$TE_group <- unlist(temp[match(survival_four$X_PATIENT,rownames(temp)),1])

## PFI ##
fit <- survfit(Surv(PFI.time, PFI) ~ TE_group, data = survival_four)

library(survminer)
ggsurvplot(fit,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           palette = c("#DC0000FF", "#3C5488FF"),
           ggtheme = My_Theme,xlab = "PFI(time)")

TCGA_COAD_TE_extract <- TCGA_COAD_TE_exp_nec_no0[rownames(TCGA_COAD_TE_exp_nec_no0) %in% tot,]
box_plot <- melt(cbind("name"=rownames(TCGA_COAD_TE_extract),TCGA_COAD_TE_extract))

box_plot$strata <- paste0("strata",unlist(temp[match(box_plot$variable,rownames(temp)),1]))
box_plot <- na.omit(box_plot)

box_plot$strata[box_plot$strata == "strata1"] <- "Poor"
box_plot$strata[box_plot$strata == "strata2"] <- "Good"

box_plot <- box_plot[box_plot$value>0,]
surv_box <- ggplot(box_plot, aes(x=strata, y=value, fill=strata)) + geom_quasirandom(alpha=0.5)+geom_boxplot(width=0.2, lwd=1) +scale_fill_manual(values=c("#3C5488FF", "#DC0000FF"))
surv_box <- surv_box + facet_wrap(~name,scales="free", ncol=3)
surv_box <- surv_box+ xlab("Survival Strata")+ylab("TE expression")+My_Theme

## 4------. Mutational Signature------ ##
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
