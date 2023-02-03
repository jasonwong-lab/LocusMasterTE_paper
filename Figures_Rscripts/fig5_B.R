My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

combined_raw <- read.delim("https://figshare.com/ndownloader/files/39038678", comment.char="#")
hg38_gencode_rmsk_indi_loc_annotate <- read.delim("https://figshare.com/ndownloader/files/39052802", header=FALSE)

## From Spanki select only TP from lasTEq
tp_org_tel <- combined_raw[which(combined_raw$Long!=0 & combined_raw$org_tel!=0),1:2]
tp_mod_tel <- combined_raw[which(combined_raw$Long!=0 & combined_raw$mod_tel!=0),1:2]
only_mod <- tp_mod_tel[!tp_mod_tel$TE_name %in% tp_org_tel$TE_name,]

only_mod$loc <- unlist(hg38_gencode_rmsk_indi_loc_annotate[match(only_mod$TE_name,hg38_gencode_rmsk_indi_loc_annotate$V4),5])

### Extract tumor sample only
TCGA_COAD_TE_exp <- get(load("https://figshare.com/ndownloader/files/39052376"))
only_mod_imp <- only_mod[only_mod$loc != "exonic",]
TCGA_COAD_TE_exp_nec_no0 <- TCGA_COAD_TE_exp[rownames(TCGA_COAD_TE_exp) %in% only_mod_imp$TE_name,]
TCGA_COAD_TE_exp_nec_no0 <- TCGA_COAD_TE_exp_nec_no0[,str_detect(colnames(TCGA_COAD_TE_exp_nec_no0),"01A")]
colnames(TCGA_COAD_TE_exp_nec_no0) <- substr(colnames(TCGA_COAD_TE_exp_nec_no0),1,12)

### run Lasso-cox regression model ###
library(glmnet)
library(survival)
survival_four_raw <- read.delim("https://figshare.com/ndownloader/files/39038849")
survival_four_raw <- survival_four_raw[str_detect(survival_four_raw$sample,"-01"),]
survival_four_raw$X_PATIENT <- gsub("-",".", survival_four_raw$X_PATIENT, fixed=TRUE)
survival_four_raw <- survival_four_raw[survival_four_raw$OS.time != 0,]
TCGA_COAD_TE_exp_nec_no0 <- TCGA_COAD_TE_exp_nec_no0[,colnames(TCGA_COAD_TE_exp_nec_no0) %in% survival_four_raw$X_PATIENT]
TCGA_COAD_TE_exp_nec_no0 <- data.frame(TCGA_COAD_TE_exp_nec_no0)

## only seelect TE with 80% expressed
threshold <- ncol(TCGA_COAD_TE_exp_nec_no0) * 0.2
TCGA_COAD_TE_exp_nec_no0 <- TCGA_COAD_TE_exp_nec_no0[rowSums(TCGA_COAD_TE_exp_nec_no0 == 0) <= threshold, ]
survival_four_raw <- survival_four_raw[survival_four_raw$X_PATIENT %in% colnames(TCGA_COAD_TE_exp_nec_no0),]
for(i in 1:nrow(TCGA_COAD_TE_exp_nec_no0)){
  survival_four_raw <- cbind(survival_four_raw, name = unlist(TCGA_COAD_TE_exp_nec_no0[i,match(survival_four_raw$X_PATIENT,colnames(TCGA_COAD_TE_exp_nec_no0))]))
}
colnames(survival_four_raw)[12:ncol(survival_four_raw)] <- rownames(TCGA_COAD_TE_exp_nec_no0)

te <- survival_four_raw[,12:ncol(survival_four_raw)]
te <- as.matrix(sapply(te, as.numeric))  

### Look at OS and PFI ###
## OS ##
for(i in c(1:100)){
  cv.fit2 = cv.glmnet(te,
                      Surv(survival_four_raw$OS.time, survival_four_raw$OS),
                      alpha = 1,
                      family = "cox")
  small.lambda.index <- which(cv.fit2$lambda == cv.fit2$lambda.min)
  small.lambda.betas <- cv.fit2$glmnet.fit$beta[,small.lambda.index]
  small.lambda.betas <- data.frame(small.lambda.betas)
  small.lambda.betas <- rownames(small.lambda.betas)[small.lambda.betas$small.lambda.betas !=0]
}

### Selected 5 TEs
tot <- c("AluSp_47159","MIRb_187527","PABL_A-int_27","AluSp_9954","AluYe5_1025")

### find optimal clusters ###
library(ConsensusClusterPlus)
TCGA_COAD_TE_extract <- TCGA_COAD_TE_exp_nec_no0[rownames(TCGA_COAD_TE_exp_nec_no0) %in% tot,]
TCGA_COAD_TE_extract = sweep(TCGA_COAD_TE_extract,1, apply(TCGA_COAD_TE_extract,1,mean,na.rm=T))
TCGA_COAD_TE_extract <- data.matrix(TCGA_COAD_TE_extract)
results = ConsensusClusterPlus(TCGA_COAD_TE_extract,maxK=10,reps=100,pFeature=1,distance="euclidean",clusterAlg="pam",plot="pdf")

temp <- data.frame(results[[2]]$consensusClass)
survival_four <- survival_four_raw[,c(1:10)]
survival_four$TE_group <- unlist(temp[match(survival_four$X_PATIENT,rownames(temp)),1])

## total 150 TEs --> FDR < 0.05, P-value < 0.05
## OS ##
fit <- survfit(Surv(OS.time, OS) ~ TE_group, data = survival_four)
library(survminer)
ggsurvplot(fit,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata", 
           linetype = "strata",
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = My_Theme,xlab = "OS(time)")

### functionally + location annotate TEs ###
tot <- data.frame(tot)
tot$loc <- unlist(hg38_gencode_rmsk_indi_loc_annotate[match(tot$tot,hg38_gencode_rmsk_indi_loc_annotate$V4),5])
tot$chr <- unlist(hg38_gencode_rmsk_indi_loc_annotate[match(tot$tot,hg38_gencode_rmsk_indi_loc_annotate$V4),1])
tot$start <- unlist(hg38_gencode_rmsk_indi_loc_annotate[match(tot$tot,hg38_gencode_rmsk_indi_loc_annotate$V4),2])
tot$end <- unlist(hg38_gencode_rmsk_indi_loc_annotate[match(tot$tot,hg38_gencode_rmsk_indi_loc_annotate$V4),3])