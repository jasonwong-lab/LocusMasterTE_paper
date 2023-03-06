My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

### prepare_extraction
TCGA_COAD_TE_exp <- get(load("https://figshare.com/ndownloader/files/39384629"))
diff_TE <- read.delim("https://figshare.com/ndownloader/files/39384632")
diff_TE <- diff_TE[diff_TE$diff > 0.1,]
only_mod_imp <- diff_TE[diff_TE$loc != "exonic",]
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

### Look at PFI ###
## PFI ##
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
tot <- c("LSU-rRNA_Hsa_401","THE1C_758")
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
