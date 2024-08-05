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

