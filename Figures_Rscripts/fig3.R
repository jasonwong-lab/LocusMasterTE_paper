My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

## Load dataset
combined <- read.delim("https://figshare.com/ndownloader/files/39628645", comment.char="#")

combined_extract <- combined[,c(1,4:11)]
f1mat <- data.frame(rep(0,3))
rownames(f1mat) <- c("Precision","Recall","F1")
for(i in c(1:6)){
  f1mat <- cbind(f1mat,rep(0,3))
}
colnames(f1mat) <- colnames(combined)[5:12]
## do it for 6 different quantification methods
for(idx in c(3:9)){
  color_temp <- "#91D1C2FF"
  scatterdf <- combined_extract[,c(1:2,idx)]
  colnames(scatterdf)[3] <- "Short"
  scatterdf_conf <- scatterdf
  
  ## confusion matrix
  scatterdf_conf$label <- "TP"
  scatterdf_conf$label[scatterdf_conf$Long != 0 & scatterdf_conf$Short ==0] <- "FN"
  scatterdf_conf$label[scatterdf_conf$Long == 0 & scatterdf_conf$Short !=0] <- "FP"
  scatterdf_conf$label[scatterdf_conf$Short == 0 & scatterdf_conf$Long ==0] <- "TN"
  
  confusion_mat <- scatterdf_conf %>% group_by(label) %>% summarise(Freq=n())
  confusion_mat
  if(idx == 4){
    precision <- as.numeric(confusion_mat[confusion_mat$label=="TP",2]) / (as.numeric(confusion_mat[confusion_mat$label=="TP",2])+0)
    
  }else{
    precision <- as.numeric(confusion_mat[confusion_mat$label=="TP",2]) / (as.numeric(confusion_mat[confusion_mat$label=="TP",2])+as.numeric(confusion_mat[confusion_mat$label=="FP",2]))
  }
  recall <- as.numeric(confusion_mat[confusion_mat$label=="TP",2]) / (as.numeric(confusion_mat[confusion_mat$label=="TP",2])+as.numeric(confusion_mat[confusion_mat$label=="FN",2]))
  f1 <- (2*precision*recall)/(precision+recall)
  
  f1mat[,idx-2] <- c(precision, recall, f1)
}

graph_list <- data.frame(c(3:9))
graph_list$color <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
graph_list$title <- c("Telescope","LocusMasterTE", "SQuIRE","Unique Counts","Best Counts","RSEM","SalmonTE")
graph_list$file_title <- c("Telescope","Our_Method", "SQuIRE","uniq","best","RSEM","salmonTE")

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
  scatterdf <- scatterdf[scatterdf$label !="TN",]
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
f1mat_plot$label <- c("Telescope","LocusMasterTE", "SQuIRE","Unique Counts","Best Counts","RSEM","SalmonTE")
f1mat_plot$label <- factor(f1mat_plot$label, levels=c("Unique Counts","Best Counts","RSEM","Telescope","LocusMasterTE", "SQuIRE","SalmonTE"))

f1mat_precision <- ggplot(aes(x = Recall, y = Precision,label=label), data = f1mat_plot) +geom_point(aes(fill=label), shape=21, size= 6, color="black")
f1mat_precision <- f1mat_precision +scale_fill_manual(values=c("#3C5488FF","#F39B7FFF","#8491B4FF","#E64B35FF","#4DBBD5FF","#00A087FF","#91D1C2FF"))
f1mat_precision <- f1mat_precision + ggtitle("Precision vs Recall")

f1mat_precision <- f1mat_precision + theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="blank",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                                           axis.title.x = element_text(colour = "black", size = 10),
                                           axis.text.x = element_text(colour = "black",size = 10),
                                           axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
f1mat_precision <- f1mat_precision + geom_label_repel(force = 50, size=3)
