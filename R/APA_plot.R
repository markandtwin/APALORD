#' Visualize the 3'end APA changes between the two samples
#' 
#' This function loads a gtf file and extract necessary information of each gene from it.
#' @param APA_data APA data generated by APA_profile function
#' @param P_cutoff adjusted P value cutoff to be considered as significant change
#' @param delta distance between two samples to be considered as significant APA events
#' @param internal_priming_exclude logic value for whether or not to exclude genes with detected internal_priming events
#' @param split.by.APA_type whether or not to split the plot by APA_type.
#' @return a data frame showing the significant APA change between sample1 and sample2
#' @export

APA_plot <- function(APA_data,P_cutoff=0.05,delta=0.1,internal_priming_exclude=F,split.by.APA_type=F){
  if(internal_priming_exclude){
    APA_table <-APA_data[is.na(internal_priming)]
  }else {
    APA_table <-APA_data
  }
  APA_table$P_adj <-p.adjust(APA_table$pvalue, method = "fdr")
  APA_table$Col <- "gray"
  APA_table$Vol_P <- (-log10(APA_table$P_adj))
  APA_table[Vol_P == Inf, Vol_P := (max(-log10(APA_table[P_adj!=0]$P_adj))+0.1)]
  APA_table$APA_trend <- "no change"
  APA_table[(P_adj < P_cutoff & APA_change > delta), Col := "red"]
  APA_table[(P_adj < P_cutoff & APA_change > delta), APA_trend := "longer"]
  APA_table[(P_adj < P_cutoff & APA_change < -delta), Col := "blue"]
  APA_table[(P_adj < P_cutoff & APA_change < -delta), APA_trend := "shorter"]
  shortening <- length(subset(APA_table,Col=="blue")$gene_id)
  lengthening <- length(subset(APA_table,Col=="red")$gene_id)
  all <- length(APA_table$gene_id)
  names <-grep("^PAUs_", colnames(APA_table), value = TRUE)
  sample_names <- sub("^PAUs_", "", names)
#  pdf("APA_plot.pdf",7,7)
  par(mar = c(6, 5, 4, 3) + 0.1)
  
#  dev.off()
  if(split.by.APA_type){
    type_list <-unique(APA_table$APA_type)
    par(mfcol=c(1,length(type_list)))
    for(n in 1:length(type_list)){
      shortening <- length(subset(APA_table[APA_type==type_list[n]],Col=="blue")$gene_id)
      lengthening <- length(subset(APA_table[APA_type==type_list[n]],Col=="red")$gene_id)
      all <- length(APA_table[APA_type==type_list[n]]$gene_id)
      plot(x=APA_table[APA_type==type_list[n]]$APA_change,  y=APA_table[APA_type==type_list[n]]$Vol_P, xlim=c(-1,1), 
           main =paste0("APA trend (",unique(APA_C$APA_type)[n], ")", sep=""),  
           pch= 20, sub = paste(shortening,"shortened","and",lengthening,"lengthened" ,"in all", all, "events",sep=" "),
           xlab = paste0("APA change (",sample_names[2], " - ", sample_names[1], ")", sep=""),  ylab = "-log10(adjusted p value)",  
           col = APA_table[APA_type==type_list[n]]$Col, cex=0.7)
      abline(v = c(-delta,delta), col = "green", lty = 2, lwd = 2)
      abline(h=(-log10(P_cutoff)), col = "green", lty = 2, lwd = 2)
    }
    par(mfcol=c(1,1))
  }else {
    plot(x=APA_table$APA_change,  y=APA_table$Vol_P, xlim=c(-1,1), 
         main =paste0("APA trend (",sample_names[2], " vs ", sample_names[1], ")", sep=""),  
         pch= 20, sub = paste(shortening,"shortened","and",lengthening,"lengthened" ,"in all", all, "events",sep=" "),
         xlab = paste0("APA change (",sample_names[2], " - ", sample_names[1], ")", sep=""),  ylab = "-log10(adjusted p value)",  col = APA_table$Col, cex=0.7)
    abline(v = c(-delta,delta), col = "green", lty = 2, lwd = 2)
    abline(h=(-log10(P_cutoff)), col = "green", lty = 2, lwd = 2)
  }
  return(APA_table[, !c("Col","Vol_P"), with = FALSE])
}