#' Check the 3'end distribution within the same PolyA site.
#' 
#' This function checks the 3'end distribution within individual PAS to detect cleavage site shift.
#' @import dplyr stringr tidyr pbmcapply data.table
#' @param reads information from RNAseq samples
#' @param PAU_data PAU information from previous analysis
#' @param control which group in the data is used as the control group
#' @param experimental which group in the data is used as the experimental group
#' @param cores number of threads used for the computation
#' @param min_reads minium reads count required at PAS for heterogenity analysis
#' @param direct_RNA whether or not the data is direct RNAseq 
#' @return a table showing the uasage of top PASs for each single gene with enough depth in the dataset of each sample 


 CSH_profile<- function(reads,PAU_data, control, experimental,cores=1, min_reads=20,  direct_RNA=F) {
  reads_dt <- as.data.table(reads)
  setkey(reads_dt, gene_id)
  PAU_table <- PAU_data
  PAU_table[,"PAS_name"] <- paste0(PAU_table$gene_id,":",PAU_table$PAS,sep="")
  PASs <- unique(PAU_table$PAS_name)
  setkey(PAU_table, PAS_name)
  
  # Vectorized PAU function
  CSH_fun <- function(PAS_ID) {
    gene <- sub(":.*", "", PAS_ID)
    position <-as.numeric(sub(".*:", "", PAS_ID))
    gene_all <- reads_dt[gene_id==gene]
    strand <- PAU_table[gene_id==gene, "strand"][1]
    # Direct RNA filtering
    if (direct_RNA) {
      gene_all <- gene_all[strand == strand]
    }
    if (strand == "+") {
      PAS_all <- gene_all[abs(chromEnd-position)<=20]
    } else {
      PAS_all <- gene_all[abs((chromStart+1)-position)<=20]
    }
    CSH_PAS <- PAU_table[PAS_name==PAS_ID][,1:5]
    if (all(table(PAS_all$treatment)>=min_reads)) {
        if (strand == "+") {
          control_3end <- PAS_all[treatment == control, chromEnd]
          experimental_3end <- PAS_all[treatment == experimental, chromEnd]
          test_ks_less <- ks.test(control_3end, experimental_3end, alternative = "less")
          test_ks_greater <- ks.test(control_3end, experimental_3end, alternative = "greater")
        } else {
          control_3end <- (PAS_all[treatment == control, chromStart]+1)
          experimental_3end <- (PAS_all[treatment == experimental, chromStart]+1)
          test_ks_less <- ks.test(experimental_3end,control_3end ,alternative = "less")
          test_ks_greater <- ks.test(experimental_3end,control_3end ,alternative = "greater")
        }
        if (test_ks_greater$statistic > test_ks_less$statistic) {
          CSH_PAS$CS_shift <- test_ks_greater$statistic
          CSH_PAS$pvalue <- test_ks_greater$p.value
        } else if (test_ks_greater$statistic < test_ks_less$statistic) {
          CSH_PAS$CS_shift <- (-test_ks_less$statistic)
          CSH_PAS$pvalue <- test_ks_less$p.value
        } else {
          CSH_PAS$CS_shift <- 0
          CSH_PAS$pvalue <- min(test_ks_less$p.value,test_ks_greater$p.value)
        }
      CSH_PAS[,c(paste0("CSH_",control,sep=""),paste0("CSH_",experimental,sep="")):= as.list(c(1-max(table(control_3end))/length(control_3end),1-max(table(experimental_3end))/length(experimental_3end)))]
        return(CSH_PAS)
      }
  }
  
  # Parallel processing
  CSH_output <- pbmclapply(PASs, CSH_fun, mc.cores = cores)
  CSH_table <- rbindlist(CSH_output,fill = T)  # Combine results
  setnames(CSH_table, old = "CS_shift", new = paste0("CS_shift (",experimental," - ", control,")", sep=""))
  return(CSH_table)
}
