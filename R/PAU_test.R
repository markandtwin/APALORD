#' Perform PAU analysis at single PAS level in each gene 
#' 
#' This function take PAU data output by PAU_by_sample function and compare PAU at single PAS level in each gene between two groups
#' @import dplyr stringr tidyr pbmcapply data.table DEXSeq
#' @param PAU_data  output by PAU_by_sample function
#' @param reads information from RNAseq samples that was used to get the PAU data
#' @param P_cutoff cutoff of adjusted P value to be considered significant
#' @return a table showing the usage of top PASs for each single gene with enough depth in the dataset of each sample 
#' @export

PAU_test <- function(PAU_data,reads, P_cutoff=0.05, FC_cutoff=1){
  groupID <- PAU_data$gene_id
  featureID <- as.character(PAU_data$PAS)
  countData <- PAU_data[, .SD, .SDcols = grep("reads$", names(PAU_data))]
  setnames(countData, gsub(" reads$", "", names(countData)))
  sample_info <-unique(reads[, .SD, .SDcols = (ncol(reads)-1):ncol(reads)])
#  control_sample <- subset(sample_info, treatment==control)$sample
#  experimental_sample <- subset(sample_info, treatment==experimental)$sample
  sampleData <- data.frame(row.names = sample_info$sample,
                           condition = sample_info$treatment)
  dxd <- DEXSeqDataSet(round(countData,0), sampleData, 
                        design= ~ sample + exon + condition:exon, 
                        featureID=featureID, groupID=groupID)
  dxd_norm <- estimateSizeFactors(dxd)
  dxd_est <- estimateDispersions(dxd_norm)
  
  dxr1 <- DEXSeq::DEXSeq(dxd)
  plotMA( dxr1, alpha=P_cutoff, cex=0.8 )
  
  dxrSig <- subset(as.data.frame(dxr1))
  PAU_data[, PAS_name := paste(gene_id, PAS, sep = ":")]
  dxrSig[, "PAS_name"] <-row.names(dxrSig)
  PAU_test_data <- merge(PAU_data[,c("gene_id","chrom","strand","gene_name","PAS","PAS_name" )], dxrSig, by = "PAS_name")
  PAU_test_data <- PAU_test_data[, !c("groupID","featureID","exonBaseMean","dispersion","genomicData"), with = FALSE]
  return(PAU_test_data)
}
