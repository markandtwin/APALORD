#' Visualize the 3'ends of all the reads in each sample for an individual gene
#' 
#' This function loads a gtf file and extract necessary information of each gene from it.
#' @import dplyr stringr tidyr pbmcapply
#' @param gene_reference information extracted from gtf file
#' @param reads information from samples
#' @param data_type if the data is direct RNA
#' @param min minium read counts required for both samples at single gene level
#' @param cores number of threads used for the computation
#' @return a table showing the APA change for each single gene with enough depth in the data
#' @export

APA_profile <- function(gene_reference, reads, min=10,cores=1,data_type="direct RNA"){
  gene_info <-gene_reference[[1]]
  Undiff_df <- reads %>%filter(treatment == "Undiff") %>%group_by(gene_id) %>%filter(n() >= min) %>%ungroup()
  Diff_df <- reads %>%filter(treatment == "Diff") %>%group_by(gene_id) %>%filter(n() >= min) %>%ungroup()
  
  genes <- intersect(Undiff_df$gene_id, Diff_df$gene_id)
  genes <- subset(genes, genes!=".")
  APA_table<- data.frame(matrix(ncol=4, nrow=(length(genes))))
  colnames(APA_table) <- c("gene_id","short_gene_id","APA_change","pvalue")
  APA_table$gene_id <- genes
  APA_table$short_gene_id <- gsub("\\.\\d+$", "", APA_table$gene_id)
  APA_table_name <- left_join (APA_table, gene_info, by = "gene_id")
  row.names(APA_table_name) <- APA_table_name$gene_id
  APA_table_name<- subset(APA_table_name, !is.na(strand))
  genes<- APA_table_name$gene_id
  
  APA_fun <-function(gene,reads,min, APA_table_name,data_type){
    gene_all <-subset(reads, gene_id==gene)
    if(data_type=="direct RNA"){
      gene_all <- subset(gene_all, strand==gene_info[gene,"strand"])
    }
    if(APA_table_name[gene,"strand"]=="+"){
      Undiff <-subset(gene_all, treatment=="Undiff")$chromEnd
      Diff <-subset(gene_all, treatment=="Diff")$chromEnd
      if((length(Undiff)>=min)&&(length(Diff)>=min)){
        test_ks_less <-ks.test(Undiff, Diff, alternative = "less")
        test_ks_greater <- ks.test(Undiff, Diff, alternative = "greater")
        APA_table_name[gene,"pvalue"] <- min(test_ks_greater$p.value,test_ks_less$p.value)
        APA_table_name[gene,"distance"] <- max(test_ks_greater$statistic,test_ks_less$statistic)
        if (test_ks_greater$p.value<test_ks_less$p.value){
          APA_table_name[gene,"APA_change"] <- APA_table_name[gene,"distance"]
        } 
        if (test_ks_greater$p.value>test_ks_less$p.value){
          APA_table_name[gene,"APA_change"] <- APA_table_name[gene,"distance"]*(-1)
        } 
        if (test_ks_greater$p.value==test_ks_less$p.value){
          APA_table_name[gene,"APA_change"] <-0
        }
        return(APA_table_name[gene,])
      }
    }
    if(APA_table_name[gene,"strand"]=="-"){
      Undiff <-subset(gene_all, treatment=="Undiff")$chromStart
      Diff <-subset(gene_all, treatment=="Diff")$chromStart
      if((length(Undiff)>=min)&&(length(Diff)>=min)){
        test_ks_less <-ks.test(Diff, Undiff, alternative = "less")
        test_ks_greater <- ks.test(Diff, Undiff, alternative = "greater")
        APA_table_name[gene,"pvalue"] <- min(test_ks_greater$p.value,test_ks_less$p.value)
        APA_table_name[gene,"distance"] <- max(test_ks_greater$statistic,test_ks_less$statistic)
        if (test_ks_greater$p.value<test_ks_less$p.value){
          APA_table_name[gene,"APA_change"] <- APA_table_name[gene,"distance"]
        } 
        if (test_ks_greater$p.value>test_ks_less$p.value){
          APA_table_name[gene,"APA_change"] <- APA_table_name[gene,"distance"]*(-1)
        } 
        if (test_ks_greater$p.value==test_ks_less$p.value){
          APA_table_name[gene,"APA_change"] <-0
        }
        return(APA_table_name[gene,])
      }
    }
  }
  output <-pbmclapply(genes, APA_fun, reads=reads, min=min, APA_table_name=APA_table_name,
                      data_type=data_type, mc.cores = cores)
  output_df <- do.call(rbind, output)
  return(output_df)
  
}