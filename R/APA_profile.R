#' Transcriptome wide APA analysis 
#' 
#' This function calls the PASs for each single gene, calculate PAUs for each PAS, and perform the statistics to get the P value and PolyA site usage distance between the selected two groups 
#' @import dplyr stringr tidyr pbmcapply
#' @param gene_reference information extracted from gtf file
#' @param reads information from samples
#' @param control which group in the data is used as the control group
#' @param experimental which group in the data is used as the experimental group
#' @param direct_RNA whether the data is direct RNAseq, TRUE or FALSE 
#' @param min_counts minium read counts required for both samples at single gene level to perform the analysis
#' @param min_reads minium read counts required for PAS calling
#' @param min_percent  minium PAU for a PAS to be considered for downstream analysis
#' @param cores number of threads used for the computation
#' @return a table showing the APA change for each single gene with enough depth in the data
#' @export

APA_profile <- function(gene_reference, reads, control,experimental, min_counts=10,min_reads=5, min_percent=1, cores=1,direct_RNA=FALSE){
  gene_info <-gene_reference
  control_df <- reads %>%filter(treatment == control) %>%group_by(gene_id) %>%filter(n() >= min_counts) %>%ungroup()
  experimental_df <- reads %>%filter(treatment == experimental) %>%group_by(gene_id) %>%filter(n() >= min_counts) %>%ungroup()
  
  genes <- intersect(control_df$gene_id, experimental_df$gene_id)
  genes <- subset(genes, genes!=".")
  APA_table<- data.frame(matrix(ncol=4, nrow=(length(genes))))
  colnames(APA_table) <- c("gene_id","short_gene_id","APA_change","pvalue")
  APA_table$gene_id <- genes
  APA_table$short_gene_id <- gsub("\\.\\d+$", "", APA_table$gene_id)
  APA_table_name <- left_join (APA_table, gene_info, by = "gene_id")
  row.names(APA_table_name) <- APA_table_name$gene_id
  APA_table_name<- subset(APA_table_name, !is.na(strand))
  genes<- APA_table_name$gene_id
  APA_table_name$APA_type<- "No APA"

  
  PAS_fun <- function (df_3end, min_percent=1,min_reads=5){
    if(length(df_3end)>=min_reads){
      freq_table <- table(df_3end)
      density_vals <- as.numeric(prop.table(freq_table))
      peak_indices <- which(density_vals>=min(sort(density_vals, decreasing = T)[1:50],na.rm = T))
      peaks_df <- data.frame(
        Value = as.numeric(names(freq_table))[peak_indices],
        Frequency = as.numeric(freq_table[peak_indices]),
        Density = 100*density_vals[peak_indices]
      )
      df <- peaks_df[order(peaks_df$Value), ]
      group <- cumsum(c(TRUE, diff(df$Value) > 20))
      
      combined_df <- do.call(rbind, lapply(split(df, group), function(group_data) {
        max_density_row <- group_data[which.max(group_data$Density), ]
        max_density_row$Frequency <- sum(group_data$Frequency)
        max_density_row$Density <- sum(group_data$Density)
        return(max_density_row)
      }))
      call_df <- combined_df[which((combined_df$Density>=min_percent)&(combined_df$Frequency>2)),]
      
      # Display the peaks
      call_df<-(call_df[order(call_df$Value, decreasing = F),])
      if(length(call_df$Frequency)!=0){
        call_df$Density <-round(as.numeric(call_df$Density), 4)
        PAS_gene_table <- c(length(call_df$Value),paste(call_df$Value, collapse = ","),
                                    paste(call_df$Frequency, collapse = ","),paste(call_df$Density, collapse = ","))
        return(PAS_gene_table)
      }
    }
  }
  
  APA_fun <-function(gene,reads,min_counts,min_reads,min_percent,APA_table_name,direct_RNA=FALSE){
    gene_all <-subset(reads, gene_id==gene)
    APA_gene <- APA_table_name[gene,]
    if(direct_RNA=="direct RNA"){
      gene_all <- subset(gene_all, strand==gene_info[gene,"strand"])
    }
    if(APA_gene[,"strand"]=="+"){
      control_3end <-subset(gene_all, treatment==control)$chromEnd
      experimental_3end <-subset(gene_all, treatment==experimental)$chromEnd
      if((length(control_3end)>=min_counts)&&(length(experimental_3end)>=min_counts)){
        test_ks_less <-ks.test(control_3end, experimental_3end, alternative = "less")
        test_ks_greater <- ks.test(control_3end, experimental_3end, alternative = "greater")
        APA_gene[,"pvalue"] <- min(test_ks_greater$p.value,test_ks_less$p.value)
        APA_gene[,"distance"] <- max(test_ks_greater$statistic,test_ks_less$statistic)
        if (test_ks_greater$p.value<test_ks_less$p.value){
          APA_gene[,"APA_change"] <- APA_gene[,"distance"]
        } 
        if (test_ks_greater$p.value>test_ks_less$p.value){
          APA_gene[,"APA_change"] <- APA_gene[,"distance"]*(-1)
        } 
        if (test_ks_greater$p.value==test_ks_less$p.value){
          APA_gene[,"APA_change"] <-0
        }
        df_3end <- c(control_3end,experimental_3end)
        APA_gene[,c("number","PASs","PASs reads","PASs fraction")]<-PAS_fun(df_3end,min_percent=min_percent,min_reads=min_reads)
        PASs_gene <- as.numeric(unlist(strsplit(APA_gene[,]$PASs, split = ",")))
        PAUs_control <-round(sapply(PASs_gene, function(x) {
          mean(abs(control_3end-x)<=20) * 100}),2)
        PAUs_experimental <-round(sapply(PASs_gene, function(x) {
          mean(abs(experimental_3end-x)<=20) * 100}),2)
        PAU_changes<-round((PAUs_experimental-PAUs_control),2)
        APA_gene[,c("PAUs_control","PAUs_experimental","PAU_changes")]<-c(paste(PAUs_control, collapse = ","),
                                                                          paste(PAUs_experimental, collapse = ","),paste(PAU_changes, collapse = ","))
        return(APA_gene)
      }
    }
    if(APA_gene[,"strand"]=="-"){
      control_3end <-subset(gene_all, treatment==control)$chromStart
      experimental_3end <-subset(gene_all, treatment==experimental)$chromStart
      if((length(control_3end)>=min_counts)&&(length(experimental_3end)>=min_counts)){
        test_ks_less <-ks.test(experimental_3end, control_3end, alternative = "less")
        test_ks_greater <- ks.test(experimental_3end, control_3end, alternative = "greater")
        APA_gene[,"pvalue"] <- min(test_ks_greater$p.value,test_ks_less$p.value)
        APA_gene[,"distance"] <- max(test_ks_greater$statistic,test_ks_less$statistic)
        if (test_ks_greater$p.value<test_ks_less$p.value){
          APA_gene[,"APA_change"] <- APA_gene[,"distance"]
        } 
        if (test_ks_greater$p.value>test_ks_less$p.value){
          APA_gene[,"APA_change"] <- APA_gene[,"distance"]*(-1)
        } 
        if (test_ks_greater$p.value==test_ks_less$p.value){
          APA_gene[,"APA_change"] <-0
        }
        df_3end <- c(control_3end,experimental_3end)
        APA_gene[,c("number","PASs","PASs reads","PASs fraction")]<-PAS_fun(df_3end,min_percent=min_percent,min_reads=min_reads)
        PASs_gene <- as.numeric(unlist(strsplit(APA_gene[,]$PASs, split = ",")))
        PAUs_control <-round(sapply(PASs_gene, function(x) {
          mean(abs(control_3end-x)<=20) * 100}),2)
        PAUs_experimental <-round(sapply(PASs_gene, function(x) {
          mean(abs(experimental_3end-x)<=20) * 100}),2)
        PAU_changes<-round((PAUs_experimental-PAUs_control),2)
        APA_gene[,c("PAUs_control","PAUs_experimental","PAU_changes")]<-c(paste(PAUs_control, collapse = ","),
                                                                          paste(PAUs_experimental, collapse = ","),paste(PAU_changes, collapse = ","))
        return(APA_gene[,])
      }
    }
  }

  output <-pbmclapply(genes, APA_fun, reads=reads, min_counts=min_counts,min_reads=min_reads, min_percent=min_percent,APA_table_name=APA_table_name,
                      direct_RNA=direct_RNA, mc.cores = cores)
  output_df <- do.call(rbind, output)

  for(gene in genes){
    PASs_gene <- sort(as.numeric(unlist(strsplit(output_df[gene,]$PASs, split = ","))))
    if(length(PASs_gene)>1){
      if((output_df[gene,"strand"]=="-")&&(PASs_gene[length(PASs_gene)]<output_df[gene,"last_exon_End"])){
        output_df[gene,"APA_type"] <- "Last_exon_tandem_APA"
      }
      else if((output_df[gene,"strand"]=="+")&&(PASs_gene[1]>output_df[gene,"last_exon_Start"])){
        output_df[gene,"APA_type"] <- "Last_exon_tandem_APA"
      }
      else {
        output_df[gene,"APA_type"] <- "Mixed_APA"
      }
    }
  }
  output_df$P_adjusted
  output_df[,1:(ncol(output_df) - 6)] <- output_df[,1:(ncol(output_df) - 6)] %>%
    mutate_if(is.character, as.factor)
  return(output_df)
}