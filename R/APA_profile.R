#' Transcriptome wide APA analysis 
#' 
#' This function calls the PASs for each single gene, calculate PAUs for each PAS, and perform the statistics to get the P value and PolyA site usage distance between the selected two groups 
#' @import pbmcapply data.table
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


APA_profile <- function(gene_reference, reads, control, experimental, 
                        min_counts=10, min_reads=5, min_percent=1, 
                        cores=1, direct_RNA=FALSE){
  
  # Filter and group control and experimental data
  control_df <- reads[treatment == control, .SD[.N >= min_counts], by = gene_id]
  experimental_df <- reads[treatment == experimental, .SD[.N >= min_counts], by = gene_id]
  
  # Genes of interest
  genes <- intersect(control_df$gene_id, experimental_df$gene_id)
  genes <- genes[genes != "."]
  
  # Pre-allocate the APA_table
  APA_table <- data.table(gene_id = genes, 
                          short_gene_id = gsub("\\.\\d+$", "", genes), 
                          APA_change = NA, pvalue = NA, APA_type = "No APA")
  
  # Merge APA_table with gene_info
  gene_info <-gene_reference
  APA_table_name <- merge(APA_table, gene_info, by = "gene_id", all.x = TRUE)
  setkey(APA_table_name, gene_id)  # Set key for faster subsetting
  
  # Define PAS_fun with data.table operations for speed
  PAS_fun <- function(df_3end, min_percent=1, min_reads=5){
    if (length(df_3end) >= min_reads) {
      freq_table <- table(df_3end)
      density_vals <- as.numeric(prop.table(freq_table))
      peak_indices <- which(density_vals >= min(sort(density_vals, decreasing = TRUE)[1:50], na.rm = TRUE))
      
      peaks_df <- data.table(Value = as.numeric(names(freq_table))[peak_indices], 
                             Frequency = as.numeric(freq_table[peak_indices]), 
                             Density = 100 * density_vals[peak_indices])
      peaks_df <- peaks_df[order(Value)]
      
      df <- peaks_df[order(peaks_df$Value), ]
      group <- cumsum(c(TRUE, diff(df$Value) > 20))  # Group peaks
      
      # Combine peak groups
      combined_df <- rbindlist(lapply(split(df, group), function(group_data) {
        max_density_row <- group_data[which.max(group_data$Density), ]
        max_density_row$Frequency <- sum(group_data$Frequency)
        max_density_row$Density <- sum(group_data$Density)
        return(max_density_row)
      }))
      
      call_df <- combined_df[combined_df$Density >= min_percent & combined_df$Frequency > 2, ]
      call_df <- call_df[order(call_df$Density, decreasing = TRUE), ]
      
      if (nrow(call_df) > 0) {
        call_df[, Density := round(Density, 2)]
        PAS_gene_table <- as.vector(c(length(call_df$Value),paste(call_df$Value, collapse = ","),
                            paste(call_df$Frequency, collapse = ","),paste(call_df$Density, collapse = ",")))
        return(PAS_gene_table)
      }
    }
  }
  
  # Define APA_fun optimized with data.table
  APA_fun <- function(gene) {
    gene_all <- reads[gene_id == gene]
    APA_gene <- APA_table_name[gene_id == gene]
    
    if (direct_RNA == T) {
      gene_all <- gene_all[strand == gene_info[gene, strand]]
    }
    
    if (all(table(gene_all$sample)>=min_counts)) {
      strand <- gene_info[gene, strand]
      if (strand == "+") {
        control_3end <- gene_all[treatment == control, chromEnd]
        experimental_3end <- gene_all[treatment == experimental, chromEnd]
        test_ks_less <- ks.test(control_3end, experimental_3end, alternative = "less")
        test_ks_greater <- ks.test(control_3end, experimental_3end, alternative = "greater")
        APA_gene$pvalue <- min(test_ks_greater$p.value, test_ks_less$p.value)
        APA_gene$distance <- max(test_ks_greater$statistic, test_ks_less$statistic)
      } else {
        control_3end <- gene_all[treatment == control, chromStart]
        experimental_3end <- gene_all[treatment == experimental, chromStart]
        test_ks_less <- ks.test(experimental_3end,control_3end ,alternative = "less")
        test_ks_greater <- ks.test(experimental_3end,control_3end ,alternative = "greater")
        APA_gene$pvalue <- min(test_ks_greater$p.value, test_ks_less$p.value)
        APA_gene$distance <- max(test_ks_greater$statistic, test_ks_less$statistic)
      }
      if (test_ks_greater$p.value < test_ks_less$p.value) {
        APA_gene$APA_change <- APA_gene$distance
      } else if (test_ks_greater$p.value > test_ks_less$p.value) {
        APA_gene$APA_change <- -APA_gene$distance
      } else {
        APA_gene$APA_change <- 0
      }
      
      df_3end <- c(control_3end, experimental_3end)
      PAS_info <- PAS_fun(df_3end, min_percent, min_reads)
      
      if (!is.null(PAS_info)) {
        PASs_gene <- as.numeric(unlist(strsplit(PAS_info[2], split = ",")))
        PAUs_control <- round(sapply(PASs_gene, function(x) mean(abs(control_3end - x) <= 20) * 100), 2)
        PAUs_experimental <- round(sapply(PASs_gene, function(x) mean(abs(experimental_3end - x) <= 20) * 100), 2)
        PAU_changes <- round(PAUs_experimental - PAUs_control, 2)
        APA_gene[,c("number_of_PAS","PAS_coordinates","PAS_read_counts","PAS_PAUs") := as.list(PAS_info)]
        APA_gene[, c("PAUs_control", "PAUs_experimental", "PAU_changes") := as.list(c(paste(PAUs_control, collapse = ","),
                                                                                      paste(PAUs_experimental, collapse = ","),paste(PAU_changes, collapse = ",")))]
        if (as.numeric(PAS_info[1]) > 1) {
          if((strand=="-")&&(PASs_gene[length(PASs_gene)]<APA_gene[,"chromEnd"])){
            APA_gene[,"APA_type"] <- "Last_exon_tandem_APA"
          }
          else if((strand=="+")&&(PASs_gene[1]>APA_gene[,"chromStart"])){
            APA_gene[,"APA_type"] <- "Last_exon_tandem_APA"
          }
          else {
            APA_gene[,"APA_type"] <- "Mixed_APA"
          }
        }
      }
      return(APA_gene)
    }
  }
  
  # Parallelize APA analysis using pbmclapply
  output <- pbmclapply(genes, APA_fun, mc.cores = cores)
  
  # Combine the results
  output_df <- rbindlist(output)
  output_df <- output_df[, !c("chromStart", "chromEnd","distal_stop_codon_Start","distal_stop_codon_End","distance"), with = FALSE]
  return(output_df)
}
