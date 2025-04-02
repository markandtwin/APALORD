#' Transcriptome wide APA analysis 
#' 
#' This function calls the PASs for each single gene, calculate PAUs for each PAS, and perform the statistics to get the P value and PolyA site usage distance between the selected two groups 
#' @import pbmcapply data.table
#' @param gene_reference information extracted from gtf file
#' @param reads information from samples
#' @param control which group in the data is used as the control group
#' @param experimental which group in the data is used as the experimental group
#' @param direct_RNA whether the data is direct RNAseq, TRUE or FALSE 
#' @param min_counts minium read counts required for both groups at single gene level to perform the analysis
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
      # Group peaks
      if (nrow(df)>1){
        group <- integer(nrow(df))
        current_group <- 1
        group[1] <- current_group
        rep_val <- df$Value[1]
        # Assign groups: start a new group when a value exceeds min_val + threshold
        for (i in 2:nrow(df)) {
          group_df <- df[which(group==current_group),]
          rep_val <- group_df[which.max(group_df$Frequency),"Value"]
          if (df$Value[i] - rep_val > 20) {
            current_group <- current_group + 1
          }
          group[i] <- current_group
        }
        
        
        # Combine peak groups
        combined_df <- rbindlist(lapply(split(df, group), function(group_data) {
          max_density_row <- group_data[which.max(group_data$Density), ]
          max_density_row$Frequency <- sum(group_data$Frequency)
          max_density_row$Density <- sum(group_data$Density)
          return(max_density_row)
        }))
      } else {
        combined_df <-df
      }
      
      call_df <- combined_df[combined_df$Density >= min_percent & combined_df$Frequency > 2, ]
      call_df <- call_df[order(call_df$Value, decreasing = FALSE), ]
      
      if (nrow(call_df)>1){
        marks<-vector()
        for (i in 1:(nrow(call_df)-1)) {
          diff <- call_df[i+1,]$Value - call_df[i,]$Value  # Since values are ordered, subtraction works directly
          if (diff < 20*2) {
            m <- ifelse(call_df[i,]$Frequency < call_df[i+1,]$Frequency, i, i+1)
            marks <- append(marks, m)
          }
        }
        if (length(marks)>0){
          call_df <- call_df[-marks,]
        }
      }
      
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
    
    if (all(table(gene_all$treatment)>=min_counts)) {
      strand <- gene_info[gene, strand]
      if (strand == "+") {
        control_3end <- gene_all[treatment == control, chromEnd]
        experimental_3end <- gene_all[treatment == experimental, chromEnd]
      } else {
        control_3end <- (gene_all[treatment == control, chromStart]+1)
        experimental_3end <- (gene_all[treatment == experimental, chromStart]+1)
      }
      
      
      df_3end <- c(control_3end, experimental_3end)
      PAS_info <- PAS_fun(df_3end, min_percent, min_reads)
      if (is.null(PAS_info)){
        n_PAS <- 0
      } else {
        n_PAS <-length(as.numeric(unlist(strsplit(PAS_info[2], ","))))
      }
      
      if ((n_PAS>=2)&(all(table(gene_all$sample)>=min_reads))) {
        PASs_gene <- as.numeric(unlist(strsplit(PAS_info[2], split = ",")))
        for (PAS in PASs_gene){
          control_3end[abs(control_3end - PAS) <= 20] <- PAS
          experimental_3end[abs(experimental_3end - PAS) <= 20] <- PAS
        }
        reads_control <- table(control_3end[control_3end %in% PASs_gene])
        reads_control <- c(reads_control, setNames(rep(0, length(PASs_gene) - length(reads_control)), PASs_gene[!PASs_gene %in% names(reads_control)]))
        reads_control <-  reads_control[order(names(reads_control))]
        PAUs_control <- round(100*reads_control/length(control_3end),2)
        reads_experimental <- table(experimental_3end[experimental_3end %in% PASs_gene])
        reads_experimental <- c(reads_experimental, setNames(rep(0, length(PASs_gene) - length(reads_experimental)), PASs_gene[!PASs_gene %in% names(reads_experimental)]))
        reads_experimental <-  reads_experimental[order(names(reads_experimental))]
        PAUs_experimental <- round(100*reads_experimental/length(experimental_3end),2)
        PAU_changes <- round(PAUs_experimental - PAUs_control, 2)
        APA_gene[,c("number_of_PAS","PAS_coordinates","PAS_read_counts","PAS_PAUs") := as.list(PAS_info)]
        APA_gene[, c(paste0("PAUs_",control,sep=""), paste0("PAUs_",experimental,sep=""), "PAU_changes") := as.list(c(paste(PAUs_control, collapse = ","),
                                                                                      paste(PAUs_experimental, collapse = ","),paste(PAU_changes, collapse = ",")))]
        APA_gene[, c(paste0("reads_",control,sep=""), paste0("reads_",experimental,sep="")) := as.list(c(paste(reads_control, collapse = ","),paste(reads_experimental, collapse = ",")))]
        if (strand == "+") {
          test_ks_less <- ks.test(control_3end, experimental_3end, alternative = "less")
          test_ks_greater <- ks.test(control_3end, experimental_3end, alternative = "greater")
        } else {
          test_ks_less <- ks.test(experimental_3end,control_3end ,alternative = "less")
          test_ks_greater <- ks.test(experimental_3end,control_3end ,alternative = "greater")
        }
        if (test_ks_greater$statistic > test_ks_less$statistic) {
          APA_gene$APA_change <- test_ks_greater$statistic
          APA_gene$pvalue <- test_ks_greater$p.value
        } else if (test_ks_greater$statistic < test_ks_less$statistic) {
          APA_gene$APA_change <- (-test_ks_less$statistic)
          APA_gene$pvalue <- test_ks_less$p.value
        } else {
          APA_gene$APA_change <- 0
          APA_gene$pvalue <- min(test_ks_less$p.value,test_ks_greater$p.value)
        }
        if((strand=="-")&&(PASs_gene[length(PASs_gene)]<APA_gene[,"last_exon_chromEnd"])){
          APA_gene[,"APA_type"] <- "Last_exon_tandem_APA"
        } else if((strand=="+")&&(PASs_gene[1]>APA_gene[,"last_exon_chromStart"])){
          APA_gene[,"APA_type"] <- "Last_exon_tandem_APA"
        } else {
          APA_gene[,"APA_type"] <- "Mixed_APA"
        }
        return(APA_gene)
      }
    }
  }
  
  # Parallelize APA analysis using pbmclapply
  output <- pbmclapply(genes, APA_fun, mc.cores = cores)
  
  # Combine the results
  output_df <- rbindlist(output, fill = T)
  output_df$number_of_PAS <- as.numeric(output_df$number_of_PAS)
  output_df <- output_df[, !c("last_exon_chromStart", "last_exon_chromEnd"), with = FALSE]
  return(output_df)
}
