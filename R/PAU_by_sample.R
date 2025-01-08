#' Call PASs at single gene level and calculate the usage of each found PAS
#' 
#' This function extract PolyA sites information of each gene and calculate the usage of each found PAS from the provided RNAseq data.
#' @import dplyr stringr tidyr pbmcapply data.table
#' @param gene_reference information extracted from gtf file
#' @param reads information from RNAseq samples
#' @param cores number of threads used for the computation
#' @param min_reads minium reads count required at single gene level for PAS calling
#' @param min_percent minium percent required for a PAS to be included (0-100)
#' @param direct_RNA whether or not the data is direct RNAseq 
#' @return a table showing the top PASs for each single gene with enough depth in the data and also a bed6 file table for all the PASs
#' @export

PAU_by_sample <- function(gene_reference, reads, min_reads=5, min_percent=1, cores=1, direct_RNA=F) {
  gene_info <- gene_reference[,c(1:2, 5:6)]
  reads_dt <- as.data.table(reads)
  setkey(reads_dt, gene_id)
  
  # Pre-filter the reads
  reads_dt <- reads_dt[unique(reads_dt[,.N, by = gene_id][N >= min_reads]$gene_id)]  # Filter genes with fewer than min_reads
  genes <- unique(reads_dt$gene_id)
  genes <- as.vector(genes[genes != "."])
  
  
  PAS_table <- gene_info[gene_id %in% genes]
  setkey(PAS_table, gene_id)  # Ensure efficient subsetting
  
  # Vectorized PAU function
  PAU_fun <- function(gene) {
    gene_all <- reads_dt[gene]
    
    # Direct RNA filtering
    if (direct_RNA) {
      strand <- PAS_table[gene, strand]
      gene_all <- gene_all[strand == strand]
    }
    
    # Handle strand-specific processing
    strand <- PAS_table[gene, strand]
    if (strand == "+") {
      df_3end <- gene_all$chromEnd
    } else {
      df_3end <- gene_all$chromStart
    }
    
    # Calculate frequency table and density
    if (length(df_3end) >= min_reads) {
      freq_table <- table(df_3end)
      density_vals <- as.numeric(prop.table(freq_table))
      top_peaks <- density_vals >= min(sort(density_vals, decreasing = TRUE)[1:50], na.rm = TRUE)
      
      peaks_df <- data.frame(
        Value = as.numeric(names(freq_table))[top_peaks],
        Frequency = as.numeric(freq_table[top_peaks]),
        Density = 100 * density_vals[top_peaks]
      )
      
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
        call_df$Density <- round(call_df$Density, 4)
        PAS_gene <- PAS_table[gene, 1:4][rep(1, nrow(call_df)),]
        PAS_gene$PAS <- call_df$Value
        
        # Parallelize per sample
        for (sample in unique(reads_dt$sample)) {
          sample_df <- gene_all[sample == sample]
          PAS_gene[, sample] <- round(sapply(call_df$Value, function(x) {
            mean(abs(sample_df$chromEnd - x) <= 20) * 100
          }), 2)
        }
        return(PAS_gene)
      }
    }
  }
  
  # Parallel processing
  PAU_output <- pbmclapply(genes, PAU_fun, mc.cores = cores)
  PAU_table <- rbindlist(PAU_output)  # Combine results
  
  return(PAU_table)
}
