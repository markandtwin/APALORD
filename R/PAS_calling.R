#' Call PASs at single gene level
#' 
#' This function extract PolyA sites information of each gene from the provided RNAseq data.
#' @import dplyr stringr tidyr pbmcapply
#' @param gene_reference information extracted from gtf file
#' @param reads information from RNAseq samples
#' @param cores number of threads used for the computation
#' @param min_reads minium read counts required at single gene level for PAS calling
#' @param min_percent minium percent required for a PAS to be included (0-100)
#' @param direct_RNA whether or not the data is direct RNAseq 
#' @return a table showing the top PASs for each single gene with enough depth in the data and also a bed6 file table for all the PASs
#' @export

PAS_calling <- function(gene_reference, reads,min_reads=5, min_percent=1,cores=1,direct_RNA=FALSE){
  gene_info <- gene_reference[,.(gene_id, strand, gene_name, gene_biotype)]
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
        call_df$Density <- round(call_df$Density, 2)
        PAS_gene <- PAS_table[gene, 1:4][rep(1, nrow(call_df)),]
        PAS_gene$PAS <- call_df$Value
        add_df<- data.frame(matrix(ncol=6, nrow=(length(call_df$Value))))
        if(strand=="+"){
          add_df<-data.frame(
            chrom = gene_all$chrom[1], # Numeric column
            start = as.numeric(call_df$Value)-1, # Numeric column
            end = as.numeric(call_df$Value),
            name = gene_info[gene,gene_name],
            density= as.numeric(call_df$Density),
            strand = gene_info[gene,"strand"])
        }
        if(strand=="-"){
          add_df<-data.frame(
            chrom = gene_all$chrom[1], # Numeric column
            start = as.numeric(call_df$Value), # Numeric column
            end = as.numeric(call_df$Value)+1,
            name = gene_info[gene,gene_name],
            density= as.numeric(call_df$Density),
            strand = strand)
        }
        return(add_df)
      }
    }
  }
  
  # Parallel processing
  PAS_output <- pbmclapply(genes, PAU_fun, mc.cores = cores)
  PAS_bed <- rbindlist(PAS_output,fill = T)  # Combine results
  
  return(PAS_bed)
}
