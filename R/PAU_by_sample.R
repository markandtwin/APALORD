#' Call PASs at single gene level and calculate the usage of each found PAS
#' 
#' This function extract PolyA sites information of each gene and calculate the usage of each found PAS from the provided RNAseq data.
#' @import dplyr stringr tidyr pbmcapply data.table 
#' @param gene_reference information extracted from gtf file
#' @param reads information from RNAseq samples
#' @param cores number of threads used for the computation
#' @param min_reads minium reads count required at single gene level for PAS calling and PAU calculation
#' @param min_percent minium percent required for a PAS to be included (0-100)
#' @param direct_RNA whether or not the data is direct RNAseq 
#' @param internal_priming whether or not PASs subject to internal priming filtering
#' @param pattern to look for internal priming sequencing surrouding PAS, either "pre", "post" or "both".
#' @param genome_file path to a genome sequence file (.fa or .fasta)
#' @return a table showing the uasage of top PASs for each single gene with enough depth in the dataset of each sample 
#' @export

PAU_by_sample <- function(gene_reference, reads, min_reads=5, min_percent=1, cores=1, direct_RNA=F,internal_priming=F,pattern="post",genome_file=NULL) {
  gene_info <- gene_reference[,c(1:2, 5:6)]
  if(internal_priming){
    genome <- FaFile(genome_file)
    open(genome)
    genome
  }
  reads_dt <- as.data.table(reads)
  setkey(reads_dt, gene_id)
  
  # Pre-filter the reads
  reads_dt <- reads_dt[gene_id%in%unique(reads_dt[,.N, by = gene_id][N >= min_reads]$gene_id)]  # Filter genes with fewer than min_reads
  genes <- unique(reads_dt$gene_id)
  genes <- as.vector(genes[genes != "."])
  
  
  PAS_table <- gene_info[gene_id %in% genes]
  setkey(PAS_table, gene_id)  # Ensure efficient subsetting
  
  # Vectorized PAU function
  PAU_fun <- function(gene) {
    gene_all <- reads_dt[gene]
    
    # Direct RNA filtering
    if (direct_RNA) {
      gene_all <- gene_all[strand == gene_info[gene, strand]]
    }
    
    # Handle strand-specific processing
    strand <- PAS_table[gene, strand]
    if (strand == "+") {
      df_3end <- gene_all$chromEnd
    } else {
      df_3end <- (gene_all$chromStart+1)
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
      call_df <- call_df[order(call_df$Value, decreasing = F), ]
      
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
        call_df$Density <- round(call_df$Density, 2)
        PAS_gene <- PAS_table[gene, 1:4][rep(1, nrow(call_df)),]
        PAS_gene$PAS <- call_df$Value
        if (direct_RNA==F&internal_priming){PAS_gene <- Internal_priming(PAS_gene,pattern=pattern)}
        if (nrow(PAS_gene)<nrow(call_df)){
#         False_df <- as.data.table(call_df)[Value!%in%PAS_gene$PAS]
          call_df <- as.data.table(call_df)[Value%in%PAS_gene$PAS]
        }
        samples <- unique(reads_dt$sample)
        
        # Parallelize per sample
        if (all(table(gene_all$sample)>=min_reads)&nrow(PAS_gene) > 0) {
          for (sample_name in samples) {
            sample_df <- gene_all[sample == sample_name]
            if (strand == "+") {
              PAS_gene[, paste(sample_name,"PAU")] <- round(sapply(call_df$Value, function(x) {
                mean(abs(sample_df$chromEnd - x) <= 20) * 100
              }), 2)
              PAS_gene[, paste(sample_name,"reads")] <- sapply(call_df$Value, function(x) {
                sum(abs(sample_df$chromEnd - x) <= 20)})
            } else {
              PAS_gene[, paste(sample_name,"PAU")] <- round(sapply(call_df$Value, function(x) {
                mean(abs(sample_df$chromStart - x) <= 20) * 100
              }), 2)
              PAS_gene[, paste(sample_name,"reads")] <- sapply(call_df$Value, function(x) {
                sum(abs(sample_df$chromStart - x) <= 20)})
            }
          }
          return(PAS_gene)
        }
      }
    }
  }
  
  # Parallel processing
  PAU_output <- pbmclapply(genes, PAU_fun, mc.cores = cores)
  PAU_table <- rbindlist(PAU_output,fill = T)  # Combine results
  
  return(PAU_table)
}
