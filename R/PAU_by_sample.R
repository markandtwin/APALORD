#' Call PASs at single gene level and calculate the usage of each found PAS
#' 
#' This function extract PolyA sites information of each gene and calculate the usage of each found PAS from the provided RNAseq data.
#' @param gene_reference information extracted from gtf file
#' @param reads information from RNAseq samples
#' @param cores number of threads used for the computation
#' @param min_reads minium reads count required at single gene level for PAS calling and PAU calculation per sample
#' @param min_percent minium percent required for a PAS to be included (0-100)
#' @param direct_RNA whether or not the data is direct RNAseq 
#' @param internal_priming whether or not PASs subject to internal priming filtering
#' @param pattern to look for internal priming sequencing surrouding PAS, either "pre", "post" or "both".
#' @param genome_file path to a genome sequence file (.fa or .fasta)
#' @param bin max distance of cleavage site from PAS
#' @return a table showing the uasage of top PASs for each single gene with enough depth in the dataset of each sample 
#' @export

PAU_by_sample <- function(gene_reference, reads, min_reads=5, min_percent=1, cores=1, direct_RNA=F,internal_priming=F,pattern="post",genome_file=NULL,bin=50) {
  gene_info <- gene_reference[,c(1:2, 5:6)]
  if(internal_priming){
    genome <- Rsamtools::FaFile(genome_file)
    open(genome)
    genome
  }
  reads$read_id <-c(1:nrow(reads))
  reads_dt <- data.table::as.data.table(reads[, !("treatment"), with = FALSE])
  data.table::setkey(reads_dt, gene_id)

  
  # Pre-filter the reads
  reads_dt <- reads_dt[gene_id%in%unique(reads_dt[,.N, by = gene_id][N >= min_reads]$gene_id)]  # Filter genes with fewer than min_reads
  genes <- unique(reads_dt$gene_id)
  genes <- as.vector(genes[genes != "."])
  samples <- unique(reads_dt$sample)
  
  
  PAS_table <- gene_info[gene_id %in% genes]
  data.table::setkey(PAS_table, gene_id)  # Ensure efficient subsetting
  
  find_match <- function(a, b_vec) {
    distances <- abs(a - b_vec)
    if (any(distances <= bin)) {
      return(b_vec[which.min(distances)])  # Return closest match
    } else {
      return(NA)  # No match within bin distance
    }
  }
  
  match_fun <- function(strand,reads_data,call_df){
    B <- call_df$Value
    if(strand=="+"){
      A <- reads_data$chromEnd
      matches <- sapply(A, find_match, b_vec = B)
      reads_data$match <-matches
    }
    if(strand=="-"){
      A <- reads_data$chromStart+1
      matches <- sapply(A, find_match, b_vec = B)
      reads_data$match <-matches
    }
    PAS_sample_counts <- table(reads_data$match)
    return(reads_data)
  }
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
          if (df$Value[i] - rep_val > bin) {
            current_group <- current_group + 1
          }
          group[i] <- current_group
        }
        
        
        # Combine peak groups
        combined_df <-  data.table::rbindlist(lapply(split(df, group), function(group_data) {
          max_density_row <- group_data[which.max(group_data$Density), ]
          max_density_row$Frequency <- sum(group_data$Frequency)
          max_density_row$Density <- sum(group_data$Density)
          return(max_density_row)
        }))
      } else {
        combined_df <-df
      }
      call_df <- combined_df[combined_df$Density >= min_percent & combined_df$Frequency >=min_reads, ]
      call_df <- call_df[order(call_df$Value, decreasing = F), ]
      
      
      if (nrow(call_df) > 0) {
        call_df$Density <- round(call_df$Density, 2)
      # Parallelize per sample
        if (all(table(factor(gene_all$sample,levels = samples))>=min_reads)&nrow(call_df) > 0) {
          PAS_sample_counts_raw <- data.table::as.data.table(match_fun(strand,reads_data =  gene_all,call_df))
          tab <- table(PAS_sample_counts_raw$match)
          keep <- names(tab[tab >= min_reads])
          PAS_sample_counts <- PAS_sample_counts_raw[match %in% keep]
          number_of_PAS <- length(names(table(PAS_sample_counts$match)))
          PAS_gene <- PAS_table[gene, 1:4][rep(1, number_of_PAS),]
          PAS_gene$PAS <- names(table(PAS_sample_counts$match))
          all_matches <- PAS_gene$PAS
          for (sample_name in samples) {
            sample_depth <-nrow(PAS_sample_counts_raw[sample == sample_name]) 
            PAS_gene[,paste(sample_name,"reads")] <- as.integer(table(factor(PAS_sample_counts[sample == sample_name]$match,levels = all_matches)))
            PAS_gene[, paste(sample_name,"PAU")] <- as.numeric(round(100*PAS_gene[,get(paste(sample_name,"reads"))]/sample_depth,3))
            PAS_gene[, paste(sample_name,"host gene reads",sep=" ")] <- sample_depth
          }
          if (direct_RNA==F&internal_priming){
            clean_PAS_gene <- Internal_priming(PAS_gene,pattern=pattern,genome=genome)
            PAS_gene$internal_priming <- ifelse(paste(PAS_gene$gene_id, PAS_gene$PAS) %in%
                                                  paste(clean_PAS_gene$gene_id, clean_PAS_gene$PAS),"no", "yes")
          }
          return(PAS_gene)
        }
      }
    }
  }
  
  # Parallel processing
  PAU_output <- pbmcapply::pbmclapply(genes, PAU_fun, mc.cores = cores)
  PAU_table <-  data.table::rbindlist(PAU_output,fill = T)  # Combine results
  PAU_table$PAS <-  as.numeric(PAU_table$PAS)
  return(PAU_table)
}
