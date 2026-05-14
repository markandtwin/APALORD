#' Call PASs at single gene level
#' 
#' This function extract PolyA sites information of each gene from the provided RNAseq data.
#' @param gene_reference information extracted from gtf file
#' @param reads information from RNAseq samples 
#' @param cores number of threads used for the computation
#' @param min_reads minimum read counts required at single gene level for PAS calling per sample
#' @param min_percent minimum percent required for a PAS to be included (0-100)
#' @param direct_RNA whether or not the data is direct RNAseq 
#' @param internal_priming whether or not PASs subject to internal priming filtering
#' @param pattern to look for internal priming sequencing surrounding PAS, either "pre", "post" or "both".
#' @param genome_file path to a genome sequence file (.fa or .fasta)
#' @param bin max distance of cleavage site from PAS
#' @return a table showing the top PASs for each single gene with enough depth in the data and also a bed6 file table for all the PASs
#' @export

PAS_calling <- function(gene_reference, reads,min_reads=5, min_percent=1,cores=1,
                        direct_RNA=FALSE,internal_priming=F,pattern = "post", genome_file=NULL,bin=50){
  if(internal_priming){
    genome <- Rsamtools::FaFile(genome_file)
    Rsamtools::open(genome)
    genome
  }
  gene_info <- gene_reference[,.(gene_id,chrom, strand, gene_name, gene_biotype)]
  reads$read_id <-c(1:nrow(reads))
  samples <- unique(reads$sample)
  reads_dt <- data.table::as.data.table(reads[, !(c("treatment")), with = FALSE])
  data.table::setkey(reads_dt, gene_id)
  
  # Pre-filter the reads
  reads_dt <- reads_dt[gene_id%in%unique(reads_dt[,.N, by = gene_id][N >= min_reads]$gene_id)]  # Filter genes with fewer than min_reads
  genes <- unique(reads_dt$gene_id)
  genes <- as.vector(genes[genes != "."])
  
  
  
  PAS_table <- gene_info[gene_id %in% genes]
  data.table::setkey(PAS_table, gene_id)  # Ensure efficient subsetting
  
  find_match <- function(a, b_vec) {
    distances <- abs(a - b_vec)
    if (any(distances <= bin)) {
      return(b_vec[which.min(distances)])  # Return closest match
    } else {
      return(NA)  # No match within bin
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
      df_3end <- gene_all$chromStart+1
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
        combined_df <- data.table::rbindlist(lapply(split(df, group), function(group_data) {
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
        PAS_sample_counts_raw <- data.table::as.data.table(match_fun(strand,reads_data =  gene_all,call_df))
        tab <- table(PAS_sample_counts_raw$match)
        tab[tab >= min_reads]
        keep <- names(tab[tab >= min_reads])
        PAS_sample_counts <- PAS_sample_counts_raw[match %in% keep]
        if(nrow(PAS_sample_counts)>0){
          PAS_counts <- data.table::as.data.table(table(PAS_sample_counts$match))
          colnames(PAS_counts)<-c("Value","reads")
          PAS_counts[, Value := as.numeric(Value)]
          call_df <- merge(call_df, PAS_counts, by = "Value", all.y = TRUE)
        }else{
          call_df[,"reads"] <- 0
        }
        call_df[, "PAU"] <- round(100*call_df[,"reads"]/nrow(gene_all),3)

        PAS_gene <- gene_info[gene, ][rep(1, nrow(call_df)),]
        PAS_gene$PAS <- call_df$Value
        if (direct_RNA==F&internal_priming){PAS_gene <- Internal_priming(PAS_gene, pattern = pattern,genome = genome)}
        if (nrow(PAS_gene)<nrow(call_df)){
#         False_df <- data.table::as.data.table(call_df)[Value!%in%PAS_gene$PAS]
          call_df <- data.table::as.data.table(call_df)[Value%in%PAS_gene$PAS]
        }
  
        
        if(nrow(call_df)>0){
          add_df<- data.frame(matrix(ncol=6, nrow=(length(call_df$Value))))
          add_df<-data.frame(
            chrom = gene_all$chrom[1], # Numeric column
            start = as.numeric(call_df$Value)-1, # Numeric column
            end = as.numeric(call_df$Value),
            name = paste(gene,gene_info[gene,gene_name]),
            density= as.numeric(call_df$PAU),
            strand = gene_info[gene,"strand"])
          return(add_df)
        }
      }
    }
  }
  
  # Parallel processing
  PAS_output <- pbmcapply::pbmclapply(genes, PAU_fun, mc.cores = cores)
  PAS_bed <- data.table::rbindlist(PAS_output,fill = T)  # Combine results
  
  return(PAS_bed)
}
