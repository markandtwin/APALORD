#' Transcriptome wide APA analysis 
#' 
#' This function calls the PASs for each single gene, calculate PAUs for each PAS, and perform the statistics to get the P value and PolyA site usage distance between the selected two groups 
#' @param gene_reference information extracted from gtf file
#' @param reads information from samples
#' @param control which group in the data is used as the control group
#' @param experimental which group in the data is used as the experimental group
#' @param direct_RNA whether the data is direct RNA-Seq, TRUE or FALSE
#' @param pattern to look for internal priming sequencing surrounding PAS, either "pre", "post" or "both".
#' @param internal_priming whether or not PASs subject to internal priming filtering
#' @param genome_file path to a genome sequence file (.fa or .fasta)
#' @param min_counts minimum read counts required for both groups at single gene level to perform the analysis
#' @param min_reads minimum read counts required for PAS calling per sample
#' @param min_percent  minimum PAU for a PAS to be considered for downstream analysis
#' @param cores number of threads used for the computation
#' @param bin max distance of cleavage site from PAS
#' @return a table showing the APA change for each single gene with enough depth in the data
#' @export


APA_profile_compare <- function(gene_reference, reads, control, experimental, 
                        min_counts=10, min_reads=5, min_percent=1, 
                        cores=1, direct_RNA=FALSE,internal_priming=F,pattern="post",genome_file=NULL,bin=50){
  if(internal_priming){
    genome <- Rsamtools::FaFile(genome_file)
    open(genome)
    genome
  }
  
  reads$read_id <-c(1:nrow(reads))
  
# reads <- data.table::copy(reads)[, !("sample"), with = FALSE]

  # Filter and group control and experimental data
  data.table::setnames(reads, "sample", "sample_col")
  control_reads <- reads[treatment== control]
  gene_counts <- control_reads[, .N, by = .(gene_id, sample_col)][N >= min_counts]
  control_df <- control_reads[gene_id %in% gene_counts$gene_id]
  experimental_reads <- reads[treatment== experimental]
  gene_counts <- experimental_reads[, .N, by = .(gene_id, sample_col)][N >= min_counts]
  experimental_df <- experimental_reads[gene_id %in% gene_counts$gene_id]
  
  # Genes of interest
  genes <- intersect(control_df$gene_id, experimental_df$gene_id)
  genes <- genes[genes != "."]
  reads <- data.table::copy(reads)[, !("sample_col"), with = FALSE]
  # Pre-allocate the APA_table
  APA_table <- data.table::data.table(gene_id = genes, 
                          short_gene_id = gsub("\\.\\d+$", "", genes), 
                          APA_change = NA, pvalue = NA, APA_type = "No APA")
  
  # Merge APA_table with gene_info
  gene_info <-gene_reference
  APA_table_name <- merge(APA_table, gene_info, by = "gene_id", all.x = TRUE)
  data.table::setkey(APA_table_name, gene_id)  # Set key for faster subsetting
  
  find_match <- function(a, b_vec) {
    distances <- abs(a - b_vec)
    if (any(distances <= bin)) {
      return(b_vec[which.min(distances)])  # Return closest match
    } else {
      return(NA)  # No match within bin
    }
  }
  
  match_fun <- function(strand,reads_data,PASs_gene){
    if(strand=="+"){
      reads_data$end <- reads_data$chromEnd
      matches <- sapply(reads_data$end, find_match, b_vec = PASs_gene)
      reads_data$match <-matches
    }
    if(strand=="-"){
      reads_data$end <- reads_data$chromStart+1
      matches <- sapply(reads_data$end, find_match, b_vec = PASs_gene)
      reads_data$match <-matches
    }
    reads_data$corrected_end <- ifelse(
      is.na(reads_data$match),
      reads_data$end,
      reads_data$match
    )
    return(reads_data)
  }
  
  # Define PAS_fun with data.table operations for speed
  PAS_fun <- function(gene_all, gene_vector,min_percent=1, min_reads=5){
    if (gene_vector$strand == "+") {
      df_3end <- gene_all$chromEnd
    } else {
      df_3end <- gene_all$chromStart+1
    }
    if (length(df_3end) >= min_reads) {
      freq_table <- table(df_3end)
      density_vals <- as.numeric(prop.table(freq_table))
      peak_indices <- which(density_vals >= min(sort(density_vals, decreasing = TRUE)[1:50], na.rm = TRUE))
      
      peaks_df <- data.table::data.table(Value = as.numeric(names(freq_table))[peak_indices], 
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
      call_df <- call_df[order(call_df$Value, decreasing = FALSE), ]
      
      if (nrow(call_df) > 0) {
        PAS_sample_counts_raw <- data.table::as.data.table(match_fun(gene_vector$strand,reads_data =  gene_all,call_df$Value))
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
        
        PAS_gene <- gene_vector[rep(1, nrow(call_df)),]
        PAS_gene$PAS <- call_df$Value
        if (direct_RNA==F&internal_priming){PAS_gene <- Internal_priming(PAS_gene, pattern = pattern,genome = genome)}
        if (nrow(PAS_gene)<nrow(call_df)){
          False_df <- data.table::as.data.table(call_df)[!Value%in%PAS_gene$PAS]
          call_df <- data.table::as.data.table(call_df)[Value%in%PAS_gene$PAS]
        }else{
          False_df <- NA
        }
        return(list(call_df,False_df,PAS_sample_counts_raw))
      }
    }
  }
  
  # Adjusted KS test for discrete positional data
  ks_tie_adjusted <- function(x, y, B = 5000, seed = NULL) {
    x <- as.numeric(x); y <- as.numeric(y)
    n1 <- length(x); n2 <- length(y)
    pooled <- c(x, y)
    vals <- sort(unique(pooled))
    
    diffs <- ecdf(x)(vals) - ecdf(y)(vals)
    i <- which.max(abs(diffs))
    
    D_obs <- abs(diffs[i])
    signedD <- diffs[i]
    max_diff_at <- vals[i]
    
    if (!is.null(seed)) set.seed(seed)
    D_perm <- replicate(B, {
      idx <- sample.int(n1 + n2, n1, replace = FALSE)
      x0 <- pooled[idx]; y0 <- pooled[-idx]
      d0 <- ecdf(x0)(vals) - ecdf(y0)(vals)
      max(abs(d0))
    })
    
    pval <- (1 + sum(D_perm >= D_obs)) / (B + 1)
    
    list(D = D_obs, signedD = signedD, max_diff_at = max_diff_at, p.value = pval,
         n1 = n1, n2 = n2)
  }
  
  # Define APA_fun optimized with data.table
  APA_fun <- function(gene) {
    gene_all <- reads[gene_id == gene]
    APA_gene <- APA_table_name[gene_id == gene]
    
    if (direct_RNA == TRUE) {
      gene_all <- gene_all[strand == (APA_gene$strand)]
    }
    
    if (all(table(factor(gene_all$treatment,levels=c(control,experimental)))>min_counts)) {
      strand <- gene_info[gene_id==gene]$strand

      gene_vector <- gene_info[gene_id==gene]
      PAS_all <- PAS_fun(gene_all,gene_vector, min_percent, min_reads)
      PAS_info <- PAS_all[[1]]
      if (direct_RNA==F&internal_priming){
        False_PAS <- PAS_all[[2]]$Value
        APA_gene[,"internal_priming"]<-ifelse(length(False_PAS)>0,paste(False_PAS, collapse = ","),as.numeric(NA))
      }
      
      if (is.null(PAS_info)){
        n_PAS <- 0
      } else {
        n_PAS <-nrow(PAS_info)
      }
      
      if (n_PAS>0) {
        PAS_sample_counts <- PAS_all[[3]]
        all_matches <- PAS_all[[1]]$Value
        APA_gene[, PAS_coordinates := paste(all_matches, collapse = ",")]
        APA_gene[, number_of_PAS := length(all_matches)]
        
        for (group in c(control,experimental)){
          group_all <- gene_all[treatment==group]
          PAS_gene <- data.table::data.table(
            gene = rep(gene, length(all_matches)),
            PAS = as.numeric(all_matches)
          )
          PAS_group_counts <- data.table::as.data.table(table(factor(PAS_sample_counts[treatment==group]$match, levels=all_matches))) 
          if(nrow(PAS_group_counts)>0){
            colnames(PAS_group_counts)<-c("PAS","reads")
            PAS_group_counts[, PAS := as.double(PAS)]
            PAS_gene <- merge(PAS_gene, PAS_group_counts, by = "PAS", all.x = TRUE)
          }else{
            PAS_gene[,"reads"] <- 0
          }
          PAS_gene[,"PAU"] <- round(100*PAS_gene[,"reads"]/nrow(group_all),3)
          APA_gene[, c(
            paste0("reads_", group),
            paste0("PAUs_", group)
          ) := as.list(c(
            paste(PAS_gene$reads, collapse = ","),
            paste(PAS_gene$PAU,   collapse = ",")
          ))]
          if(group==control){
            control_3end <- PAS_sample_counts[treatment==group]$corrected_end
            control_PAU <- PAS_gene$PAU
          }else if (group==experimental){
            experimental_3end <- PAS_sample_counts[treatment==group]$corrected_end
            experimental_PAU <- PAS_gene$PAU
          }
        }
        APA_gene[,"PAU_changes"] <- paste(round((experimental_PAU-control_PAU),3),collapse = ",")
        test_ks <- ks_tie_adjusted(control_3end, experimental_3end)
        APA_gene$pvalue <-test_ks$p.value
        if (strand == "+") {
          APA_gene$APA_change <- (test_ks$signedD)
        } else {
          APA_gene$APA_change <- (-test_ks$signedD)
        }
        if(APA_gene$number_of_PAS>1){
          if((strand=="-")&&(max(PAS_gene$PAS) < APA_gene[,"last_exon_chromEnd"])){
            APA_gene[,"APA_type"] <- "Last_exon_tandem_APA"
          } else if((strand=="+")&&(min(PAS_gene$PAS) > APA_gene[,"last_exon_chromStart"])){
            APA_gene[,"APA_type"] <- "Last_exon_tandem_APA"
          } else {
            APA_gene[,"APA_type"] <- "Mixed_APA"
          }
        }else{
          APA_gene[,"APA_type"] <- "Single PAS"
        }
        return(APA_gene)
      }
    }
  }
  
  # Parallelize APA analysis using pbmclapply
  output <- pbmcapply::pbmclapply(genes, APA_fun, mc.cores = cores)
  
  # Combine the results
  output_df <- data.table::rbindlist(output, fill = T)
  output_df$number_of_PAS <- as.numeric(output_df$number_of_PAS)
  output_df <- output_df[, !c("last_exon_chromStart", "last_exon_chromEnd"), with = FALSE]
  return(output_df)
  gc()
}
