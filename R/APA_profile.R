#' Transcriptome wide APA analysis 
#' 
#' This function calls the PASs for each single gene, calculate PAUs for each PAS, and perform the statistics to get the P value and PolyA site usage distance between the selected two groups 
#' @param gene_reference information extracted from gtf file
#' @param reads information from samples
#' @param control which group in the data is used as the control group
#' @param experimental which group in the data is used as the experimental group
#' @param direct_RNA whether the data is direct RNAseq, TRUE or FALSE
#' @param pattern to look for internal priming sequencing surrouding PAS, either "pre", "post" or "both".
#' @param internal_priming whether or not PASs subject to internal priming filtering
#' @param genome_file path to a genome sequence file (.fa or .fasta)
#' @param min_counts minium read counts required for both groups at single gene level to perform the analysis
#' @param min_reads minium read counts required for PAS calling
#' @param min_percent  minium PAU for a PAS to be considered for downstream analysis
#' @param cores number of threads used for the computation
#' @return a table showing the APA change for each single gene with enough depth in the data
#' @export


APA_profile <- function(gene_reference, reads, control, experimental, 
                        min_counts=10, min_reads=5, min_percent=1, 
                        cores=1, direct_RNA=FALSE,internal_priming=F,pattern="post",genome_file=NULL){
  if(internal_priming){
    genome <- Rsamtools::FaFile(genome_file)
    open(genome)
    genome
  }
  
  reads$read_id <-c(1:nrow(reads))
  reads <- data.table::copy(reads)[, !("sample"), with = FALSE]

  # Filter and group control and experimental data
  control_reads <- reads[treatment== control]
  gene_counts <- control_reads[, .N, by = gene_id][N >= min_counts]
  control_df <- control_reads[gene_id %in% gene_counts$gene_id]
  experimental_reads <- reads[treatment== experimental]
  gene_counts <- experimental_reads[, .N, by = gene_id][N >= min_counts]
  experimental_df <- experimental_reads[gene_id %in% gene_counts$gene_id]
  
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
  
  find_match <- function(a, b_vec) {
    distances <- abs(a - b_vec)
    if (any(distances <= 20)) {
      return(b_vec[which.min(distances)])  # Return closest match
    } else {
      return(NA)  # No match within 20
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
  PAS_fun <- function(df_3end, gene_vector,min_percent=1, min_reads=5){
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
      
      call_df <- combined_df[combined_df$Density >= min_percent & combined_df$Frequency >=5, ]
      call_df <- call_df[order(call_df$Value, decreasing = FALSE), ]
      
      False_df <- data.table()
      if (nrow(call_df) > 0) {
        PAS_gene <- gene_vector[rep(1, nrow(call_df)),]
        PAS_gene$PAS <- call_df$Value
        if (direct_RNA==F&internal_priming){PAS_gene <- Internal_priming(PAS_gene,pattern=pattern,genome=genome)}
        if (nrow(PAS_gene)<nrow(call_df)){
          False_df <- as.data.table(call_df)[!Value%in%PAS_gene$PAS]
##          call_df <- as.data.table(call_df)[Value%in%PAS_gene$PAS]
        }

        PAS_gene_table <- as.vector(c(length(call_df$Value),paste(call_df$Value, collapse = ",")))
        return(list(PAS_gene_table,False_df))
      }
    }
  }
  
  # Adjusted KS test for discrete positional data
  ks_tie_adjusted <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    n1 <- length(x)
    n2 <- length(y)
    N <- n1 + n2
    
    vals <- sort(unique(c(x, y)))
    ecdf_x <- ecdf(x)
    ecdf_y <- ecdf(y)
    diffs <- ecdf_x(vals) - ecdf_y(vals)
    
    i <- which.max(abs(diffs))
    D <- abs(diffs[i])
    signedD <- diffs[i]
    max_diff_at <- vals[i]
    
    # frequencies for tie adjustment
    freq <- table(c(x, y)) / N
    tie_var <- sum(freq * (1 - freq))
    n_eff <- n1 * n2 / (n1 + n2)
    
    # adjusted z-score
    z <- D / sqrt(tie_var / n_eff)
    
    # two-sided p-value from normal approximation
    pval <- 2 * (1 - pnorm(z))
    pval <- min(pval, 1)
    
    list(D = D, signedD = signedD, max_diff_at = max_diff_at, p.value = pval)
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
        df_3end <- gene_all$chromEnd
      } else {
        df_3end <- gene_all$chromStart+1
      }
      
      
      gene_vector <- APA_table_name[gene, c("gene_id","chrom","strand")]
      PAS_all <- PAS_fun(df_3end,gene_vector, min_percent, min_reads)
      PAS_info <- PAS_all[[1]]
      if (direct_RNA==F&internal_priming){
        False_PAS <- PAS_all[[2]]$Value
        APA_gene[,"internal_priming"]<-ifelse(length(False_PAS)>0,paste(False_PAS, collapse = ","),as.numeric(NA))
      }
      
      if (is.null(PAS_info)){
        n_PAS <- 0
      } else {
        n_PAS <-length(as.numeric(unlist(strsplit(PAS_info[2], ","))))
      }
      
      if ((n_PAS>0)&(all(table(gene_all$sample)>=min_reads))) {
        APA_gene[,c("number_of_PAS","PAS_coordinates") := as.list(c(PAS_info[1],PAS_info[2]))]

        for (group in c(control,experimental)){
          group_all <- gene_all[treatment==group]
          PAS_gene <- data.table(
            gene = rep(gene, length(as.numeric(unlist(strsplit(PAS_info[2], split = ","))))),
            PAS = as.numeric(unlist(strsplit(PAS_info[2], split = ",")))
          )
          group_corrected <- as.data.table(match_fun(strand,reads_data =  group_all,PAS_gene$PAS))
          PAS_group_counts <- as.data.table(table(group_corrected$match)) 
          if(nrow(PAS_group_counts)>0){
            colnames(PAS_group_counts)<-c("PAS","reads")
            PAS_group_counts[, PAS := as.double(PAS)]
            PAS_gene <- merge(PAS_gene, PAS_group_counts, by = "PAS", all.x = TRUE)
            PAS_gene[is.na(reads), reads := 0L]
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
            control_3end <- group_corrected$corrected_end
            control_PAU <- PAS_gene$PAU
          }else if (group==experimental){
            experimental_3end <- group_corrected$corrected_end
            experimental_PAU <- PAS_gene$PAU
          }
        }
        APA_gene[,"PAU_changes"] <- paste(((experimental_PAU-control_PAU)),collapse = ",")
        test_ks <- ks_tie_adjusted(control_3end, experimental_3end)
        APA_gene$pvalue <-test_ks$p.value
        if (strand == "+") {
          APA_gene$APA_change <- (test_ks$signedD)
        } else {
          APA_gene$APA_change <- (-test_ks$signedD)
        }

        if((strand=="-")&&(max(PAS_gene$PAS) < APA_gene[,"last_exon_chromEnd"])){
          APA_gene[,"APA_type"] <- "Last_exon_tandem_APA"
        } else if((strand=="+")&&(min(PAS_gene$PAS) > APA_gene[,"last_exon_chromStart"])){
          APA_gene[,"APA_type"] <- "Last_exon_tandem_APA"
        } else {
          APA_gene[,"APA_type"] <- "Mixed_APA"
        }
        return(APA_gene)
      }
    }
  }
  
  # Parallelize APA analysis using pbmclapply
  output <- pbmcapply::pbmclapply(genes, APA_fun, mc.cores = cores)
  
  # Combine the results
  output_df <- rbindlist(output, fill = T)
  output_df$number_of_PAS <- as.numeric(output_df$number_of_PAS)
  output_df <- output_df[, !c("last_exon_chromStart", "last_exon_chromEnd"), with = FALSE]
  return(output_df)
}
