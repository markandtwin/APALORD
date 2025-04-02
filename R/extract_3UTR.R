#' Load a gtf file as reference to extract 3'end annotation
#' 
#' This function loads a gtf file and extract 3'end annotation of each transcript from it.
#' @import  stringr pbmcapply dplyr data.table 
#' @param gtf_file path to the input gtf file
#' @param cores number of threads used for the computation
#' @return a table including transcript 3'end annotation.



extract_3UTR <- function(gtf_file, cores = 1) {
  # Read data using fread for speed
  gtf_df <- fread(gtf_file, header = FALSE, sep = "\t")
  
  
  # Process transcripts
  tx_info <- gtf_df[V3 == "transcript"]
  tx_info[, gene_id := gsub('"','',str_match(V9, "gene_id ([^;]+)")[,2])]
  tx_info[, transcript_id := gsub('"','',str_match(V9, "transcript_id ([^;]+)")[,2])]
  tx_info[, gene_name := gsub('"','',str_match(V9, "gene_name ([^;]+)")[,2])]
  tx_info[is.na(gene_name), gene_name := gsub('"','',str_match(V9, "gene_symbol ([^;]+)")[,2])]
  tx_info[, strand := V7]
  tx_info[, chrom := V1]
  tx_info[, chromStart := V4]
  tx_info[, chromEnd := V5]
  
  tx_info <- unique(tx_info[, .(gene_id,transcript_id,gene_name, chrom, chromStart, chromEnd, strand)])
  setkey(tx_info, transcript_id)
  
  # Extract unique gene IDs from exons
  txs <- unique(tx_info$transcript_id)
  
  # Define the function to extract distal stop codon information
  extract_fun <- function(tx,  tx_info) {
    strand <- tx_info[tx, strand][1]
    single_tx_info <-tx_info[transcript_id == tx][,c(1:4,7)]
    if (strand == "+") {
      single_tx_info[,"PAS"] <- tx_info[transcript_id == tx]$chromEnd
    } else if (strand == "-") {
      single_tx_info[,"PAS"] <- (tx_info[transcript_id == tx]$chromStart+1)
    }
    return(single_tx_info)
  }
  rm(gtf_df)
  gc()
  # Apply function using parallel processing
  extract_3UTR_df <- pbmclapply(txs, extract_fun, 
                               tx_info = tx_info,
                               mc.cores = cores)
  
  # Combine results
  tx_reference <- rbindlist(extract_3UTR_df)
  tx_reference[, c("chrom", "strand") := lapply(.SD, as.factor), .SDcols = c("chrom", "strand")]
  tx_reference[, "PAS" := lapply(.SD, as.numeric), .SDcols = "PAS"]
  setkey(tx_reference, transcript_id)
  
  return(tx_reference)
}
