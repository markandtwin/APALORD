#' Load a gtf file for reference
#' 
#' This function loads a gtf file and extract necessary information of each gene from it.
#' @import  stringr pbmcapply dplyr data.table
#' @param infile path to the input gtf file
#' @param cores number of threads used for the computation
#' @return a table including gene annotation, stop codon and last exon from the gtf file
#' @export


load_gtf <- function(infile, cores = 1) {
  # Read data using fread for speed
  gtf_df <- fread(infile, header = FALSE, sep = "\t")
  
  # Filter and process 'gene' annotations
  gene_annotation <- gtf_df[V3 == "gene"]
  
  # Efficiently extract gene information
  gene_annotation[, c("gene_id", "gene_name", "gene_biotype") := 
                    list(
                      str_match(V9, "gene_id ([^;]+)")[,2],
                      str_match(V9, "gene_name ([^;]+)")[,2],
                      str_match(V9, "gene_biotype ([^;]+)")[,2]
                    )]
  
  # If gene_name is missing, use gene_symbol
  gene_annotation[is.na(gene_name), gene_name := str_match(V9, "gene_symbol ([^;]+)")[,2]]
  gene_annotation[, strand := V7]
  
  # If gene_name is missing, use gene_symbol
  gene_annotation[is.na(gene_biotype), gene_biotype := str_match(V9, "gene_type ([^;]+)")[,2]]
  gene_annotation[, strand := V7]
  
  # Create gene_info table
  gene_info <- unique(gene_annotation[, .(gene_id, gene_name, strand)])
  setkey(gene_info, gene_id)
  
  # Process stop codons
  gene_stop_codon <- gtf_df[V3 == "stop_codon"]
  gene_stop_codon[, c("gene_id", "transcript_id") := 
                    list(
                      str_match(V9, "gene_id ([^;]+)")[,2],
                      str_match(V9, "transcript_id ([^;]+)")[,2]
                    )]
  gene_stop_codon[, gene_name := str_match(V9, "gene_name ([^;]+)")[,2]]
  gene_stop_codon[is.na(gene_name), gene_name := str_match(V9, "gene_symbol ([^;]+)")[,2]]
  gene_stop_codon[, strand := V7]
  gene_stop_codon[, chrom := V1]
  gene_stop_codon[, chromStart := V4]
  gene_stop_codon[, chromEnd := V5]
  
  # Process exons
  gene_exon <- gtf_df[V3 == "exon"]
  gene_exon[, c("gene_id", "transcript_id") := 
              list(
                str_match(V9, "gene_id ([^;]+)")[,2],
                str_match(V9, "transcript_id ([^;]+)")[,2]
              )]
  gene_exon[, gene_name := str_match(V9, "gene_name ([^;]+)")[,2]]
  gene_exon[is.na(gene_name), gene_name := str_match(V9, "gene_symbol ([^;]+)")[,2]]
  gene_exon[, strand := V7]
  gene_exon[, chrom := V1]
  gene_exon[, chromStart := V4]
  gene_exon[, chromEnd := V5]
  
  exon_info <- unique(gene_exon[, .(gene_id, chrom, chromStart, chromEnd, strand, gene_name)])
  
  # Extract unique gene IDs from exons
  exon_genes <- unique(exon_info$gene_id)
  
  # Define the function to extract distal stop codon information
  extract_fun <- function(gene, stop_codon_info, exon_info, gene_info) {
    strand <- gene_info[gene, strand]
    
    if (strand == "+") {
      df <- exon_info[gene_id == gene]
      gene_exon_info <- df[which.max(chromEnd),]
      stop_df <- stop_codon_info[gene_id == gene]
      if (nrow(stop_df) > 0) {
        gene_exon_info[, c("distal_stop_codon_Start", "distal_stop_codon_End") := 
                         stop_df[which.max(chromStart), .(chromStart, chromEnd)]]
      } else {
        gene_exon_info[, c("distal_stop_codon_Start", "distal_stop_codon_End") := NA]
      }
    } else if (strand == "-") {
      df <- exon_info[gene_id == gene]
      gene_exon_info <- df[which.min(chromStart),]
      stop_df <- stop_codon_info[gene_id == gene]
      if (nrow(stop_df) > 0) {
        gene_exon_info[, c("distal_stop_codon_Start", "distal_stop_codon_End") := 
                         stop_df[which.min(chromStart), .(chromStart, chromEnd)]]
      } else {
        gene_exon_info[, c("distal_stop_codon_Start", "distal_stop_codon_End") := NA]
      }
    }
    
    return(gene_exon_info)
  }
  
  # Apply function using parallel processing
  extract_distal <- pbmclapply(exon_genes, extract_fun, 
                               stop_codon_info = gene_stop_codon, 
                               exon_info = exon_info, 
                               gene_info = gene_info, 
                               mc.cores = cores)
  
  # Combine results
  gene_reference <- rbindlist(extract_distal)%>% mutate_if(is.character, as.factor)
  gene_reference$gene_id <- gsub('"', '', gene_reference$gene_id)
  gene_reference$gene_name <- gsub('"', '', gene_reference$gene_name)
  setDT(gene_reference)
  setkey(gene_reference, gene_id)

  return(gene_reference)
}
