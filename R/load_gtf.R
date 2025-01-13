#' Load a gtf file for reference
#' 
#' This function loads a gtf file and extract necessary information of each gene from it.
#' @import  stringr pbmcapply dplyr data.table 
#' @param gtf_file path to the input gtf file
#' @param cores number of threads used for the computation
#' @return a table including gene annotation, stop codon and last exon from the gtf file
#' @export


load_gtf <- function(gtf_file, cores = 1) {
  # Read data using fread for speed
  gtf_df <- fread(gtf_file, header = FALSE, sep = "\t")
  
  # Filter and process 'gene' annotations
  gene_info <- gtf_df[V3 == "gene"]
  
  # Efficiently extract gene information
  gene_info[, c("gene_id", "gene_name", "gene_biotype") := 
                    list(
                      str_match(V9, "gene_id ([^;]+)")[,2],
                      str_match(V9, "gene_name ([^;]+)")[,2],
                      str_match(V9, "gene_biotype ([^;]+)")[,2]
                    )]
  
  # If gene_name is missing, use gene_symbol
  gene_info[is.na(gene_name), gene_name := str_match(V9, "gene_symbol ([^;]+)")[,2]]
  gene_info[, strand := V7]
  
  # If gene_name is missing, use gene_symbol
  gene_info[is.na(gene_biotype), gene_biotype := str_match(V9, "gene_type ([^;]+)")[,2]]
  gene_info[, strand := V7]
  
  # Create gene_info table
  gene_info <- unique(gene_info[, .(gene_id, gene_name, gene_biotype)])
  setkey(gene_info, gene_id)
  
  
  
  # Process exons
  exon_info <- gtf_df[V3 == "exon"]
  exon_info[, c("gene_id", "transcript_id") := 
              list(
                str_match(V9, "gene_id ([^;]+)")[,2],
                str_match(V9, "transcript_id ([^;]+)")[,2]
              )]
  exon_info[, gene_name := str_match(V9, "gene_name ([^;]+)")[,2]]
  exon_info[is.na(gene_name), gene_name := str_match(V9, "gene_symbol ([^;]+)")[,2]]
  exon_info[, strand := V7]
  exon_info[, chrom := V1]
  exon_info[, chromStart := V4]
  exon_info[, chromEnd := V5]
  
  exon_info <- unique(exon_info[, .(gene_id, chrom, chromStart, chromEnd, strand)])
  setkey(exon_info, gene_id)
  
  # Extract unique gene IDs from exons
  exon_genes <- unique(exon_info$gene_id)
  
  # Define the function to extract distal stop codon information
  extract_fun <- function(gene,  exon_info) {
    strand <- exon_info[gene, strand][1]
    
    if (strand == "+") {
      gene_exon_info <- exon_info[gene_id == gene][which.max(chromEnd), ]
    } else if (strand == "-") {
      gene_exon_info <- exon_info[gene_id == gene][which.min(chromStart), ]
    }
    return(gene_exon_info)
  }
  rm(gtf_df)
  gc()
  # Apply function using parallel processing
  extract_distal <- pbmclapply(exon_genes, extract_fun, 
                               exon_info = exon_info,
                               mc.cores = cores)
  
  # Combine results
  exon_reference <- rbindlist(extract_distal)%>% mutate_if(is.character, as.factor)
  gene_reference <- merge(exon_reference, gene_info, by = "gene_id", all.x =  TRUE)
  gene_reference$gene_id <- gsub('"', '', gene_reference$gene_id)
  gene_reference$gene_name <- gsub('"', '', gene_reference$gene_name)
  gene_reference$gene_biotype<- gsub('"', '', gene_reference$gene_biotype)
  setnames(gene_reference, c("chromStart", "chromEnd"), c("last_exon_chromStart", "last_exon_chromEnd"))
  setkey(gene_reference, gene_id)

  return(gene_reference)
}
