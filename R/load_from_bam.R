#' Load two groups of samples from bam files directly
#' 
#' This function loads the reads information from bam files for each group and extract necessary information for the analysis.
#' @rdname load_from_bam
#' @import  dplyr data.table GenomicFeatures Rsamtools GenomicAlignments purrr dplyr bambu pbmcapply data.table
#' @param  infile1  path to the input IsoQuant output files of group 1 samples(Control) 
#' @param  infile2  path to the input IsoQuant output files of group 2 samples(Treated)
#' @param  group1 name or condition of samples imported from infile1  
#' @param  group2 name or condition of samples imported from infile2
#' @param  gtf_file path to annotation file to assign the read to a gene in the bam file
#' @param  bambu whether or not to turn off bambu when assign read to gene. It's highly recommended to keep it on to make sure the assignment is accurate although it takes time to run it
#' @param  genome_file path to genome sequence file (.fa or .fasta) that was used to align the reads in sample bam files
#' @param  cores number of threads used to process the data
#' @return a table including all the reads from the two groups of samples
#' @export
load_from_bam <- function(infile1, infile2, group1="group1", group2="group2",gtf_file,bambu=T,genome_file=NULL,cores=1){
  read_data_all <- data.table()
  txdb <- makeTxDbFromGFF(gtf_file)
  genes <- genes(txdb)
  
  if (bambu){
    se.multiSample<- bambu(reads = c(infile1,infile2), annotations = annotations, genome = genome.file, trackReads = T, ncore = cores)
    genes <- intersect(unique(rowData(se.multiSample)$GENEID),genes$gene_id)
  }
  
  assign_gene <- function(gene){
    txs <- which(rowData(se.multiSample)$GENEID == gene)
    meta_tbl <- se.multiSample@metadata$readToTranscriptMaps[[i]]
    
    filtered_tbl <- meta_tbl %>%
      filter(
        map_lgl(equalMatches, ~ !is.null(.x) && any(.x %in% txs)) |
          map_lgl(compatibleMatches, ~ !is.null(.x) && any(.x %in% txs))
      )
    
    if (nrow(filtered_tbl) > 0) {
      data.table(read_id = filtered_tbl$readId, gene_id = gene)
    } else {
      NULL
    }
  }
  
  
  
  for (i in 1:length(infile1)){
    reads <- readGAlignments(infile1[i], use.names = TRUE)
    # Read ID
    read_ids <- names(reads)
    
    
    # Chromosome / seqname
    seqnames <- as.character(seqnames(reads))
    
    # Strand
    strand_info <- as.character(strand(reads))
    
    # Compute 5' and 3' end based on strand
    five_prime <- ifelse(strand_info == "+", end(reads), start(reads))
    three_prime <- ifelse(strand_info == "+", end(reads), start(reads))
    
    # Combine into a data frame
    read_data <- data.table(
      read_id = read_ids,
      strand = strand_info,
      gene_id = "NA",
      seqname = seqnames,
      chromStart = five_prime,
      chromEnd = three_prime,
      sample = infile1[i],
      treatment = group1
    )
    read_data <-unique(read_data, by = "read_id")
    # Assign to gene
    
    if(bambu){
      results_list <- pbmclapply(genes, assign_gene, mc.cores = cores)
      
      # Combine results
      read_gene_map <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
      read_gene_map <- unique(read_gene_map, by = "read_id")
      read_data[read_id %in% read_gene_map$read_id, gene_id := read_gene_map$gene_id[match(read_id, read_gene_map$read_id)]]
    }else{
      overlaps <- findOverlaps(reads, genes)
      read_data$gene_id[queryHits(overlaps)] <- names(genes)[subjectHits(overlaps)]
    }
    read_data_all <- rbind(read_data_all, read_data)
  }
  
  for (i in 1:length(infile2)){
    reads <- readGAlignments(infile2[i], use.names = TRUE)
    # Read ID
    read_ids <- names(reads)
    
    # Chromosome / seqname
    seqnames <- as.character(seqnames(reads))
    
    # Strand
    strand_info <- as.character(strand(reads))
    
    # Compute 5' and 3' end based on strand
    five_prime <- ifelse(strand_info == "+", end(reads), start(reads))
    three_prime <- ifelse(strand_info == "+", end(reads), start(reads))
    
    # Combine into a data frame
    read_data <- data.table(
      read_id = read_ids,
      strand = strand_info,
      gene_id = "NA",
      seqname = seqnames,
      chromStart = five_prime,
      chromEnd = three_prime,
      sample = infile2[i],
      treatment = group2
    )
    read_data <-unique(read_data, by = "read_id")
    # Assign to gene
    i <- i+length(infile1)
    if(bambu){
      results_list <- pbmclapply(genes, assign_gene, mc.cores = cores)
      
      # Combine results
      read_gene_map <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
      read_gene_map <- unique(read_gene_map, by = "read_id")
      read_data[read_id %in% read_gene_map$read_id, gene_id := read_gene_map$gene_id[match(read_id, read_gene_map$read_id)]]
    }else{
      overlaps <- findOverlaps(reads, genes)
      read_data$gene_id[queryHits(overlaps)] <- names(genes)[subjectHits(overlaps)]
    }
    read_data_all <- rbind(read_data_all, read_data)
  }
    
  read_data_all <- subset(read_data_all,gene_id!="NA")
  return(read_data_all)
}