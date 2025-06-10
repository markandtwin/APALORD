#' Load two groups of samples from bam files directly
#' 
#' This function loads the reads information from bam files for each group and extract necessary information for the analysis.
#' @rdname load_from_bam
#' @import  dplyr data.table GenomicFeatures Rsamtools GenomicAlignments
#' @param  infile1  path to the input IsoQuant output files of group 1 samples(Control) 
#' @param  infile2  path to the input IsoQuant output files of group 2 samples(Treated)
#' @param  group1 name or condition of samples imported from infile1  
#' @param  group2 name or condition of samples imported from infile2
#' @param  gtf annotation file to assign the read to a gene in the bam file.
#' @return a table including all the reads from the two groups of samples
#' @export
load_from_bam <- function(infile1, infile2, group1="group1", group2="group2",gtf_file=gtf_file){
  read_data_all <- data.table()
  txdb <- makeTxDbFromGFF(gtf_file)
  genes <- genes(txdb) 
  for (sample in infile1){
    reads <- readGAlignments(sample, use.names = TRUE)
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
      gene_id = NA,
      seqname = seqnames,
      chromStart = five_prime,
      chromEnd = three_prime,
      sample = sample,
      treatment = group1
    )
    # Assign to gene
    overlaps <- findOverlaps(reads, genes)
    read_data$gene_id[queryHits(overlaps)] <- names(genes)[subjectHits(overlaps)]
    read_data_all <- rbind(read_data_all, read_data)
  }
  
  for (sample in infile2){
    reads <- readGAlignments(sample, use.names = TRUE)
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
      gene_id = NA,
      seqname = seqnames,
      chromStart = five_prime,
      chromEnd = three_prime,
      sample = sample,
      treatment = group2
    )
    # Assign to gene
    overlaps <- findOverlaps(reads, genes)
    read_data$gene_id[queryHits(overlaps)] <- names(genes)[subjectHits(overlaps)]
    read_data_all <- rbind(read_data_all, read_data)
  }
  read_data_all <-unique(read_data_all, by = "read_id")
  read_data_all <-subset(read_data_all,!is.na(gene_id))
  return(read_data_all)
}