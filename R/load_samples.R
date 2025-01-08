#' Load two groups of samples
#' 
#' This function loads the output files for each group and extract necessary information for the analysis.
#' @rdname load_samples
#' @import dplyr stringr tidyr
#' @param  infile1  path to the input IsoQuant output files of group 1 samples(Control) 
#' @param  infile2  path to the input IsoQuant output files of group 2 samples(Treated)
#' @param  group1 name or condition of samples imported from infile1  
#' @param  group2 name or condition of samples imported from infile2
#' @return a table including all the reads from the two groups of samples
#' @export

load_samples <- function(infile1, infile2,group1="group1",group2="group2"){
  reads_all <-data.frame(read_id=character(),
                         strand=character(),
                         gene_id=character(),
                         chrom=character(),
                         chromStart=numeric(),
                         chromEnd=numeric(),
                         sample=character(),
                         treatment=character(),stringsAsFactors = TRUE)
  for (sample in infile1){
    sample_gene_reads <- read.csv(paste0(sample, "/OUT.read_assignments.tsv"), header = T, sep = "\t",skip = 2)[,c(1,3,5)]
    sample_gene_reads <- sample_gene_reads[!duplicated(sample_gene_reads$X.read_id), ]
    sample_reads_bed <- read.csv(paste0(sample, "/OUT.corrected_reads.bed"), header = T, sep = "\t",skip = 0)[,1:4]
    sample_reads_bed <- sample_reads_bed[!duplicated(sample_reads_bed$name), ]
    sample_reads <- left_join(sample_gene_reads, sample_reads_bed, by = c("X.read_id" = "name"))
    colnames(sample_reads) <-c("read_id","strand","gene_id","chrom","chromStart","chromEnd")
    sample_reads$sample <- sample
    sample_reads$treatment <-group1 
    sample_reads[c(1:4,7:8)] <- lapply(sample_reads[c(1:4,7:8)], as.factor)
    reads_all <-rbind(reads_all,sample_reads)
  }
  
  for (sample in infile2){
    sample_gene_reads <- read.csv(paste0(sample, "/OUT.read_assignments.tsv"), header = T, sep = "\t",skip = 2)[,c(1,3,5)]
    sample_gene_reads <- sample_gene_reads[!duplicated(sample_gene_reads$X.read_id), ]
    sample_reads_bed <- read.csv(paste0(sample, "/OUT.corrected_reads.bed"), header = T, sep = "\t",skip = 0)[,1:4]
    sample_reads_bed <- sample_reads_bed[!duplicated(sample_reads_bed$name), ]
    sample_reads <- left_join(sample_gene_reads, sample_reads_bed, by = c("X.read_id" = "name"))
    colnames(sample_reads) <-c("read_id","strand","gene_id","chrom","chromStart","chromEnd")
    sample_reads$sample <- sample
    sample_reads$treatment <- group2
    sample_reads[c(1:4,7:8)] <- lapply(sample_reads[c(1:4,7:8)], as.factor)
    reads_all <-rbind(reads_all,sample_reads)
  }
  
  return(reads_all)
}