#' Load a group of sample(s)
#' 
#' This function loads the output files from the same group of sample(s)  
#' @rdname load_data
#' @import dplyr stringr tidyr
#' @param  infile  path to the input IsoQuant output files of a group of sample(s)
#' @param  group name or condition of sample(s) imported from infile  
#' @return a table including all the reads from this group of sample(s)
#' @export

load_data <- function(infile,group=group){
  reads_all <-data.frame(read_id=character(),
                         strand=character(),
                         gene_id=character(),
                         chrom=character(),
                         chromStart=numeric(),
                         chromEnd=numeric(),
                         sample=character(),
                         treatment=character(),stringsAsFactors = TRUE)
  for (sample in infile){
    sample_gene_reads <- read.csv(paste0(sample, "/OUT.read_assignments.tsv"), header = T, sep = "\t",skip = 2)[,c(1,3,5)]
    sample_gene_reads <- sample_gene_reads[!duplicated(sample_gene_reads$X.read_id), ]
    sample_reads_bed <- read.csv(paste0(sample, "/OUT.corrected_reads.bed"), header = T, sep = "\t",skip = 0)[,1:4]
    sample_reads_bed <- sample_reads_bed[!duplicated(sample_reads_bed$name), ]
    sample_reads <- left_join(sample_gene_reads, sample_reads_bed, by = c("X.read_id" = "name"))
    colnames(sample_reads) <-c("read_id","strand","gene_id","chrom","chromStart","chromEnd")
    sample_reads$sample <- sample
    sample_reads$treatment <- group
    sample_reads[c(1:4,7:8)] <- lapply(sample_reads[c(1:4,7:8)], as.factor)
    reads_all <-rbind(reads_all,sample_reads)
  }
  return(reads_all)
}