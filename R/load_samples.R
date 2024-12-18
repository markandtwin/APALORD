#' Load samples
#' 
#' This function loads the output files for each sample and extract necessary information for the analysis.
#' @import dplyr stringr tidyr
#' @param  infile1  path to the input IsoQuant output files of sample 1 (Control) 
#' @param  infile2  path to the input IsoQuant output files of sample 2 (Treated)
#' @return a S4 object of two lists including PASs information and APA changes.
#' @export

load_samples <- function(infile1, infile2){
  Undiff_gene_reads <- read.csv(paste0(infile1, "/OUT.read_assignments.tsv"), header = T, sep = "\t",skip = 2)[,c(1,3,5)]
  Undiff_gene_reads <- Undiff_gene_reads[!duplicated(Undiff_gene_reads$X.read_id), ]
  Undiff_reads_bed <- read.csv(paste0(infile1, "/OUT.corrected_reads.bed"), header = T, sep = "\t",skip = 0)[,1:4]
  Undiff_reads_bed <- Undiff_reads_bed[!duplicated(Undiff_reads_bed$name), ]
  Undiff_reads <- left_join(Undiff_gene_reads, Undiff_reads_bed, by = c("X.read_id" = "name"))
  colnames(Undiff_reads) <-c("read_id","strand","gene_id","chrom","chromStart","chromEnd")
  Undiff_reads[1:4] <- lapply(Undiff_reads[1:4], as.factor)
  Undiff_reads$treatment <- "Undiff"
  
  Diff_gene_reads <- read.csv(paste0(infile2, "/OUT.read_assignments.tsv"), header = T, sep = "\t",skip = 2)[,c(1,3,5)]
  Diff_gene_reads <- Diff_gene_reads[!duplicated(Diff_gene_reads$X.read_id), ]
  Diff_reads_bed <- read.csv(paste0(infile2, "/OUT.corrected_reads.bed"), header = T, sep = "\t",skip = 0)[,1:4]
  Diff_reads_bed <- Diff_reads_bed[!duplicated(Diff_reads_bed$name), ]
  Diff_reads <- left_join(Diff_gene_reads, Diff_reads_bed, by = c("X.read_id" = "name"))
  colnames(Diff_reads) <-c("read_id","strand","gene_id","chrom","chromStart","chromEnd")
  Diff_reads[1:4] <- lapply(Diff_reads[1:4], as.factor)
  Diff_reads$treatment <- "Diff"
  
  reads_all <-rbind(Undiff_reads, Diff_reads)
  reads_all$treatment <-as.factor(reads_all$treatment)
  return(reads_all)
}