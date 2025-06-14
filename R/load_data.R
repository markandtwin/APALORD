#' Load a group of sample(s)
#' 
#' This function loads the output files from the same group of sample(s)  
#' @rdname load_data
#' @import dplyr data.table
#' @param  infile  path to the input IsoQuant output files of a group of sample(s)
#' @param  group name or condition of sample(s) imported from infile  
#' @return a table including all the reads from this group of sample(s)
#' @export

load_data <- function(infile,group=group){
  print(paste("Collecting data from ", group, " group"))
  sample <-infile
  file1 <- list.files(path = sample, pattern="OUT.read_assignments.tsv",full.names = T)
  sample_gene_reads <- fread(file1, header = TRUE, sep = "\t", skip = 2)
  colnames(sample_gene_reads)[1] <- "read_id"
  sample_gene_reads <- unique(sample_gene_reads[,c(1,3,5)], by = "read_id")
  
  # Read bed file and remove duplicates
  file2 <- list.files(path = sample, pattern="OUT.corrected_reads.bed",full.names = T)
  sample_reads_bed <- fread(file2, header = TRUE, sep = "\t")
  colnames(sample_reads_bed)[1] <- "chrom"
  sample_reads_bed <- unique(sample_reads_bed[,1:4], by = "name")
  
  # Merge data
  sample_reads <- merge(sample_gene_reads, sample_reads_bed, by.x = "read_id", by.y = "name", all.x = TRUE)
  
  # Add columns for sample and treatment
  sample_reads[, sample := sample]
  sample_reads[, treatment := group]
  reads_all <- sample_reads
  reads_all[, (setdiff(names(reads_all), c("chromStart", "chromEnd")) ) := lapply(.SD, as.factor), .SDcols = setdiff(names(reads_all), c("chromStart", "chromEnd"))]
  reads_all <-subset(reads_all,gene_id!=".")
  return(reads_all)
}