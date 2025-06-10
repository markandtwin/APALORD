#' Load two groups of samples
#' 
#' This function loads the IsoQuant output files for each group and extract necessary information for the analysis.
#' @rdname load_samples
#' @import  dplyr data.table 
#' @param  infile1  path to the input IsoQuant output files of group 1 samples(Control) 
#' @param  infile2  path to the input IsoQuant output files of group 2 samples(Treated)
#' @param  group1 name or condition of samples imported from infile1  
#' @param  group2 name or condition of samples imported from infile2
#' @return a table including all the reads from the two groups of samples
#' @export


load_samples <- function(infile1, infile2, group1="group1", group2="group2") {
  
  # Initialize an empty list to store results
  reads_list <- list()
  
  # Define a helper function to read and process files for a given sample
  process_sample <- function(sample, group) {
    gc()
    # Read gene reads and remove duplicates
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
    return(sample_reads)
  }
  
  # Process samples from infile1 (group1)
  for (sample in infile1) {
    print(paste0("Collecting data from ", group1, " group"))
    sample_data <- process_sample(sample, group1)
    reads_list[[length(reads_list) + 1]] <- sample_data
  }
  
  # Process samples from infile2 (group2)
  for (sample in infile2) {
    print(paste("Collecting data from ", group2, " group"))
    sample_data <- process_sample(sample, group2)
    reads_list[[length(reads_list) + 1]] <- sample_data
  }
  
  # Combine all processed data into a single data.table
  reads_all <- rbindlist(reads_list)
  reads_all[, (setdiff(names(reads_all), c("chromStart", "chromEnd")) ) := lapply(.SD, as.factor), .SDcols = setdiff(names(reads_all), c("chromStart", "chromEnd"))]
  reads_all <-subset(reads_all,gene_id!=".")
  return(reads_all)
}
