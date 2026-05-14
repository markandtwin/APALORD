#' Load  samples
#' 
#' This function loads the IsoQuant output file including all the samples for each group and extract necessary information for the analysis.
#' @rdname load_matrix
#' @param  infile  path to the input IsoQuant output files of group 1 samples(Control) 
#' @param  design  experiment design showing the conditions of each sample file
#' @return a table including all the reads from the two groups of samples
#' @export


load_matrix <- function(infile, design = design) {
  
  # Initialize an empty list to store results
  reads_list <- list()
  
  # Define a helper function to read and process files for a given sample
  process_sample <- function(sample, design) {
    base::gc()
    # Read gene reads and remove duplicates
    file1 <- list.files(path = sample, pattern="OUT.read_assignments.tsv",full.names = T)
    sample_gene_reads <- data.table::fread(file1, header = TRUE, sep = "\t", skip = 2)
    colnames(sample_gene_reads)[1] <- "read_id"
    sample_gene_reads <- unique(sample_gene_reads[,c(1,3,5,10)], by = "read_id")
    
    # Read bed file and remove duplicates
    file2 <- list.files(path = sample, pattern="OUT.corrected_reads.bed",full.names = T)
    sample_reads_bed <- data.table::fread(file2, header = TRUE, sep = "\t")
    colnames(sample_reads_bed)[1] <- "chrom"
    sample_reads_bed <- unique(sample_reads_bed[,1:4], by = "name")
    
    # Merge data
    sample_reads <- merge(sample_gene_reads, sample_reads_bed, by.x = "read_id", by.y = "name", all.x = TRUE)
    
    # Add columns for sample and treatment
    sample_reads$sample <- sample_reads$groups
    sample_data <- merge(sample_reads[,c(1:3,5:8)], design, by = "sample")
    return(sample_data)
  }
  
  # Process samples from infile (group1)
  for (sample in infile) {
    message("Collecting data from ", sample)
    sample_data <- process_sample(sample, design)
    reads_list[[length(reads_list) + 1]] <- sample_data
  }
  

  
  # Combine all processed data into a single data.table
  reads_all <- data.table::rbindlist(reads_list)
  reads_all[, (base::setdiff(names(reads_all), c("chromStart", "chromEnd")) ) := lapply(.SD, as.factor), .SDcols = setdiff(names(reads_all), c("chromStart", "chromEnd"))]
  reads_all <-subset(reads_all,gene_id!=".")
#  reads_all$read_id <-seq_len(nrow(reads_all))
  return(reads_all[,c(2:7,1,8)])
}
