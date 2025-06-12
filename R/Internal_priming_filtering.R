#' Call internal priming events
#' 
#' This function extract PolyA sites information from the reference genome to call internal priming originate from  A-rich stretches
#' @import dplyr stringr tidyr GenomicAlignments Biostrings Rsamtools
#' @param PAS_data date table with columns of PAS information including chrom, strand and PAS
#' @param pattern to look for internal priming sequencing surrounding PAS, either "pre", "post" or "both".
#' @param genome genome sequence loaded from a .fasta file.
#' @return the same PAS_data without internal_primed "PAS".
#' @export

Internal_priming <- function(PAS_data, pattern="post", genome = genome){
  PAS_data$A_enrich <- as.numeric(NA)
  for(n in 1:nrow(PAS_data)) {
    gene_PAS <-PAS_data[n,]
    if(gene_PAS$strand=="+"){
      region <- GRanges(gene_PAS$chrom, ranges = IRanges(start = gene_PAS$PAS-9, end = gene_PAS$PAS+10), strand=gene_PAS$strand)
      PAS_sequence <-as.character(getSeq(genome, region))
    } else {
      region <- GRanges(gene_PAS$chrom, ranges = IRanges(start = gene_PAS$PAS-10, end = gene_PAS$PAS+9), strand=gene_PAS$strand)
      PAS_sequence <-as.character(getSeq(genome, region))
    }
    pre_fraction <- sum(strsplit(PAS_sequence, "")[[1]][1:10]=="A")/10
    post_fraction <- sum(strsplit(PAS_sequence, "")[[1]][11:20]=="A")/10
    if(pattern=="post"){
      PAS_data[n,]$A_enrich <- as.numeric(post_fraction)
    }else if(pattern=="pre"){
      PAS_data[n,]$A_enrich <- as.numeric(pre_fraction)
    }else if(pattern=="both"){
      PAS_data[n,]$A_enrich <- max(as.numeric(pre_fraction),as.numeric(post_fraction))  
      }
  }
  
  PAS_data$internal_priming<-ifelse(PAS_data$A_enrich >=0.7, TRUE, FALSE)
  PAS_data <- subset(PAS_data,internal_priming==F)[,!c("A_enrich","internal_priming")]
  return(PAS_data)
}