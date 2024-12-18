#' Load a gtf file for reference
#' 
#' This function loads a gtf file and extract necessary information of each gene from it.
#' @import dplyr stringr tidyr readr
#' @param infile path to the input gtf file
#' @return a list of two lists including gene annotation and stop codon from the gtf file
#' @export

load_gtf <- function(infile){
  gtf_df <- read.delim(infile, header = F, sep = "\t", comment="#")
  gene_annotation <- subset(gtf_df, V3=="gene")
  gene_annotation_expand <- gene_annotation %>%
    mutate(V9 = as.character(V9)) %>%
    separate_wider_regex(V9,patterns = c("gene_id ", gene_id=".\\S*","; gene_"),
                         too_few = "align_start")
  gene_annotation_expand$gene_name <- str_match(gene_annotation$V9, "gene_name ([^;]+)")[,2]
  if(length(na.omit(gene_annotation_expand$gene_name))==0){
    gene_annotation_expand$gene_name <- str_match(gene_annotation$V9, "gene_symbol ([^;]+)")[,2]
  }
  gene_annotation_expand$gene_biotype<- str_match(gene_annotation$V9, "gene_biotype ([^;]+)")[,2]
  if(length(na.omit(gene_annotation_expand$gene_biotype))==0){
    gene_annotation_expand$gene_biotype<- str_match(gene_annotation$V9, "gene_type ([^;]+)")[,2]
  }
  gene_annotation_expand$strand <- gene_annotation_expand$V7
  
  gene_info <- as.data.frame(gene_annotation_expand[,c("gene_id","gene_name","strand")])
  rownames(gene_info) <- gene_annotation_expand$gene_id
  
  gene_stop_codon <- subset(gtf_df, V3=="stop_codon")
  gene_stop_codon_expand <- gene_stop_codon %>%
    mutate(V9 = as.character(V9)) %>%
    separate_wider_regex(V9,patterns = c("gene_id ", gene_id=".\\S*","; transcript_id ",
                                         transcript_id=".\\S*","; gene_type "),
                         too_few = "align_start")
  gene_stop_codon_expand$gene_name <- str_match(gene_stop_codon$V9, "gene_name ([^;]+)")[,2]
  gene_stop_codon_expand$strand <- gene_stop_codon_expand$V7
  gene_stop_codon_expand$chrom <- gene_stop_codon_expand$V1
  gene_stop_codon_expand$chromStart <- gene_stop_codon_expand$V4
  gene_stop_codon_expand$chromEnd <- gene_stop_codon_expand$V5
  stop_codon_info <- as.data.frame(gene_stop_codon_expand[,9:15])
  
  stop_codon_genes <- unique(stop_codon_info$gene_id)
  for(gene in stop_codon_genes){
    if(gene_info[gene,"strand"]=="+"){
      df <- subset(stop_codon_info,gene_id==gene)
      gene_info[gene, "distal_stop_codon_Start"] <- max(df$chromStart)
      gene_info[gene, "distal_stop_codon_End"] <- max(df$chromEnd)
    }
    if(gene_info[gene,"strand"]=="-"){
      df <- subset(stop_codon_info,gene_id==gene)
      gene_info[gene, "distal_stop_codon_Start"] <- min(df$chromStart)
      gene_info[gene, "distal_stop_codon_End"] <- min(df$chromEnd)
    }
  }
  return(list(gene_info,stop_codon_info))
}