#' Load a gtf file for reference
#' 
#' This function loads a gtf file and extract necessary information of each gene from it.
#' @import dplyr stringr tidyr readr
#' @param infile path to the input gtf file
#' @param cores number of threads used for the computation
#' @return a table including gene annotation, stop codon and last exon from the gtf file
#' @export

load_gtf <- function(infile,cores=1){
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
  if(length(na.omit(gene_stop_codon_expand$gene_name))==0){
    gene_stop_codon_expand$gene_name <- str_match(gene_stop_codon$V9, "gene_symbol ([^;]+)")[,2]
  }
  gene_stop_codon_expand$strand <- gene_stop_codon_expand$V7
  gene_stop_codon_expand$chrom <- gene_stop_codon_expand$V1
  gene_stop_codon_expand$chromStart <- gene_stop_codon_expand$V4
  gene_stop_codon_expand$chromEnd <- gene_stop_codon_expand$V5
  stop_codon_info <- as.data.frame(gene_stop_codon_expand[,9:15])
  
  
  gene_exon <- subset(gtf_df, V3=="exon")
  gene_exon_expand <- gene_exon %>%
    mutate(V9 = as.character(V9)) %>%
    separate_wider_regex(V9,patterns = c("gene_id ", gene_id=".\\S*","; transcript_id ",
                                         transcript_id=".\\S*","; gene_type "),
                         too_few = "align_start")
  gene_exon_expand$gene_name <- str_match(gene_exon$V9, "gene_name ([^;]+)")[,2]
  if(length(na.omit(gene_exon_expand$gene_name))==0){
    gene_exon_expand$gene_name <- str_match(gene_exon$V9, "gene_symbol ([^;]+)")[,2]
  }
  gene_exon_expand$strand <- gene_exon_expand$V7
  gene_exon_expand$chrom <- gene_exon_expand$V1
  gene_exon_expand$chromStart <- gene_exon_expand$V4
  gene_exon_expand$chromEnd <- gene_exon_expand$V5
  exon_info <- as.data.frame(gene_exon_expand[,9:15])
  
  exon_genes <- unique(exon_info$gene_id)
  extract_fun<- function(gene,stop_codon_info, exon_info, gene_info){
    if(gene_info[gene,"strand"]=="+"){
      df <- subset(exon_info,gene_id==gene)
      gene_exon_info <- df[which.max(df[, "chromEnd"]),-c(2)]
      stop_df<-subset(stop_codon_info,gene_id==gene)
      if(length(stop_df$chromStart)>0){
        gene_exon_info[,c("distal_stop_codon_Start","distal_stop_codon_End")] <- stop_df[which.max(stop_df[,"chromStart"]),c("chromStart","chromEnd")]
      }else{
        gene_exon_info[,c("distal_stop_codon_Start","distal_stop_codon_End")]<-NA
      }
    }
    if(gene_info[gene,"strand"]=="-"){
      df <- subset(exon_info,gene_id==gene)
      gene_exon_info <- df[which.min(df[, "chromStart"]),-c(2)]
      stop_df<-subset(stop_codon_info,gene_id==gene)
      if(length(stop_df$chromStart)>0){
        gene_exon_info[,c("distal_stop_codon_Start","distal_stop_codon_End")] <- stop_df[which.min(stop_df[,"chromStart"]),c("chromStart","chromEnd")]
      }else{
        gene_exon_info[,c("distal_stop_codon_Start","distal_stop_codon_End")]<-NA
      }
    }
    return(gene_exon_info)
  }
  extract_distal <-pbmclapply(exon_genes, extract_fun, stop_codon_info=stop_codon_info,exon_info=exon_info, gene_info=gene_info, mc.cores = cores)
  gene_reference <- do.call(rbind, extract_distal)
  row.names(gene_reference) <- gene_reference$gene_id
  return(gene_reference)
}