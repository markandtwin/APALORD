#' Call PASs at single gene level and calculate the usage of each found PAS
#' 
#' This function extract PolyA sites information of each gene and calculate the usage of each found PAS from the provided RNAseq data.
#' @import dplyr stringr tidyr pbmcapply
#' @param gene_reference information extracted from gtf file
#' @param reads information from RNAseq samples
#' @param cores number of threads used for the computation
#' @param min_reads minium reads count required at single gene level for PAS calling
#' @param min_percent minium percent required for a PAS to be included (0-100)
#' @param direct_RNA whether or not the data is direct RNAseq 
#' @return a table showing the top PASs for each single gene with enough depth in the data and also a bed6 file table for all the PASs
#' @export

PAU_by_sample <- function(gene_reference, reads,min_reads=5, min_percent=1,cores=1,direct_RNA=F){
  gene_info <-gene_reference[1:4]
  reads_df <- reads  %>%group_by(gene_id) %>%filter(n() >= min_reads) %>%ungroup()
  
  genes <- levels(reads_df$gene_id)
  genes <- subset(genes, genes!=".")
  PAS_table <- as.data.frame(subset(gene_info,gene_id%in%genes))
  rownames(PAS_table) <- PAS_table$gene_id
  
  PAU_fun <- function (gene, reads, PAS_table,min_percent=1,min_reads=5,direct_RNA=F){
    gene_all <-subset(reads, gene_id==gene)
    if(direct_RNA){
      gene_all <- subset(gene_all, strand==gene_info[gene,"strand"])
    }
    if(PAS_table[gene,"strand"]=="+"){
      df_3end <-gene_all$chromEnd}
    if(PAS_table[gene,"strand"]=="-"){
      df_3end<-gene_all$chromStart}
    
    if(length(df_3end)>=min_reads){
      freq_table <- table(df_3end)
      density_vals <- as.numeric(prop.table(freq_table))
      peak_indices <- which(density_vals>=min(sort(density_vals, decreasing = T)[1:50],na.rm = T))
      peaks_df <- data.frame(
        Value = as.numeric(names(freq_table))[peak_indices],
        Frequency = as.numeric(freq_table[peak_indices]),
        Density = 100*density_vals[peak_indices]
      )
      df <- peaks_df[order(peaks_df$Value), ]
      group <- cumsum(c(TRUE, diff(df$Value) > 20))
      
      combined_df <- do.call(rbind, lapply(split(df, group), function(group_data) {
        max_density_row <- group_data[which.max(group_data$Density), ]
        max_density_row$Frequency <- sum(group_data$Frequency)
        max_density_row$Density <- sum(group_data$Density)
        return(max_density_row)
      }))
      call_df <- combined_df[which((combined_df$Density>=min_percent)&(combined_df$Frequency>2)),]
      
      # Display the peaks
      call_df<-(call_df[order(call_df$Density, decreasing = T),])
      if(length(call_df$Frequency)!=0){
        call_df$Density <-round(as.numeric(call_df$Density), 4)
        PAS_gene <- PAS_table[gene,1:4][rep(1,nrow(call_df)),]
        PAS_gene$PAS<- call_df$Value
        for (sample in levels(reads$sample)){
          sample_df <- gene_all[which(gene_all$sample==sample),]
          if(PAS_table[gene,"strand"]=="+"){
            PAS_gene[,sample] <-round(sapply(call_df$Value, function(x) {
              mean(abs(sample_df$chromEnd-x)<=20) * 100}),2)
          }
          if(PAS_table[gene,"strand"]=="-"){
            PAS_gene[,sample] <-round(sapply(call_df$Value, function(x) {
              mean(abs(sample_df$chromStart-x)<=20) * 100}),2)
          }
        }
        return(PAS_gene)
      }
    }
  }
  
  PAU_output <-pbmclapply(genes, PAU_fun, reads=reads, PAS_table, min_percent=min_percent, min_reads=min_reads,mc.cores = cores)
  PAU_table <- do.call(rbind, PAU_output)%>% mutate_if(is.character, as.factor)
  return(PAU_table)
}