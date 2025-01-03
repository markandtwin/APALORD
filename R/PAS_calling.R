#' Call PASs at single gene level
#' 
#' This function extract PolyA sites information of each gene from the provided RNAseq data.
#' @import dplyr stringr tidyr pbmcapply
#' @param gene_reference information extracted from gtf file
#' @param reads information from RNAseq samples
#' @param cores number of threads used for the computation
#' @param counts minium reads count required at single gene level for PAS calling
#' @param min minium fraction required for a PAS to be included
#' @param direct_RNA whether or not the data is direct RNAseq 
#' @return a table showing the top PASs for each single gene with enough depth in the data and also a bed6 file table for all the PASs
#' @export

PAS_calling <- function(gene_reference, reads,counts=5, min=0.01,cores=1,direct_RNA=TRUE){
  gene_info <-gene_reference
  reads_df <- reads  %>%group_by(gene_id) %>%filter(n() >= counts) %>%ungroup()
  
  genes <- levels(reads_df$gene_id)
  genes <- subset(genes, genes!=".")
  PAS_table <- as.data.frame(subset(gene_info,gene_id%in%genes))
  rownames(PAS_table) <- PAS_table$gene_id
  
  PAS_fun <- function (gene, reads, PAS_table,min=0.01,counts=5,direct_RNA=TRUE){
    gene_all <-subset(reads, gene_id==gene)
    if(direct_RNA){
      gene_all <- subset(gene_all, strand==gene_info[gene,"strand"])
    }
    if(PAS_table[gene,"strand"]=="+"){
      df_3end <-gene_all$chromEnd}
    if(PAS_table[gene,"strand"]=="-"){
      df_3end<-gene_all$chromStart}
    
    if(length(df_3end)>=counts){
      freq_table <- table(df_3end)
      density_vals <- as.numeric(prop.table(freq_table))
      peak_indices <- which(density_vals>=min(sort(density_vals, decreasing = T)[1:50],na.rm = T))
      peaks_df <- data.frame(
        Value = as.numeric(names(freq_table))[peak_indices],
        Frequency = as.numeric(freq_table[peak_indices]),
        Density = density_vals[peak_indices]
      )
      df <- peaks_df[order(peaks_df$Value), ]
      group <- cumsum(c(TRUE, diff(df$Value) > 20))
      
      combined_df <- do.call(rbind, lapply(split(df, group), function(group_data) {
        max_density_row <- group_data[which.max(group_data$Density), ]
        max_density_row$Frequency <- sum(group_data$Frequency)
        max_density_row$Density <- sum(group_data$Density)
        return(max_density_row)
      }))
      call_df <- combined_df[which((combined_df$Density>=min)&(combined_df$Frequency>2)),]
      
      # Display the peaks
      call_df<-(call_df[order(call_df$Density, decreasing = T),])
      if(length(call_df$Frequency)!=0){
        add_df<- data.frame(matrix(ncol=6, nrow=(length(call_df$Value))))
        if(PAS_table[gene,"strand"]=="+"){
          add_df<-data.frame(
            chrom = gene_all$chrom[1], # Numeric column
            start = as.numeric(call_df$Value)-1, # Numeric column
            end = as.numeric(call_df$Value),
            name = gene_info[gene,"gene_name"],
            density= as.numeric(call_df$Density),
            strand = gene_info[gene,"strand"])
        }
        if(PAS_table[gene,"strand"]=="-"){
          add_df<-data.frame(
            chrom = gene_all$chrom[1], # Numeric column
            start = as.numeric(call_df$Value), # Numeric column
            end = as.numeric(call_df$Value)+1,
            name = gene_info[gene,"gene_name"],
            density= as.numeric(call_df$Density),
            strand = gene_info[gene,"strand"])
        }
        call_df$Density <-round(as.numeric(call_df$Density), 4)
        PAS_table[gene,c("number","PASs","PASs reads","PASs fraction")] <- c(length(call_df$Value),paste(call_df$Value, collapse = ","),
                                                                             paste(call_df$Frequency, collapse = ","),paste(call_df$Density, collapse = ","))
        return(list(PAS_table[gene,],add_df))
      }
    }
  }
  
  output <-pbmclapply(genes, PAS_fun, reads=reads, PAS_table, min=min, counts=counts,mc.cores = cores)
  bed_df <- do.call(rbind, lapply(output, `[[`, 2))
  PAS_list <- do.call(rbind, lapply(output, `[[`, 1))
  return(list(PAS_list,bed_df))
}