#' Visualize the 3'ends of all the reads in each sample for an individual gene and the PAUs for each polyA site from each group/sample
#' @import data.table
#' @param gene_reference information extracted from gtf file
#' @param reads information from samples
#' @param gene_list name or id of a list of genes of interest
#' @param APA_table output from APA_profile analysis
#' @param direct_RNA whether the data is direct RNAseq, TRUE or FALSE 
#' @return plots showing the 3'end distribution of a given gene
#' @export

gene_explore <- function(gene_reference, reads, gene_list, APA_table, direct_RNA=FALSE){
  names <-grep("^PAUs_", colnames(APA_table), value = TRUE)
  sample_names <- sub("^PAUs_", "", names)
  gene_info <-gene_reference
  gene_info$short_gene_id <- gsub("\\.\\d+$", "", gene_info$gene_id)
  for (gene in gene_list){
    i <- gene_info[apply(gene_info, 1, function(x) any(grepl(gene, x)))]$gene_id
    for (id in i){
      gene_all <- reads[as.character(gene_id)==id]
      if(direct_RNA){
        gene_all <- gene_all[strand == gene_info[id, strand]]
      }
      if (gene_info[id, strand] == "+") {
        control_3end <- subset(gene_all, treatment == sample_names[1])$chromEnd
        experimental_3end <- subset(gene_all, treatment == sample_names[2])$chromEnd
      }
      if (gene_info[id, strand] == "-") {
        control_3end <- subset(gene_all, treatment == sample_names[1])$chromStart
        experimental_3end <- subset(gene_all, treatment == sample_names[2])$chromStart
      }
      if(nrow(APA_table[gene_id==id])>0){
        PAU_table <-APA_table[gene_id==id,..names]
        APA_change_table<- data.frame(PAS_position=as.numeric(unlist(strsplit(APA_table[gene_id==id]$PAS_coordinates, split = ","))),
                                      PAU_changes=as.numeric(unlist(strsplit(APA_table[gene_id==id]$PAU_changes, split = ","))),
                                      control_PAUs=as.numeric(unlist(strsplit(as.character(PAU_table[,1]), split = ","))),
                                      experiment_PAUs=as.numeric(unlist(strsplit(as.character(PAU_table[,2]), split = ","))))
        APA_change_table[,"Colors"] <- ifelse(APA_change_table$PAU_changes > 0 & APA_change_table$PAU_changes == max(APA_change_table$PAU_changes), "red",
                                              ifelse(APA_change_table$PAU_changes < 0 & APA_change_table$PAU_changes == min(APA_change_table$PAU_changes), "blue", "gray"))
        gene_annotation <-paste(APA_table[gene_id==id]$gene_name,"on",APA_table[gene_id==id]$chrom,"(",APA_table[gene_id==id]$strand, ")",sep=" ")
        APA_change_table <- APA_change_table[order(APA_change_table$PAS_position, decreasing = F), ]
        par(mar = c(6, 5, 4, 2) + 0.1,mfrow = c(3, 1))
        layout(matrix(c(1, 2, 3),3, 1, byrow = TRUE))
        plot(ecdf(control_3end),xlim=c(min(APA_change_table$PAS_position)-500,max(APA_change_table$PAS_position)+500),
             main="Cumulative Curves",xlab=gene_annotation, ylab="Fraction")
        legend(x = "bottomright", box.col = "white", 
               bg ="white", box.lwd = 2 , title=" ", legend=c(sample_names[1], sample_names[2]),  fill = c("black","pink"))
        plot(ecdf(experimental_3end), add=T, col="pink")
        plot(density(control_3end,bw=1),xlim=c(min(APA_change_table$PAS_position)-500,max(APA_change_table$PAS_position)+500),
             ylim=c(0, max(c(density(control_3end,bw=1)$y,density(experimental_3end,bw=1)$y))),
             main="PolyA sites",xlab=gene_annotation, ylab="Density")
        legend(x = "topright", box.col = "white", 
               bg ="white", box.lwd = 2 , title=" ", legend=c(sample_names[1], sample_names[2]),  fill = c("black","pink"))
        lines(density(experimental_3end,bw=1),col="pink")
        barplot(APA_change_table$PAU_changes, 
                names.arg = APA_change_table$PAS_position,
                main = paste("PAU difference at each PolyA site (", sample_names[2], "-", sample_names[1],")", sep=" "),
                xlab = "PAS coordinates",
                ylab = "Delta_PAU (%)",
                col = APA_change_table$Colors,
                ylim = c(2*min(APA_change_table$PAU_changes,0),2*max(APA_change_table$PAU_changes,0)),
                sub = gene_annotation)
        par(mfrow = c(1, 1))
      }
    }
  }
}