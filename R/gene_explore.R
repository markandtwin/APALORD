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
        variance <- max(APA_change_table$PAS_position)-min(APA_change_table$PAS_position)
        min_x <-min(APA_change_table$PAS_position)-variance*0.2
        max_x <- max(APA_change_table$PAS_position)+variance*0.2
        min_3end <-min(control_3end,experimental_3end)-variance*0.2
        max_3end <- max(control_3end,experimental_3end)+variance*0.2
        plot(ecdf(control_3end),xlim=c(min_x, max_x),col=rgb(0, 0, 0),
             main="Cumulative Curves",xlab=gene_annotation, ylab="Fraction")
        legend(x = "bottomright", lty=1,lwd=2,box.col = NA,
               bg ="white",  title=" ", legend=c(sample_names[1], sample_names[2]),  col = c(rgb(0, 0, 0),rgb(1, 0, 0.8, 0.5)))
        plot(ecdf(experimental_3end), add=T, col=rgb(1, 0, 0.8, 0.5))
        
        h_c <- hist(control_3end, breaks = max(control_3end)-min(control_3end), plot = FALSE)
        h_e <- hist(experimental_3end, breaks = max(experimental_3end)-min(experimental_3end), plot = FALSE)
        star_x <- APA_change_table$PAS_position  
        star_y <- (-0.02)
        
        plot(h_c$mids, h_c$density, type = "l", col =rgb(0, 0, 0), lwd=2,
             main="PolyA sites",xlab=gene_annotation, ylab="Density", 
             xlim = c(min_x,max_x), ylim = c(-0.05, max(h_c$density,h_e$density) * 1.2))
        legend(x = "topright",  lty=1,lwd=2,box.col = NA,
               bg =NA,  title=" ", legend=c(sample_names[1], sample_names[2]),  col = c(rgb(0, 0, 0),rgb(1, 0, 0.8, 0.5)))
        lines(h_e$mids,h_e$density,col=rgb(1, 0, 0.8, 0.5),lwd=2)
        text(star_x, star_y, labels = "PAS", cex = 0.7, col = "black")  
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