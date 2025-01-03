#' Visualize the 3'ends of all the reads in each sample for an individual gene and the PAUs for each polyA site from each group/sample
#' @param gene_reference information extracted from gtf file
#' @param reads information from samples
#' @param gene_list name or id of a list of genes of interest
#' @param APA_table output from APA_profile analysis
#' @param control which group in the data is used as the control group
#' @param experimental which group in the data is used as the experimental group
#' @param direct_RNA whether the data is direct RNAseq, TRUE or FALSE 
#' @return plots showing the 3'end distribution of a given gene
#' @export

gene_explore <- function(gene_reference, reads, gene_list, APA_table, control,experimental, direct_RNA=FALSE){
  gene_info <-gene_reference[,1:3]
  gene_info$short_gene_id <- gsub("\\.\\d+$", "", gene_info$gene_id)
  for (gene in gene_list){
    i <- na.omit(gene_info[apply(gene_info, 1, function(row) any(row == gene)),"gene_id"])
    gene_all <- subset(reads,gene_id==i)
    if(direct_RNA=="direct RNA"){
      gene_all <- subset(gene_all, strand==gene_info[gene,"strand"])
    }
    if (subset(gene_info, gene_id==i)$strand == "+") {
      control_3end <- subset(gene_all, treatment == control)$chromEnd
      experimental_3end <- subset(gene_all, treatment == experimental)$chromEnd
    }
    if (subset(gene_info, gene_id==i)$strand == "-") {
      control_3end <- subset(gene_all, treatment == control)$chromStart
      experimental_3end <- subset(gene_all, treatment == experimental)$chromStart
    }
    APA_change_table<- data.frame(PAS_position=as.numeric(unlist(strsplit(APA_table[i,]$PASs, split = ","))),
                                  PAU_changes=as.numeric(unlist(strsplit(APA_table[i,]$PAU_changes, split = ","))),
                                  control_PAUs=as.numeric(unlist(strsplit(APA_table[i,]$PAUs_control, split = ","))),
                                  experiment_PAUs=as.numeric(unlist(strsplit(APA_table[i,]$PAUs_experimental, split = ","))))
    colors <- ifelse(APA_change_table$PAU_changes == max(APA_change_table$PAU_changes), "red", 
                     ifelse(APA_change_table$PAU_changes == min(APA_change_table$PAU_changes), "blue", "gray"))
    gene_annotation <-paste(APA_table[i,"gene_name"],"on",APA_table[i,"chrom"],"(",APA_table[i,"strand"], ")",sep=" ")
    par(mar = c(6, 5, 4, 2) + 0.1,mfrow = c(2, 2))
    layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
    plot(ecdf(control_3end),xlim=c(min(c(control_3end,experimental_3end)),max(c(control_3end,experimental_3end))),
         main="Cumulative Curves",xlab=gene_annotation, ylab="Fraction")
    legend(x = "bottomright", box.col = "white", 
           bg ="white", box.lwd = 2 , title=" ", legend=c(control, experimental),  fill = c("black","pink"))
    plot(ecdf(experimental_3end), add=T, col="pink")
    plot(density(control_3end,bw=10),xlim=c(min(c(control_3end,experimental_3end)),max(c(control_3end,experimental_3end))),
         ylim=c(0, max(c(density(control_3end,bw=10)$y,density(experimental_3end,bw=10)$y))),
         main="PolyA sites",xlab=gene_annotation, ylab="Density")
    legend(x = "topright", box.col = "white", 
           bg ="white", box.lwd = 2 , title=" ", legend=c(control, experimental),  fill = c("black","pink"))
    lines(density(experimental_3end,bw=10),col="pink")
    barplot(APA_change_table$PAU_changes, 
            names.arg = APA_change_table$PAS_position,
            main = paste("PAU difference at each PolyA site (", experimental, "-", control,")", sep=" "),
            xlab = "PAS coordinates",
            ylab = "Delta_PAU (%)",
            col = colors,
            ylim = c(2*min(APA_change_table$PAU_changes),2*max(APA_change_table$PAU_changes)),
            sub = gene_annotation)
    par(mfrow = c(1, 1))
  }
}