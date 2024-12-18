#' Visualize the 3'ends of all the reads in each sample for an individual gene
#' @param gene_reference information extracted from gtf file
#' @param reads information from samples
#' @param gene_list name or id of a list of genes of interest
#' @return plots showing the 3'end distribution of a given gene
#' @export

gene_explore <- function(gene_reference, reads, gene_list){
  gene_info <-gene_reference[[1]][,1:3]
  gene_info$short_gene_id <- gsub("\\.\\d+$", "", gene_info$gene_id)
  for (gene in gene_list){
    i <- na.omit(gene_info[apply(gene_info, 1, function(row) any(row == gene)),"gene_id"])
    gene_all <- subset(reads,gene_id==i)
    if (subset(gene_info, gene_id==i)$strand == "+") {
      gene_all <- subset(gene_all, strand == "+")
      Undiff <- subset(gene_all, treatment == "Undiff")$chromEnd
      Diff <- subset(gene_all, treatment == "Diff")$chromEnd
    }
    
    if (subset(gene_info, gene_id==i)$strand == "-") {
      gene_all <- subset(gene_all, strand == "-")
      Undiff <- subset(gene_all, treatment == "Undiff")$chromStart
      Diff <- subset(gene_all, treatment == "Diff")$chromStart
    }
    par(mfrow = c(1, 2))
    plot(ecdf(Undiff),xlim=c(min(c(Undiff,Diff)),max(c(Undiff,Diff))),main="Cumulative Curves",xlab=gene, ylab="Fraction")
    legend(x = "bottomright", box.col = "white", 
           bg ="white", box.lwd = 2 , title=" ", legend=c("sample1", "sample2"),  fill = c("black","pink"))
    plot(ecdf(Diff), add=T, col="pink")
    plot(density(Undiff,bw=10),xlim=c(min(c(Undiff,Diff)),max(c(Undiff,Diff))),ylim=c(0, max(c(density(Undiff,bw=10)$y,density(Diff,bw=10)$y))),
         main="PolyA sites",xlab=gene, ylab="Density")
    legend(x = "topright", box.col = "white", 
           bg ="white", box.lwd = 2 , title=" ", legend=c("sample1", "sample2"),  fill = c("black","pink"))
    lines(density(Diff,bw=10),col="pink")
    par(mfrow = c(1, 1))
  }
}