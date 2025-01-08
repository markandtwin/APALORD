library(dplyr)
library(stringr)
library(tidyr)
library(pbmcapply)
library(readr)
library(data.table)
library(devtools)
library(roxygen2)
devtools::install("../KSAPA")
library(KSAPA)
#detach("package:KSAPA", unload = TRUE)

extdata_path <- system.file("extdata",package = "KSAPA")
gtf.file <- paste0(extdata_path,"/hg38_chr20.gtf")
gtf.file <- "~/Desktop/Human_annotation/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"
gene_reference <- load_gtf(gtf.file,cores = 7)


sample1 <- paste0(extdata_path,"/D0_isoquant")
sample2 <- paste0(extdata_path,"/D7_isoquant")
sample1 <- "~/Desktop/KSAPA/WT_D0_dRNA.no_sec.isoquant/OUT/"
sample2 <- "~/Desktop/KSAPA/WT_D7_dRNA_4.isoquant/"
reads <- load_samples(sample1,sample2)

PAS_data <- PAS_calling(gene_reference,reads,min=0.01,cores=7)
PAS_df<- PAS_data[[1]]
write.table(PAS_data[[2]],file="./PASs_hES_0.01_D7_D0.bed",quote = F,col.names = F, row.names = F,sep = "\t")


APA_data <- APA_profile(gene_reference,reads,cores = 7)
APA_table <- left_join(APA_data,PAS_df[,c("gene_id","number","Tandem_APA")],by="gene_id")
APA_gene_table <- APA_plot(APA_data)
TAPA_gene_table <- APA_plot(subset(APA_table,number>1))
#HAPA_table <- (subset(APA_table,number==1))

gene_explore(gene_reference, reads,c("EIF4B","ENSG00000101146"))
