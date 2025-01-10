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
gtf.file <- paste0(extdata_path,"/hg38_chr20.gtf.gz")
#gtf.file <- "~/Desktop/Human_annotation/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"
gene_reference <- load_gtf(gtf.file,cores = 7)


sample1 <- paste0(extdata_path,"/D0_isoquant")
sample2 <- paste0(extdata_path,"/D7_isoquant")
#sample1 <- "~/Desktop/KSAPA/WT_D0_dRNA.no_sec.isoquant/"
#sample2 <- "~/Desktop/KSAPA/WT_D7_dRNA_4.isoquant/"
reads <- load_samples(sample1,sample2, group1="D0",group2="D7")

#PAS_data <- PAS_calling(gene_reference,reads,min=0.01,cores=7)
PAU_data <- PAU_by_sample(gene_reference,reads,cores=7,direct_RNA = T)

write.table(PAU_data,file="../../PAU_by_sample_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")


APA_data <- APA_profile(gene_reference,reads,control="D0",experimental="D7", cores = 7, direct_RNA = T)

APA_gene_table <- APA_plot(APA_data)

gene_explore(gene_reference, reads,c("JAG1","ENSG00000088766","SLC2A10"),
             control="D0",experimental="D7", APA_table = APA_data, direct_RNA = T)
write.table(APA_data,file="../../APA_data_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")
write.table(APA_gene_table,file="../../APA_gene_table_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")

