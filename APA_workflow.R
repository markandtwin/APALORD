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
gtf_file <- paste0(extdata_path,"/hg38_chr21.gtf.gz")
gene_reference <- load_gtf(gtf_file,cores = 5)


sample1 <- paste0(extdata_path,"/D0_isoquant")
sample2 <- paste0(extdata_path,"/D7_isoquant")
reads <- load_samples(sample1,sample2, group1="D0",group2="D7")

PAS_data <- PAS_calling(gene_reference,reads,cores=5,direct_RNA = T)
write.table(PAS_data,file="../../PAS_bed_hES_0.01_D7_D0.bed",quote = F,col.names = F, row.names = F,sep = "\t")


PAU_data <- PAU_by_sample(gene_reference,reads,cores=7,direct_RNA = T)
write.table(PAU_data,file="../../PAU_by_sample_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")


APA_data <- APA_profile(gene_reference,reads,control="D0",experimental="D7", cores = 5, direct_RNA = T)


#########now the APA_plot function will save the plot to the directory insteading of showing in the plot panel#########
APA_gene_table <- APA_plot(APA_data)

gene_explore(gene_reference, reads,c("GART","ENSG00000185658","ZBTB21"),
             control="D0",experimental="D7", APA_table = APA_data, direct_RNA = T)
write.table(APA_data,file="../../APA_data_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")
write.table(APA_gene_table,file="../../APA_gene_table_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")

