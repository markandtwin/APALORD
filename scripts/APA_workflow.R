
devtools::install("../KSAPA")
library(KSAPA)
#detach("package:KSAPA", unload = TRUE)

extdata_path <- system.file("extdata",package = "KSAPA")
gtf_file <- paste0(extdata_path,"/hg38_chr21.gtf.gz")
gene_reference <- load_gtf(gtf_file,cores = 5)



sample1 <- c(paste0(extdata_path,"/D0/rep1"),paste0(extdata_path,"/D0/rep2"),paste0(extdata_path,"/D0/rep3"))
sample2 <- c(paste0(extdata_path,"/D7/rep1"),paste0(extdata_path,"/D7/rep2"),paste0(extdata_path,"/D7/rep3"))
reads <- load_samples(sample1,sample2, group1="D0",group2="D7")

PAS_data <- PAS_calling(gene_reference,reads,cores=5,direct_RNA = T)
write.table(PAS_data,file="../test_data/PAS_bed_hES_0.01_D7_D0.bed",quote = F,col.names = F, row.names = F,sep = "\t")


PAU_data <- PAU_by_sample(gene_reference,reads,cores=7,direct_RNA = T)
write.table(PAU_data,file="../test_data/PAU_by_sample_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")

#####################################APA analysis####################################################
PAU_test_data <- PAU_test(PAU_data, reads, P_cutoff = 0.05)
write.table(PAU_test_data,file="../test_data/hESC_PAU_test_data_all.tsv", quote = F,col.names = T, row.names = F,sep = "\t")

dPAU_test_data <- KSAPA::end_PAS_examine(PAU_data, reads, P_cutoff = 0.2,
                                         control = "D0", experimental = "D7", position = "distal", type = "FC")
pPAU_test_data <- KSAPA::end_PAS_examine(PAU_data, reads, P_cutoff = 0.2,
                                         control = "D0", experimental = "D7", position = "proximal", type = "delta")


APA_data <- APA_profile(gene_reference,reads,control="D0",experimental="D7", cores = 5, direct_RNA = T)


#########now the APA_plot function will save the plot to the directory insteading of showing in the plot panel#########
APA_gene_table <- APA_plot(APA_data)
write.table(APA_data,file="../test_data/APA_data_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")
write.table(APA_gene_table,file="../test_data/APA_gene_table_hES_0.01_D7_D0.tsv",quote = F,col.names = T, row.names = F,sep = "\t")


#################explore single gene APA change##########################
gene_explore(gene_reference, reads,c("GART","ENSG00000185658","ZBTB21"), APA_table = APA_data, direct_RNA = T)

####################Within PAS CSH shift######################
CSH_data <- CSH_profile(reads, PAU_data,control = "D0",experimental="D7", cores = 5, direct_RNA = T)
CSH_test <- CSH_plot(CSH_data)
write.table(CSH_test,file="../test_data/CSH_gene_table_hES_D0_D7.tsv",quote = F,col.names = T, row.names = F,sep = "\t")



