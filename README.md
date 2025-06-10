# APALORD: Alternative Polyadenylation Analysis of LOng-ReaDs (APALORD) RNA-seq

`APALORD` is an R package for analyzing Alternative Polyadenylation (APA) in long-read RNA-seq data. It provides tools to identify polyadenylation sites (PAS), 
quantify polyadenylation site usage 
(PAU), perform differential APA analysis, profile cleavage site heterogeneity (CSH), and visualize APA changes across conditions. Optimized for multi-core 
processing, `APALORD` is ideal for APA 
analysis utilizing long read RNA-seq across conditions.

## Features

- **PAS Calling**: Detect polyadenylation sites from RNA-seq reads.
- **PAU Quantification**: Measure polyadenylation usage per sample.
- **Differential APA Analysis**: Identify significant PAU changes and analyze distal/proximal PAS shifts.
- **APA Profiling**: Generate and visualize APA profiles between conditions.
- **Single Gene Exploration**: Examine APA changes for specific genes.


## Installation

### Prerequisites

- R (>= 4.0)
- R package: `devtools` (for installation)
- Input files:
  - GTF annotation file (e.g., `hg38_chr21.gtf.gz`)
  - Preprocessed long read RNA-seq data (e.g., output of IsoQuant) for two groups


### Install APALORD

1. **Install from GitHub**:

   ```R
   install.packages("devtools")
   devtools::install_github("markandtwin/APALORD",build_vignettes = TRUE)
   ```

2. **Verify Installation and Browse the Vignettes**:

   ```R
   library(APALORD)
   browseVignettes(package = "APALORD")
   ```

Alternatively, clone the repository and install locally:

```bash
git clone https://github.com/markandtwin/APALORD.git
cd APALORD
```
And then install it in R and browse the Vignettes:

   ```R
   devtools::install("PATH to /APALORD", build_vignettes = TRUE)
   library(APALORD)
   browseVignettes(package = "APALORD")
   ```
## Usage

The following example demonstrates the `APALORD` workflow using the package’s example data for hESC (D0) and hESC-derived neuron (D7) samples. Each step is 
presented as a separate R code chunk, which 
can be combined into a script (`scripts/APA_analysis.R`) for a pipeline execution. 

### Example Workflow

To run the full workflow, create a script (`scripts/APA_analysis.R`) by combining the chunks below, or execute them individually in an R session.

#### Step 0: Initialize Environment

```R
# Initialize APALORD workflow
# - Load package and create output directory
library(APALORD)

# Set output directory
out_dir <- "./results"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
```

#### Step 1: Load GTF Annotation

```R
# Load GTF annotation
# - Use hg38_chr21.gtf.gz for gene annotations
# - Validate file existence and use 5 cores for parallel processing
extdata_path <- system.file("extdata", package = "APALORD")
gtf_file <- paste0(extdata_path, "/hg38_chr21.gtf.gz")
gene_reference <- load_gtf(gtf_file, cores = 5)
```

#### Step 2: Load Samples

```R
# Load samples
# - Define BAM file paths for D0 (control) and D7 (experimental) groups
# - Validate file existence for all replicates
sample1 <- c(paste0(extdata_path, "/D0/rep1"), paste0(extdata_path, "/D0/rep2"), paste0(extdata_path, "/D0/rep3"))
sample2 <- c(paste0(extdata_path, "/D7/rep1"), paste0(extdata_path, "/D7/rep2"), paste0(extdata_path, "/D7/rep3"))
reads <- load_samples(sample1, sample2, group1 = "D0", group2 = "D7")
```

#### Step 3: PAS Calling (optional)

```R
# PAS calling
# - Identify PASs from reads data
# - Use direct_RNA=TRUE for direct RNA-seq

PAS_data <- PAS_calling(gene_reference, reads, cores = 5, direct_RNA = TRUE)
pas_file <- file.path(out_dir, "PAS_D7_D0.bed")
write.table(PAS_data, file = pas_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
```
This PAS calling output file is in .bed format, so users can load it into IGV or UCSC genome browser to check it. It's not required by the downstream analysis, so 
this step is optional. 

#### Step 4: PAU Quantification
##### 4.1: PAU Quantification

```R
# PAU quantification
# - Quantify polyadenylation usage per sample
# - Use multiple cores for efficiency
message("Quantifying PAU by sample...")
PAU_data <- PAU_by_sample(gene_reference, reads, cores = 7, direct_RNA = TRUE)
pau_file <- file.path(out_dir, "PAU_by_sample_D7_D0.tsv")
write.table(PAU_data, file = pau_file, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
```

##### 4.2: Differential PAU Analysis

```R
# Differential PAU analysis
# - Test for significant PAU changes with P_cutoff=0.05

PAU_test_data <- PAU_test(PAU_data, reads, P_cutoff = 0.05)
pau_test_file <- file.path(out_dir, "hESC_PAU_test_data_all.tsv")
write.table(PAU_test_data, file = pau_test_file, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

```
This step is using DEXSeq to call differentially expressed PAS between the two groups

##### 4.3: Distal/Proximal PAS Shifts (optional)

```R
# Distal/proximal PAS shifts
# - Examine shifts in distal and proximal PAS usage

dPAU_test_data <- APALORD::end_PAS_examine(PAU_data, reads, P_cutoff = 0.2,
                                         control = "D0", experimental = "D7", position = "distal", type = "FC")
pPAU_test_data <- APALORD::end_PAS_examine(PAU_data, reads, P_cutoff = 0.2,
                                         control = "D0", experimental = "D7", position = "proximal", type = "delta")
```

#### Step 5: APA Profiling

```R
# APA profiling
# - Generate APA profiles and visualize changes
message("Profiling APA changes...")
APA_data <- APA_profile(gene_reference, reads, control = "D0", experimental = "D7", cores = 5, direct_RNA = TRUE)
APA_gene_table <- APA_plot(APA_data)
# The APA_plot function generates a Vocanol plot showing the transcriptome-wide APA changes.

apa_file <- file.path(out_dir, "APA_data_hES_D7_D0.tsv")
apa_table_file <- file.path(out_dir, "APA_gene_table_hES_D7_D0.tsv")
write.table(APA_data, file = apa_file, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(APA_gene_table, file = apa_table_file, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

```

#### Visualization: Single Gene Exploration

```R
# Single gene exploration
# - Explore APA changes for specific genes
message("Exploring APA changes for genes: GART, ENSG00000185658, ZBTB21")
gene_explore(gene_reference, reads, c("GART", "ENSG00000185658", "ZBTB21"), APA_table = APA_data, direct_RNA = TRUE)
```



### Running the Workflow

To execute the full workflow, combine the chunks into `scripts/APA_analysis.R` and run:

```bash
Rscript scripts/apa_analysis.R
```

Alternatively, run each chunk individually in an R session, ensuring dependencies (e.g., `gene_reference`, `reads`, `PAU_data`, `APA_data`) are available from 
previous steps.


## Output Files

The workflow generates files in the `results/` directory:

- `PAS_D7_D0.bed`: Polyadenylation sites in BED format.
- `PAU_by_sample_D7_D0.tsv`: PAU per sample.
- `hESC_PAU_test_data_all.tsv`: Differential PAU test results.
- `APA_data_hES_D7_D0.tsv`: APA profiles.
- `APA_gene_table_hES_0.01_D7_D0.tsv`: APA gene table from plots.
- Plots from `APA_plot` and `CSH_plot` in `Plots`.



## Notes

- **Custom Data**: Use Custom own IsoQuant output files by updating paths. Ensure that the `OUT.corrected_reads.bed` and `OUT.read_assignments.tsv` files could be 
found in the PATH provided for each 
sample:

  ```R
  sample1 <- c("/path/to/Condition A rep1 IsoQuant OUT", "/path/to/Condition A rep2 IsoQuant OUT", "/path/to/Condition A rep3 IsoQuant OUT/")
  sample2 <- c("/path/to/Condition B rep1 IsoQuant OUT", "/path/to/Condition B rep2 IsoQuant OUT", "/path/to/Condition B rep3 IsoQuant OUT")
  ```

- **IsoQuant Output File Preparation**: Ensure BAM files are sorted and indexed. For example: 

  ```bash
  samtools sort sample1.bam -o sample1.sorted.bam
  samtools index sample1.sorted.bam
  isoquant.py --out sample1.isoquant --genedb /path/to/.annotation.gtf --reference /path/to/genome.fasta --complete_genedb --data_type nanopore --bam 
sample1.sorted.bam --threads 7 
--gene_quantification unique_only --transcript_quantification unique_only --splice_correction_strategy default_ont --no_secondary
  ```


- **Gene IDs**: Verify gene IDs in `gene_explore` (e.g., `GART`, `ENSG00000185658`, `ZBTB21`) match the GTF file. 

- **Direct RNA Sequencing**: The `direct_RNA = TRUE` parameter is set for direct RNA-seq. Default set to `FALSE` for other long-read protocols.

## Troubleshooting

- **Installation Issues**: If `devtools::install_github()` fails, verify `DESCRIPTION` and `NAMESPACE`. Run `R CMD check .` for diagnostics.

- **File Paths**: Check GTF and read file paths:

  ```R
  file.exists(gtf_file)
  ```


- **Plot Output**: `APA_plot` or `CSH_plot` does not save plots. To save the plot:

  ```R
  pdf("APA_plot.pdf",7,7)
  APA_gene_table <- APA_plot(APA_data)
  dev.off()  
  ```


## Citation

If you use `APALORD` in your research, please cite [INSERT CITATION OR PLACEHOLDER].

