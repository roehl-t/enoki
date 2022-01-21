# enoki
RNA-seq analysis pipeline (HISAT2, StringTie, Ballgown, heatmaps, bar plots, Gene Ontology) used during enoki transcriptome research

External Requirements:
FastQC
fqtrim
trimmomatic
HISAT2
StringTie
samtools
gffread
Ballgown
BLAST+ suite
Docker
BUSCO
PantherScore

Inputs:
demultiplexed Illumina sequences
reference genome
csv file with phenotypic data for each sample

Outputs:
PCA plots
named lists of differentially expressed genes
heatmaps
PANTHER mappings for GO analysis


Overview:
Before running the main pipeline, run the BASH script fastqc_batch.sh to perform sequence trimming to remove primers and low-quality sequences. Requires FastQC, fqtrim, and trimmomatic. Calls interleave_pairs.py (python3).

Prior to running the main script, a relevant rRNA database and a relevant protein database must be saved to the local environment. Docker must also be set up with an instance of BUSCO.

After sequence trimming, all external programs and files in the enoki repository must be listed in the rnaseq_pipeline.config.sh file except: FastQC, fqtrim, trimmomatic, interleave_pairs.py, fastqc_batch.sh, and expression_barplot.R. The locations of required files and base filenames of included samples must also be updated.

The main script for this pipeline is rnaseq_pipeline_merge.sh, a BASH shell script, and follows this order of events:
1) First, it checks for the presence of all required files. It then runs HISAT2 and StringTie and performs appropriate file conversions.
2) Prior to abundance estimation in StringTie, the script removes files from the analysis. Update these lists with any files you wish to remove (matches substrings in filenames). The remainder of the pipeline is then run separately on each specified subset of files.
3) StringTie abundance estimation is carried out and the transcriptome is created.
4) BLASTN searches the transcriptome against a ribosomal database. Any transcripts that match ribosomal genes are removed from analysis.
5) BUSCO analysis assesses transcriptome quality.
6) Ballgown analysis (run in the R environment) identifies differentially expressed genes. During this step, the rnaseq_ballgown.R script excludes samples with low numbers of detected transcripts. Adjust this value as needed. PCA plots and lists of differentially expressed genes are created during this step. The desired PCA comparisons must be specified manually in rnaseq_ballgown.R using names found in the pneotype data csv.
7) Next, BLASTX searhces a protein database to find names for the genes in the lists of differentially expressed genes.
8) Heatmaps are created to display differentially expressed genes from each sample. Heatmaps for all genes as well as user-specified subsets (based on gene, sample, or phenotypic data) can be created.
9) GO enrichment is performed. First, the transcriptome is aligned to the PANTHER hmms by PantherScore. Next, the sequences for proteins matching differentially expressed genes are retreived using eFetch from NCBI. The retreived sequences are mapped to PANTHER hmms by PantherScore. The resulting mappings are saved to a file which can be uploaded for GO analysis.

expression_barplot.R can be used to manually create bar graphs of specific genes. Subsets can be created based on gene or phenotypic data. If you need to display expression patterns of more than 6 genes, consider using a heatmap instead. heatmap.R can be configured manually as well.

PCA plots can be manually created after the pipeline is run by using manual controls in the rnaseq_ballgown.R file.
