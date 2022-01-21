# from Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095
# modified by Todd Osmundson and Thomas Roehl 2021

## Configuration file for rnaseq_pipeline.sh
##
## Place this script in a working directory and edit it accordingly.
##
## The default configuration assumes that the user unpacked the 
## chrX_data.tar.gz file in the current directory, so all the input
## files can be found in a ./chrX_data sub-directory

#how many CPUs to use on the current machine?
NUMCPUS=6

#### Program paths ####

## optional BINDIR, using it here because these programs are installed in a common directory
#BINDIR=/usr/local/bin
#HISAT2=$BINDIR/hisat2
#STRINGTIE=$BINDIR/stringtie
#SAMTOOLS=$BINDIR/samtools

#if these programs are not in any PATH directories, please edit accordingly:
#HISAT2=$(which hisat2)
#STRINGTIE=$(which stringtie)
#SAMTOOLS=$(which samtools)
HISAT2=/Applications/hisat2-2.2.1/hisat2
STRINGTIE=/Applications/stringtie-1.3.6.OSX_x86_64/stringtie
SAMTOOLS=/Applications/samtools-1.8/samtools
BALLGOWN=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/rnaseq_ballgown.R
GFFREAD=/Applications/gffread-0.12.7.OSX_x86_64/gffread
BLASTXAPP=/Applications/ncbi-blast-2.10.0+/bin/blastx
BLASTNAPP=/Applications/ncbi-blast-2.10.0+/bin/blastn
MATCHGENES=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/match_mstrng_to_gene.py
NAMEGENES=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/name_genes2.R
REMOVERRNA=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/remove_rRNA.R
FIXUNIPROT=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/extract_uniprot_title_info.R
HEATMAP=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/heatmap.R
UNIPROTFASTA=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/map_UniProtKB_FASTA.py
PANTHERSCORE=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/panther/pantherScore2.2/pantherScore2.2.pl
PANTHERFPKM=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/map_panther_fpkm.R
#### File paths for input data
### Full absolute paths are strongly recommended here.
## Warning: if using relatives paths here, these will be interpreted 
## relative to the  chosen output directory (which is generally the 
## working directory where this script is, unless the optional <output_dir>
## parameter is provided to the main pipeline script)

## Optional base directory, if most of the input files have a common path
# BASEDIR="/home/johnq/RNAseq_protocol/chrX_data"
#BASEDIR=$(pwd -P)/HiSat2_analysis_Dpoly
BASEDIR="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq"
FASTQLOC="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_qc_done"
#FASTQLOC="$BASEDIR/samples"
GENOMEIDX="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/genome_assemblies_all_files/ncbi-genomes-2021-11-01/GCA_011800155.1_ASM1180015v1/fvIndex"
#GTFFILE="$BASEDIR/genes/Dpolymorpha_annotations.gtf"
PHENODATA="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/fv_sampling_data.csv"
BLASTDIR="/Applications/ncbi-blast-2.10.0+/bin"
UNIPROTDIR="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/uniprot"
RRNADIR="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/rRNA"

TEMPLOC="$BASEDIR/tmpdirec" #this will be relative to the output directory

## list of samples 
## (only paired reads, must follow _1.*/_2.* file naming convention)
reads1=(${FASTQLOC}/*_fwd.fastq)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_fwd/_rev}")
#reads1=(${FASTQLOC}/*_R1_*)
#reads1=("${reads1[@]##*/}")
#reads2=("${reads1[@]/_R1_/_R2_}")
unpaired1=("${reads1[@]/_pairs_fwd/_unpaired}")
unpaired2=(${FASTQLOC}/*_1U.fqtrimmed.fq)
unpaired2=("${unpaired2[@]##*/}")
unpaired3=("${unpaired2[@]/_1U/_2U}")

## list of data for BLAST
uniprot_file="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/uniprot/uniprot_agaricales_taxonomy_5338.fasta"
rrna_file="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/rRNA/fungi_18S_28S_ITS.fna"
