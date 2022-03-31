# from Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095
# modified by Todd Osmundson and Thomas Roehl 2021

## Configuration file for rnaseq_pipeline.sh
##
## Place this script in a working directory and edit it accordingly.
##
## The default configuration assumes that the user unpacked the 
## chrX_data.tar.gz file in the current directory, so all the input
## files can be found in a ./chrX_data sub-directory


### add this file as an argument when running rnaseq_pipeline.sh

# data directory -- location of the .fastq.gz files. All data files should be in a single directory. Please begin with the .fastq.gz file format or edit rnaseq_pipeline.sh to match the one you have.
DATADIR=

# sequencing batches -- if data was sequenced multiple times (such that each sample has multiple files for forward and reverse reads) they will need to be concatenated. Add the filename pattern that distinguishes each batch here.
SEQBATCHES=("_001" "_002")
# use this if there is only one batch of data: SEQBATCHES=("none")

# forward/reverse labels -- the filename patterns that distinguish forward and reverse samples, list the forward pattern first
FWDREV=("_R1" "_R2")

# working directory -- where to put files that are being worked on currently
WRKDIR=

# final directory -- where the files should go for storage
DESTDIR=

#how many CPUs to use on the current machine?
NUMCPUS=6


#### Program paths ####
#if these programs are not in any PATH directories, please edit accordingly:

FASTQC=/Applications/FastQC.app/Contents/MacOS/fastqc
FQTRIM=/Applications/fqtrim-0.9.7/fqtrim
TRIMMOMATIC=/Applications/Trimmomatic-0.38/trimmomatic-0.38.jar
INTERLEAVE=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_concat/fqtrim_output/interleave_pairs.py
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
BASEDIR="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq"
TEMPLOC="$BASEDIR/tmpdirec" #this will be relative to the output directory
################################### temploc should be removed, but is referenced in rnaseq_pipeline.sh

TRIMMOMATICADAPTERS="/Applications/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa" # pick the file that has your adapters
FASTQLOC="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_qc_done"
GENOMEIDX="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/genome_assemblies_all_files/ncbi-genomes-2021-11-01/GCA_011800155.1_ASM1180015v1/fvIndex"
PHENODATA="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/fv_sampling_data.csv"
BLASTDIR="/Applications/ncbi-blast-2.10.0+/bin"
UNIPROTDIR="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/uniprot"
RRNADIR="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/rRNA"
PANTHERLIBDIR="/Volumes/RAID_5_data_arra/Todd/Thomas_Roehl_RNASeq/panther/target4/famlib/rel/PANTHER16.0_altVersion/ascii/PANTHER16.0" # download and extract any one of the .tgz files (all contain same data) from http://data.pantherdb.org/ftp/panther_library/current_release/


## list of samples 
## (only paired reads, must follow _1.*/_2.* file naming convention)
############################## needs update -- these will be generated dynamically
reads1=(${FASTQLOC}/*_fwd.fastq)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_fwd/_rev}")
unpaired1=("${reads1[@]/_pairs_fwd/_unpaired}")
unpaired2=(${FASTQLOC}/*_1U.fqtrimmed.fq)
unpaired2=("${unpaired2[@]##*/}")
unpaired3=("${unpaired2[@]/_1U/_2U}")


## sample subset declarations
  # list samples you want to remove from analysis below
  # you can run the same analyses on different sets by creating a new removelist with and specifying all the samples to remove
  # then, add the new removelist in quotes to the REMOVELISTS variable
  # to identify and separate each set of files, add a unique name for each set in the SETNAMES variable
SETNAMES=("all" "300" "5k" "culnor" "priyou")
removelist1="-51 -48 -39"
removelist2="-51 -50 -48 -45 -39 -21"
removelist3="-12 -14 -17 -21 -22 -29 -30 -31 -33 -34 -38 -39 -40 -42 -45 -48 -50 -51 -52 -53 -56 -57"
removelist4="-10 -11 -13 -18-B -21 -25 -27 -28 -30 -31 -32 -34 -38 -39 -43 -44 -45 -48 -51 -52 -54-B -55 -56 -57"
removelist5="-51 -48 -39 -12 -14 -15 -16 -17 -19 -20 -22 -23 -24 -26 -29 -33 -35 -36 -37 -40 -41 -42 -46 -47 -49 -50 -53"
REMOVELISTS=("${removelist1}" "${removelist2}" "${removelist3}" "${removelist4}" "${removelist5}")


## list of data for BLAST
uniprot_file="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/uniprot/uniprot_agaricales_taxonomy_5338.fasta"
rrna_file="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/rRNA/fungi_18S_28S_ITS.fna"
