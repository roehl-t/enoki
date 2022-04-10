# from Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095
# modified by Todd Osmundson and Thomas Roehl 2021

## Configuration file for rnaseq_pipeline.sh
##
## Place this script in a working directory and edit it accordingly.
##
## To run the pipeline, run the rnaseq_pipeline.sh script with this file as the first argument
##    ex: username:~$ bash ./rnaseq_pipeline.sh ./rnaseq_pipeline.config.sh


#### User Options ####

# location for all log files
LOGLOC=

# data directory -- location of the .fastq.gz files. All data files should be in a single directory. Please begin with the .fastq.gz file format or edit rnaseq_pipeline.sh to match the one you have.
DATADIR=

# sequencing batches -- if data was sequenced multiple times (such that each sample has multiple files for forward and reverse reads) they will need to be concatenated. Add the filename pattern that distinguishes each batch here.
SEQBATCHES=("_001" "_002")
# use this if there is only one batch of data: SEQBATCHES=("none")

# forward/reverse labels -- the filename patterns that distinguish forward and reverse samples, list the forward pattern first
FWDREV=("_R1" "_R2")

# use unpaired reads? -- either Y to use the unpaired reads or N to ignore them
USEUNPAIRED=Y

# working directory -- where to put files that are being worked on currently
  # hint: use an SSD for speed improvement
WRKDIR=

# final directory -- where the files should go for storage
  # hint: use a large drive
  ## must have a different name than the working directory
DESTDIR=

# how many CPUs to use on the current machine?
NUMCPUS=6


#### Program paths ####
#if these programs are not in any PATH directories, please edit accordingly:

## paths for programs not in enoki
FASTQC=/Applications/FastQC.app/Contents/MacOS/fastqc
FQTRIM=/Applications/fqtrim-0.9.7/fqtrim
TRIMMOMATIC=/Applications/Trimmomatic-0.38/trimmomatic-0.38.jar
HISAT2=/Applications/hisat2-2.2.1/hisat2
STRINGTIE=/Applications/stringtie-1.3.6.OSX_x86_64/stringtie
SAMTOOLS=/Applications/samtools-1.8/samtools
BALLGOWN=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/rnaseq_ballgown.R
GFFREAD=/Applications/gffread-0.12.7.OSX_x86_64/gffread
BLASTXAPP=/Applications/ncbi-blast-2.10.0+/bin/blastx
BLASTNAPP=/Applications/ncbi-blast-2.10.0+/bin/blastn
PANTHERSCORE=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/panther/pantherScore2.2/pantherScore2.2.pl

## paths for programs included in enoki
scripts=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts
# all these should be in the same folder
  # if they are, simply update scripts= with the name of the folder
  # if scripts are in separate folders, update all the variables below
INTERLEAVE=${scripts}/interleave_pairs.py
REMOVERRNA=${scripts}/remove_rRNA.R
TRANSCRIPTCOUNTER=${scripts}/transcript_counter.R
MATCHGENES=${scripts}/match_mstrng_to_gene.py
NAMEGENES=${scripts}/name_genes2.R
FIXUNIPROT=${scripts}/extract_uniprot_title_info.R
HEATMAP=${scripts}/heatmap.R
UNIPROTFASTA=${scripts}/map_UniProtKB_FASTA.py
PANTHERFPKM=${scripts}/map_panther_fpkm.R


#### File paths for input data
### Full absolute paths are strongly recommended here.
## Warning: if using relatives paths here, these will be interpreted 
## relative to the  chosen output directory (which is generally the 
## working directory where this script is, unless the optional <output_dir>
## parameter is provided to the main pipeline script)

TRIMMOMATICADAPTERS="/Applications/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa" # pick the file that has your adapters
FASTQLOC="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_qc_done"
GENOME="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/genome_assemblies_all_files/ncbi-genomes-2021-11-01/GCA_011800155.1_ASM1180015v1/GCA_011800155.1_ASM1180015v1_genomic.fna"
GENOMEIDX="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/genome_assemblies_all_files/ncbi-genomes-2021-11-01/GCA_011800155.1_ASM1180015v1/fvIndex"
# notes for pheotype data file (see example file)
    # do not use spaces in your column names!
    # must include a column labeled "ids" that matches the sample IDs used in the filenames
PHENODATA="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/fv_sampling_data.csv"
BLASTDIR="/Applications/ncbi-blast-2.10.0+/bin"
PANTHERLIBDIR="/Volumes/RAID_5_data_arra/Todd/Thomas_Roehl_RNASeq/panther/target4/famlib/rel/PANTHER16.0_altVersion/ascii/PANTHER16.0" # download and extract any one of the .tgz files (all contain same data) from http://data.pantherdb.org/ftp/panther_library/current_release/


## list of data for BLAST
RRNAFILE="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/rRNA/fungi_18S_28S_ITS.fna"
UNIPROTFILE="/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/uniprot/uniprot_agaricales_taxonomy_5338.fasta"


## sample subset declarations
  # list samples you want to remove from analysis below, separated by a space
  # you can run the same analyses on different sets by creating a new removelist and specifying all the samples to remove
    # samples will be identified by partial matches to the filenames, so you must include at a minimum a short unique pattern to identify the sample
  # then, add the new removelist in quotes to the REMOVELISTS variable
  # to identify and separate each set of files, add a unique name for each set in the SETNAMES variable
  # to automatically remove samples containing below a certain threshold number of sequences, add "auto-####" to SETNAMES
  # to run the analysis on all samples, add "auto-all" to SETNAMES
  # do not add a removelist for "auto-###" or "auto-all" setnames -- we suggest putting all auto-##### setnames at the end of the SETNAMES list to be sure your manual removelists are in the correct order
SETNAMES=("auto-all" "300" "5k" "culnor" "priyou")
removelist2="-51 -50 -48 -45 -39 -21"
removelist3="-12 -14 -17 -21 -22 -29 -30 -31 -33 -34 -38 -39 -40 -42 -45 -48 -50 -51 -52 -53 -56 -57"
removelist4="-10 -11 -13 -18-B -21 -25 -27 -28 -30 -31 -32 -34 -38 -39 -43 -44 -45 -48 -51 -52 -54-B -55 -56 -57"
removelist5="-51 -48 -39 -12 -14 -15 -16 -17 -19 -20 -22 -23 -24 -26 -29 -33 -35 -36 -37 -40 -41 -42 -46 -47 -49 -50 -53"
REMOVELISTS=("${removelist2}" "${removelist3}" "${removelist4}" "${removelist5}")
  # to remove a sample from every set, add it here
    # some files may be too small for StringTie to handle, so list them here as you encounter the errors (although, rnaseq_pipeline.sh should automatically find and remove them for you)
REMOVEALWAYS="-51 -48 -39"


#### Other Options

# BUSCO dataset name -- what taxonomic group should be used as the BUSCO dataset (see instructions BUSCO's website)?
BUSCODATASET="agaricales_odb10"

# Ballgown options
    # to work with more than one set, specify the covariates and confounding variables separately for each set
    # make sure you specify an equal number of setnames and covariates/confounding variables
    # only one covariate should be specified for each set
    # mutliple confounding variables may be specified for each set
  # column name in the phenotype data file to be used as the covariate by Ballgown
covset1="tissue"
covset2="tissue"
covset3="tissue"
covset4="type"
covset5="type"
COVARIATES=("${covset1}" "${covset2}" "${covset3}" "${covset4}" "${covset5}")
  # column names in the phenotype data file to be used as confounding varaibles by Ballgown
adjvars1=("type")
adjvars2=("type")
adjvars3=("type")
adjvars4=("tissue")
adjvars5=("tissue")
ADJVARSETS=("${adjvars1}" "${adjvars2}" "${adjvars3}" "${adjvars4}" "${adjvars5}")
  # for PCA, specify pairs (or singles) of comparisons to make
pcapair1=("batch" "jar")
pcapair2=("type" "tissue")
pcapair3=("position")
PCAPAIRS="${pcapair1};${pcapair2};${pcapair3}"

# NCBI API key (if none, type "")
NCBIAPI=""

# for creating heatmaps and bar plots
  # these lists need to be separated by commas
  # will partial match to the gene, organism, or sample name
genelist1="cytochrome"
genelist2="hydrophobin,psh"
genelist3="mitogen,mapk,map"
genelist4="elongation"
MAPGENELISTS=("${genelist1}" "${genelist2}" "${genelist3}" "${genelist4}")
orglist1="Flammulina"
MAPORGLISTS=("${orglist1}")
samplelist1=("nor,cul")
MAPSAMPLELISTS=("${samplelist1}")
