# from Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095
# modified by Todd Osmundson and Thomas Roehl 2021

#!/usr/bin/env bash
usage() {
    NAME=$(basename $0)
    cat <<EOF
Usage:
  ${NAME} [output_dir]
Wrapper script for HISAT2/StringTie RNA-Seq analysis protocol.
In order to configure the pipeline options (input/output files etc.)
please copy and edit a file rnaseq_pipeline.config.sh which must be
placed in the current (working) directory where this script is being launched.

Output directories "hisat2" and "ballgown" will be created in the 
current working directory or, if provided, in the given <output_dir>
(which will be created if it does not exist).

EOF
}

OUTDIR="./output"
if [[ "$1" ]]; then
 if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit 1
 fi
 OUTDIR=$1
fi

## load variables
if [[ ! -f ./rnaseq_pipeline.config.sh ]]; then
 usage
 echo "Error: configuration file (rnaseq_pipeline.config.sh) missing!"
 exit 1
fi

source ./rnaseq_pipeline.config.sh
WRKDIR=$(pwd -P)
errprog=""
if [[ ! -x $SAMTOOLS ]]; then
    errprog="samtools"
fi
if [[ ! -x $HISAT2 ]]; then
    errprog="hisat2"
fi
if [[ ! -x $STRINGTIE ]]; then
    errprog="stringtie"
fi
if [[ ! -x $BLASTNAPP ]]; then
    errprog="BLASTN"
fi
if [[ ! -x $BLASTXAPP ]]; then
    errprog="BLASTX"
fi
if [[ ! -f $BALLGOWN ]]; then
    errprog="ballgown"
fi
if [[ ! -f $MATCHGENES ]]; then
    errprog="python MSTRNG/gene matching"
fi
if [[ ! -f $NAMEGENES ]]; then
    errprog="R gene name matching"
fi
if [[ ! -f $REMOVERRNA ]]; then
    errprog="R rRNA removal"
fi
if [[ ! -f $FIXUNIPROT ]]; then
    errprog="R UniProt BLAST rewriting"
fi
if [[ ! -f $HEATMAP ]]; then
    errprog="R heat map"
fi
if [[ ! -f $UNIPROTFASTA ]]; then
    errprog="Python3 Uniprot to FASTA conversion"
fi
if [[ ! -f $PANTHERSCORE ]]; then
    errprog="PANTHER HMM scoring"
fi
if [[ ! -f $PANTHERFPKM ]]; then
    errprog="R PANTHER to FPKM mapping"
fi



if [[ "$errprog" ]]; then    
  echo "ERROR: $errprog program not found or not executable; please edit the configuration script."
  exit 1
fi

#determine samtools version
newsamtools=$( ($SAMTOOLS 2>&1) | grep 'Version: 1\.')

set -e
#set -x

if [[ $OUTDIR != "." ]]; then
  mkdir -p $OUTDIR
  cd $OUTDIR
fi

SCRIPTARGS="$@"
ALIGNLOC=./hisat2
BALLGOWNLOC=${OUTDIR}/ballgown
BALLGOWNLOC1=${BALLGOWNLOC}

LOGFILE=./run.log

for d in "$TEMPLOC" "$ALIGNLOC" "$BALLGOWNLOC" ; do
 if [ ! -d $d ]; then
    mkdir -p $d
 fi
done

export PATH="$PATH:${BLASTDIR}"

# main script block
mojo() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $SCRIPTARGS
touch ${ALIGNLOC}/mergelist.txt
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    #sample="${sample%_R1*}"
    sample="${sample%_L001*}"
    echo "Sample name: " ${sample}
    
    # in case of interruption, skip finished files
    if [[ -f ${ALIGNLOC}/${sample}.gtf ]]; then
        echo "${sample} complete in HISAT2 and StringTie. Skipping..."
    else
        
        unpairedBoth=${FASTQLOC}/${unpaired1[$i]},${FASTQLOC}/${unpaired2[$i]},${FASTQLOC}/${unpaired3[$i]}
        stime=`date +"%Y-%m-%d %H:%M:%S"`
        echo "[$stime] Processing sample: $sample"
        echo [$stime] "   * Alignment of reads to genome (HISAT2)"
        $HISAT2 -p $NUMCPUS --dta -x ${GENOMEIDX} \
         -1 ${FASTQLOC}/${reads1[$i]} \
         -2 ${FASTQLOC}/${reads2[$i]} \
         -S ${TEMPLOC}/${sample}.sam 2>${ALIGNLOC}/${sample}.alnstats \
         -U ${unpairedBoth} \
         --seed 12345 --un-gz ${BASEDIR}/unmapped \
         --un-conc ${BASEDIR}/unmapped
         #defaults kept for strandedness, alignment, and scoring options

        echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Alignments conversion (SAMTools)"
        if [[ "$newsamtools" ]]; then
         $SAMTOOLS view -S -b ${TEMPLOC}/${sample}.sam | \
          $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${sample}.bam -
        else
         $SAMTOOLS view -S -b ${TEMPLOC}/${sample}.sam | \
          $SAMTOOLS sort -@ $NUMCPUS - ${ALIGNLOC}/${sample}
        fi
        $SAMTOOLS index ${ALIGNLOC}/${sample}.bam
        #$SAMTOOLS flagstat ${ALIGNLOC}/${sample}.bam
        
        echo "..removing intermediate files"
        rm ${TEMPLOC}/${sample}.sam
        #rm ${TEMPLOC}/${sample}.unsorted.bam

        echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
        $STRINGTIE -p $NUMCPUS -o ${ALIGNLOC}/${sample}.gtf \
         -l ${sample} ${ALIGNLOC}/${sample}.bam \
         -c 1.5
         # -m kept at default of 0.01
    fi
done

## run rest of pipeline on various sample sets
setnames=("all" "300" "5k" "culnor" "priyou")
removelist1="-51 -48 -39"
removelist2="-51 -50 -48 -45 -39 -21"
removelist3="-12 -14 -17 -21 -22 -29 -30 -31 -33 -34 -38 -39 -40 -42 -45 -48 -50 -51 -52 -53 -56 -57"
removelist4="-10 -11 -13 -18-B -21 -25 -27 -28 -30 -31 -32 -34 -38 -39 -43 -44 -45 -48 -51 -52 -54-B -55 -56 -57"
removelist5="-51 -48 -39 -12 -14 -15 -16 -17 -19 -20 -22 -23 -24 -26 -29 -33 -35 -36 -37 -40 -41 -42 -46 -47 -49 -50 -53"
removelists=("${removelist1}" "${removelist2}" "${removelist3}" "${removelist4}" "${removelist5}")

for ((j=0; j<="${#setnames[@]}"-1; j++ )); do
    setname=${setnames[j]}
    removelist=${removelists[j]}
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Beginning analysis of sample set ${setname}"
    
    # make folder and update BALLGOWNLOC
    if [[ ! -d ${BALLGOWNLOC1}/${setname} ]]; then
        mkdir ${BALLGOWNLOC1}/${setname}
    fi
    BALLGOWNLOC=${BALLGOWNLOC1}/${setname}

    ## merge transcript file
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Merge selected transcripts (StringTie)"
    ls -1 ${ALIGNLOC}/*.gtf > ${ALIGNLOC}/${setname}_mergelist.txt
        
    # remove lines from transcript list that match samples in the removelist
    for removestr in ${removelist}; do
        grep -v -- ${removestr} ${ALIGNLOC}/${setname}_mergelist.txt > ${ALIGNLOC}/temp.txt
        mv ${ALIGNLOC}/temp.txt ${ALIGNLOC}/${setname}_mergelist.txt
    done

    $STRINGTIE --merge -p $NUMCPUS \
        -o ${BALLGOWNLOC}/stringtie_merged.gtf ${ALIGNLOC}/${setname}_mergelist.txt

    ## estimate transcript abundance
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Estimate abundance for each sample (StringTie)"
    for ((i=0; i<=${#reads1[@]}-1; i++ )); do
        sample="${reads1[$i]%%.*}"
        dsample="${sample%_L001*}"
        sample="${sample%_L001*}"
        # if the sample name matches anything in the current removelist, skip that sample
        remove="F"
        for removestr in ${removelist}; do
            if [[ "${sample}" == *"${removestr}"* ]]; then
                remove="T"
            fi
        done
        if [ "${remove}" == "F" ]; then
            if [ ! -d ${BALLGOWNLOC}/${dsample} ]; then
                mkdir -p ${BALLGOWNLOC}/${dsample}
            fi
            $STRINGTIE -e -B -p $NUMCPUS -G ${BALLGOWNLOC}/stringtie_merged.gtf \
            -o ${BALLGOWNLOC}/${dsample}/${dsample}.gtf ${ALIGNLOC}/${sample}.bam
        fi
    done

    # create transcript FASTA file
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Generate fasta transcript file from .gtf"
    ## note: hard to find genome using index location, so whole path is typed out
    $GFFREAD -w ${BALLGOWNLOC}/fv_transcriptome.fa -g /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/genome_assemblies_all_files/ncbi-genomes-2021-11-01/GCA_011800155.1_ASM1180015v1/GCA_011800155.1_ASM1180015v1_genomic.fna ${BALLGOWNLOC}/stringtie_merged.gtf

    ### BLAST against ribosomal database

    # make rRNA database
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Create rRNA database (makeblastdb)"

    $BLASTDIR/makeblastdb -in ${rrna_file} -out ${RRNADIR}/ncbi_rrna_fungi -dbtype nucl

    # BLAST rRNA database
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Query rRNA database with transcroptome (BLASTN)"

    $BLASTNAPP -query ${BALLGOWNLOC}/fv_transcriptome.fa \
        -db ${RRNADIR}/ncbi_rrna_fungi \
        -out ${RRNADIR}/blastn_results.csv \
        -evalue 1e-6 -num_threads ${NUMCPUS} \
        -num_alignments 1 -outfmt "6 qseqid stitle sacc evalue pident bitscore length qstart qend sstart send"

    ## remove matched transcripts from ctab files
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Removing rRNA..."

    # save a copy of the original tables in a new folder
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Save a copy of original expression tables"

    for direc in ${BALLGOWNLOC}/ROEHL-FV-*/; do
        for table in ${direc}*.ctab; do
            if [[ ! -d ${direc}original_tables ]]; then
                mkdir ${direc}original_tables
            fi
            basename=${table##*/}
            # restore original files -- do on every rerun
            if [[ -f ${direc}original_tables/${basename} ]]; then
                cp ${direc}/original_tables/${basename} ${table}
            fi
            # save original files -- do on only first run
            if [[ ! -f ${direc}original_tables/${basename} ]]; then
                cp ${table} ${direc}original_tables/${basename}
            fi
        done
    done

    # use R script to remove rRNA and overwrite the .ctab files (leaving originals in subdirectories)
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Remove rRNA from expression tables (R)"

    if [ ! -d ${OUTDIR}/R_logs/ ]; then
        mkdir ${OUTDIR}/R_logs
    fi

    Rscript ${REMOVERRNA} ${RRNADIR}/blastn_results.csv ${BALLGOWNLOC} ${OUTDIR}


    ### busco completeness test
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Testing Busco completeness..."

    busconame="busco_output_fv"

    if [[ -d ${BALLGOWNLOC}/${busconame} ]]; then
        busconame="busco_output_fv_"`date +"%Y-%m-%d_%H-%M-%S"`
    fi

    docker run -u $(id -u) -v ${BALLGOWNLOC}/:/busco_wd/ ezlabgva/busco:v4.1.4_cv1 busco \
        --mode transcriptome \
        --in /busco_wd/fv_transcriptome.fa \
        --out ${busconame} \
        --lineage_dataset agaricales_odb10
    # note: -v [specify complete file path where you want busco to live]:/busco_wd/
    # note: the working directory (specified in -v) must be /busco_wd/
    # note: all file paths (except in -v) must be relative to the image root directory. In this case, -v maps the image root (/busco_wd/) to ${BALLGOWNLOC}/ . As a result, the transcriptome file must be referenced as /busco_wd/fv_transcriptome.fa because it is found in ${BALLGOWNLOC}/fv_transcriptome.fa
    # note: --out file must have a name that has never been used for a busco run. Ex: busco was previously run on the same machine with "--out busco_output" so this run must use a new name, "--out busco_output_fv"
    # note: if script must be run more than once, delete the output directories before running OR make sure to run the busconame if statement above

    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Generate the DE tables (Ballgown)"

    if [[ ! -d ${BALLGOWNLOC}/bg_output ]]; then
        mkdir ${BALLGOWNLOC}/bg_output
    fi
    if [[ ! -d ${BALLGOWNLOC}/bg_output/lists ]]; then
        mkdir ${BALLGOWNLOC}/bg_output/lists
    fi
    if [[ ! -d ${BALLGOWNLOC}/bg_output/seq_lists ]]; then
        mkdir ${BALLGOWNLOC}/bg_output/seq_lists
    fi
    if [[ ! -d ${BALLGOWNLOC}/bg_output/PCA ]]; then
        mkdir ${BALLGOWNLOC}/bg_output/PCA
    fi
    if [[ ! -d ${BALLGOWNLOC}/bg_output/named ]]; then
        mkdir ${BALLGOWNLOC}/bg_output/named
    fi

    Rscript ${BALLGOWN} ${PHENODATA} ${BALLGOWNLOC}


    ## BLAST DEGs

    # create BLAST database
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Create gene database (makeblastdb)"
    $BLASTDIR/makeblastdb -in ${uniprot_file} -out ${UNIPROTDIR}/uniprot_agaricales -dbtype prot

    for deglist in ${BALLGOWNLOC}/bg_output/lists/*.txt; do
        currentDEGList=${deglist##*/}
        currentDEGs=${currentDEGList%.txt}
        echo "Current DEG list: ${currentDEGs}"

        # biopython: match MSTRNG numbers from DEGs and .fa sequences
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Match MSTRNGs to genes for DEGs (Biopython)"
        filename=${BALLGOWNLOC}/bg_output/seq_lists/${currentDEGs}.fa
        if [[ -f ${filename} ]]; then
            rm ${filename}
        fi

        # 3 arguments:
        # output file name (fasta format)
        # fasta transcriptome with MSTRNG/sequence pairs
        # .txt list of DEGs separated by line
        python3 ${MATCHGENES} ${filename} ${BALLGOWNLOC}/fv_transcriptome.fa ${deglist}

        # search BLAST database
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Query BLAST database for gene names (BLASTX)"
        
        $BLASTXAPP -query ${BALLGOWNLOC}/bg_output/seq_lists/${currentDEGs}.fa \
            -db ${UNIPROTDIR}/uniprot_agaricales \
            -out ${UNIPROTDIR}/${setname}_${currentDEGs}_blastx_results.csv \
            -evalue 1e-3 -num_threads ${NUMCPUS} \
            -num_alignments 1 -outfmt "6 qseqid stitle sacc evalue pident bitscore length qstart qend sstart send"
        
        # extract infromation from UniProt stitle
        Rscript ${FIXUNIPROT} ${UNIPROTDIR}/${setname}_${currentDEGs}_blastx_results.csv ${UNIPROTDIR}/${setname}_${currentDEGs}_blastx_results_readable.csv ${UNIPROTDIR}
        
        # match protein names to output files
        for bgresults in ${BALLGOWNLOC}/bg_output/${currentDEGs}*csv; do
            resultname=${bgresults##*/}
            Rscript ${NAMEGENES} ${UNIPROTDIR}/${setname}_${currentDEGs}_blastx_results_readable.csv ${bgresults} ${BALLGOWNLOC}/bg_output/named/${resultname} ${OUTDIR}
        done
    done


    # BLAST transcriptome for PANTHER reference
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Query BLAST database with transcriptome (BLASTX)"
    
    $BLASTXAPP -query ${BALLGOWNLOC}/fv_transcriptome.fa \
        -db ${UNIPROTDIR}/uniprot_agaricales \
        -out ${UNIPROTDIR}/${setname}_tome_blastx_results.csv \
        -evalue 1e-3 -num_threads ${NUMCPUS} \
        -num_alignments 1 -outfmt "6 qseqid stitle sacc evalue pident bitscore length qstart qend sstart send"
    
    # extract infromation from UniProt stitle
    Rscript ${FIXUNIPROT} ${UNIPROTDIR}/${setname}_tome_blastx_results.csv ${UNIPROTDIR}/${setname}_tome_blastx_results_readable.csv ${UNIPROTDIR}


     ## generate heat maps
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Creating Heat Maps"

    if [[ ! -d ${BALLGOWNLOC}/bg_output/heatmaps ]]; then
        mkdir ${BALLGOWNLOC}/bg_output/heatmaps
    fi

    for resultfile in ${BALLGOWNLOC}/bg_output/named/*_expr.csv; do
        basename=${resultfile##*/}
        base=${basename%.csv}
       
        # map all genes
        Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${BALLGOWNLOC}/bg_output/heatmaps/${base}_heatmap_all.pdf "all" "all" "all" ${BALLGOWNLOC}/bg_output/
       
        # map specific genes
        #genelist=("cytochrome" "hydrophobin,psh" "mitogen,mapk,map" "elongation")
        #for ((k=0; k<=${#genelist[@]}-1; k++)); do
               # Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${BALLGOWNLOC}/bg_output/heatmaps/${base}_heatmap_genes_${k}.pdf ${genelist[k]} "all" "all" ${BALLGOWNLOC}/bg_output/
        #done
       
        # map genes from specific organisms
        #orglist=("Flammulina")
        #for ((k=0; k<=${#orglist[@]}-1; k++)); do
                #Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${BALLGOWNLOC}/bg_output/heatmaps/${base}_heatmap_orgs_${k}.pdf "all" "all" ${orglist[k]} ${BALLGOWNLOC}/bg_output/
        #done
       
        # map specific samples
        #samplelist=("nor,cul")
        #for ((k=0; k<=${#samplelist[@]}-1; k++)); do
             #   Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${BALLGOWNLOC}/bg_output/heatmaps/${base}_heatmap_samples_${k}.pdf "all" ${samplelist[k]} "all" ${BALLGOWNLOC}/bg_output/
        #done
    done
    
    
    ## GO analysis
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Generating tables for GO analysis..."

    if [[ ! -d ${BALLGOWNLOC}/panther ]]; then
        mkdir ${BALLGOWNLOC}/panther
    fi
    if [[ ! -d ${BALLGOWNLOC}/panther/output ]]; then
        mkdir ${BALLGOWNLOC}/panther/output
    fi

    # for each DEG list and reference list...
    
    # PantherScore won't work unless the working directory is the directory that contains PantherScore and the program's lib folder
    cd /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/panther/pantherScore2.2

    # for each named protein, get the NCBI Protein database sequence and add it to a FASTA file
    python3 ${UNIPROTFASTA} ${BALLGOWNLOC}/panther/tome_ncbi_list.csv ${BALLGOWNLOC}/panther/tome_fasta.fa "" ${UNIPROTDIR}/${setname}_tome_blastx_results_readable.csv
    # score FASTA file against PANTHER HMMs to map to PANTHER IDs
    perl ${PANTHERSCORE} -l ${PANTHERLIBDIR} -D B -V -i ${BALLGOWNLOC}/panther/tome_fasta.fa -o ${BALLGOWNLOC}/panther/output/tome_panther_mapping.txt -n -s

    for deglist in ${BALLGOWNLOC}/bg_output/named/*_expr.csv; do
        base=${deglist##*/}
        basename=${base%.csv}
        echo "Current list: ${basename}"

        # for each named protein, get the NCBI Protein database sequence and add it to a FASTA file
        python3 ${UNIPROTFASTA} ${BALLGOWNLOC}/panther/${basename}_ncbi_list.csv ${BALLGOWNLOC}/panther/${basename}_fasta.fa "" ${deglist}
        # score FASTA file against PANTHER HMMs to map to PANTHER IDs
        perl ${PANTHERSCORE} -l ${PANTHERLIBDIR} -D B -V -i ${BALLGOWNLOC}/panther/${basename}_fasta.fa -o ${BALLGOWNLOC}/panther/${basename}_panther_mapping.csv -n -s
            # output is tab-delimited: sequence ID, panther acc, panther family/subfamily, HMM e-value, HMM bitscore, alignment range
        
        # map PANTHER IDs to fpkm tables
        Rscript ${PANTHERFPKM} ${deglist} ${BALLGOWNLOC}/panther/${basename}_panther_mapping.csv ${PHENODATA} ${BALLGOWNLOC}/panther/output/${basename} ${BALLGOWNLOC}
    done
    
    cd ${WRKDIR}
    
done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

mojo 2>&1 | tee $LOGFILE
