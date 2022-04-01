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


##### NOTE
##### The following pattern is used to check the output log to see if specific events were completed and to skip the completed sections (if the user has opted to resume):

#skip="N"
#chklog "completion_text"
#if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
#    skip="Y"
#fi
#if [[ ${skip} == "N" ]]; then
#    do_normal_function
#    echo "completion_text"
#fi

# to remove this feature, find and delete all instances of this pattern



#### METHODS DECLARATIONS



### convenience method for emptying folders
    # pass folder you want to empty as a parameter when calling this function
    # specified folder will not be deleted, but all its contents will be removed
emptydir() {
    folder=$1
    rm -rf ${folder}/*
}



### convenience method for checking the log
    # pass paramater of the entry you want to check
    # stores result as "T" or "F" in variable chkresult
chklog() {
    nmatch=$(grep -c $1 ${LOGFILE})
    chkresult="F"
    if [[ nmatch > 0 ]]; then
        chkresult="T"
    fi
}



### check programs required for each function block
    # when used, add an argument for block number
chkprog() {

    errprog=""
    
    # programs for first block (initqc)
    if [[ $1 == 1 ]]; then

        if [[ ! -x $FASTQC ]]; then
            errprog="FastQC"
        fi
        if [[ ! -x $FQTRIM ]]; then
            errprog="FQTrim"
        fi
        if [[ ! -f $TRIMMOMATIC ]]; then
            errprog="Trimmomatic"
        fi
        if [[ ! -x $INTERLEAVE ]]; then
            errprog="Interleave Pairs"
        fi
    fi
    
    # programs for second block (initproc)
    if [[ $1 == 2 ]]; then
        if [[ ! -x $SAMTOOLS ]]; then
            errprog="samtools"
        else
            #determine samtools version
            newsamtools=$( ($SAMTOOLS 2>&1) | grep 'Version: 1\.')
        fi
        if [[ ! -x $HISAT2 ]]; then
            errprog="hisat2"
        fi
        if [[ ! -x $STRINGTIE ]]; then
            errprog="stringtie"
        fi
    fi
    
    # programs for third block (mojo)
    if [[ $1 == 3 ]]; then

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
    fi
    
    if [[ "$errprog" ]]; then    
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ERROR: $errprog program not found or not executable; please edit the configuration script."
        exit 1
    else
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> All programs found."
    fi

}



### initial QC block
initqc() {
    
    
    ## FASTQC analysis of raw files
    skip="N"
    chklog "raw_fastqc_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Beginning initial quality control."

        # add directory for fastqc output
        if [[ ! -d ${output}/fastqc ]]; then
            mkdir ${output}/fastqc
        fi
        # add subdirectory for raw file analysis
        if [[ ! -d ${output}/fastqc/raw ]]; then
            mkdir ${output}/fastqc/raw
        fi

        ## generate FastQC reports from raw data
        for seqfile in ${DATADIR}/*.fastq.gz; do
            skip="N"
            chklog "${output}/fastqc/raw/${seqfile##*/}"
            if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                ${fastqc_app} -t ${NCPUS} -o ${output}/fastqc/raw ${seqfile}
                echo ${output}/fastqc/raw/${seqfile##*/}
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> raw_fastqc_complete."
    fi
    
    
    ## concatenate files from multiple runs if needed
    newdatadir=${DATADIR}
    if [[ ${#SEQBATCHES[@]} > 1 ]]; then
        newdatadir=${output}/concatenate
    fi
    
    skip="N"
    chklog "concatenation_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        if [[ ${#SEQBATCHES[@]} > 1 ]]; then

            # create new folder for concatenation
            if [[ ! -d ${newdatadir} ]]; then
                mkdir ${newdatadir}
            fi

            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Concatenating files by sample (preserving forward/reverse)..."

            for seqfile in ${DATADIR}/*${SEQBATCHES[0]}*.fastq.gz; do
                # remove 0th seq batch pattern to create the base name for the concatenated file
                nameending=${seqfile##*/}
                basename=${nameending/"${SEQBATCHES[0]}"/}

                skip="N"
                chklog "${newdatadir}/${basename}"
                if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                    skip="Y"
                fi
                if [[ ${skip} == "N" ]]; then
                    for (i = 0; i < ${#SEQBATCHES[@]}; i++); do
                        # find 0th seq batch pattern and replace with ith seq batch pattern to retrieve next file in the sample, then append to concatenation file
                        nextfile=${seqfile/"${SEQBATCHES[0]}"/"${SEQBATCHES[i]}"}
                        cat ${nextfile} >> ${newdatadir}/${basename}
                    done
                    
                    echo ${newdatadir}/${basename}
                fi
            done
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> concatenation_complete."
        fi
    fi
    
    
    ## Adapter and 3' quality trimming using trimmomatic
    skip="N"
    chklog "trimmomatic_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Adapter and 3' quality trimming using Trimmomatic..."

        # add directory for trimmomatic output
        if [[ ! -d ${output}/trimmomatic ]]; then
            mkdir ${output}/trimmomatic
        fi

        # I don't know why the following line was included -- test to see if this section works without copying the files into the working directory
        #cp ${trimmomatic_adapter_dir}/TruSeq3-PE-2.fa ./TruSeq3-PE-2.fa
        for filename in ${newdatadir}/*${FWDREV[0]}*.fastq.gz; do
            nameending=${filename##*/}
            name=${nameending%.fastq.gz}
            basename=${name/"${FWDREV[0]}"/}

            skip="N"
            chklog "${output}/trimmomatic/${basename}_2U.fq"
            if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                 java -jar ${TRIMMOMATIC} PE -threads ${NCPUS} -trimlog ${LOGLOC}/trimmomatic_trimlog.txt -basein ${filename} -baseout ${output}/trimmomatic/${basename}.fq ILLUMINACLIP:${TRIMMOMATICADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
                 
                 echo ${output}/trimmomatic/${basename}_2U.fq
             fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> trimmomatic_complete"
    fi
    
    
    ## Remove low-complexity sequences using fqtrim
    skip="N"
    chklog "fqtrim_1p_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Removing low-complexity sequences and poly-A/poly-T sequences using fqtrim..."

        # make fqtrim1 folder
        if [[ ! -d ${output}/fqtrim1 ]]; then
            mkdir ${output}/fqtrim1
        fi

        # trimmomatic outputs paired sequences as ./*_1P.fq and ./*_2P.fq
        for file in ${output}/trimmomatic/*_1P.fq; do
            readpair=$(echo $file,${file/_1P/_2P} | tr -d '\n')

            skip="N"
            chklog "${readpair}"
            if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
            
                ${FQTRIM} -l 30 -p ${NCPUS} -D -o fqtrimmed.fq --outdir ${output}/fqtrim1 -r ${LOGLOC}/fqtrimlog1p.txt ${readpair}
                
                echo ${readpair}
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> fqtrim_1p_complete"
    fi
    
    
    ## incorporate unpaired sequences, if option is Y
    skip="N"
    chklog "fqtrim_1u_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        if [[ ${USEUNPAIRED} == "Y" ]]; then
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Running fqtrim on unpaired reads..."

            # trimmomatic outputs unpaired sequences as ./*_1U.fq and ./*_2U.fq
            for file in ${output}/trimmomatic/*U.fq; do
                skip="N"
                chklog "${output}/fqtrim1/${file}"
                if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                    skip="Y"
                fi
                if [[ ${skip} == "N" ]]; then

                    ${FQTRIM} -l 30 -p ${NCPUS} -D -o fqtrimmed.fq --outdir ${output}/fqtrim1 -r ${LOGLOC}/fqtrimlog1u.txt ${file}

                    echo ${output}/fqtrim1/${file}
                fi
            done
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> fqtrim_1u_complete"
        fi
    fi
    
    
    ## run FastQC on trimmed reads
    skip="N"
    chklog "trim1_fastqc_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Conducting FastQC quality check of processed reads..."

        # create folder for trimmed FastQC
        if [[ ! -d ${output}/fastqc/trim1 ]]; then
            mkdir ${output}/fastqc/trim1
        fi

        for seqfile in ${output}/fqtrim1/*.fqtrimmed.fq; do
            skip="N"
            chklog "${output}/fastqc/trim1/${seqfile##*/}"
            if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
            
                ${FASTQC} -t ${NCPUS} -o ${output}/fastqc/trim1 ${seqfile}
                echo ${output}/fastqc/trim1/${seqfile##*/}
                
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> trim1_fastqc_complete"
    fi
    
    
    ## rerun fqtrim on paired sequences, but ignore pairs
        # to preserve pair order, fqtrim will leave single-nucleotide sequences in the files
        # rerunning fqtrim is required to remove these junk sequences
    skip="N"
    chklog "fqtrim_2p_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Rerunning fqtrim on paired sequences..."

        # make new fqtrim output folder
        if [[ ! -d ${output}/fqtrim2 ]]; then
            mkdir ${output}/fqtrim2
        fi

        # run fqtrim on paired sequences
        for seqfile in ${output}/fqtrim1/*P.fqtrimmed.fq; do
            skip="N"
            chklog "${output}/fqtrim2/${seqfile##*/}"
            if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then

                ${FQTRIM} -l 30 -p ${NCPUS} -D -o 2.fq --outdir ${output}/fqtrim2 -r ${LOGLOC}/fqtrimlog2p.txt ${seqfile}
                echo ${output}/fqtrim2/${seqfile##*/}
                
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> fqtrim_2p_complete"
    fi
    
    
    ## run FastQC on retrimmed reads
    skip="N"
    chklog "trim2_fastqc_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Conducting FastQC quality check of fqtrim reruns..."

        # create folder for rerun FastQC
        if [[ ! -d ${output}/fastqc/trim2 ]]; then
            mkdir ${output}/fastqc/trim2
        fi

        for seqfile in ${output}/fqtrim2/*.fq; do
            skip="N"
            chklog "${output}/fastqc/trim2/${seqfile##*/}"
            if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
            
                ${FASTQC} -t ${NCPUS} -o ${output}/fastqc/trim2 ${seqfile}
                echo ${output}/fastqc/trim2/${seqfile##*/}
                
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> trim2_fastqc_complete"
    fi
    
    
    ## reorder paired sequences
        # paired sequences need to be in the same order in each pair file for the next step
    skip="N"
    chklog "interleave_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Performing sequence pair matching..."

        # create folder for interleave step
        if [[ ! -d ${output}/interleave ]]; then
            mkdir ${output}/interleave
        fi

        # copy over unpaired reads from fqtrim1
        for seqfile in ${output}/fqtrim1/*U.fqtrimmed.fq; do
            basename=${seqfile##*/}
            echo ${basename}

            cp ${seqfile} ${output}/interleave/${basename}
        done

        # perform interleave
        for file in ${output}/fqtrim2/*_1P.fqtrimmed.2.fq; do
            reverse=$(echo ${file/_1P/_2P} | tr -d '\n')
            basename=${file##*/}
            base=${output}/interleave/${basename%_*}
            
            skip="N"
            chklog "${output}/interleave/${base}_out_unpaired.fastq"
            if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then

                python3 ${INTERLEAVE} ${file} ${reverse} ${base}
                # output files are ./*_out_pairs_fwd.fastq ./*_out_pairs_rev.fastq and ./*_out_unpaired.fastq
                echo ${output}/interleave/${base}_out_unpaired.fastq
                
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> interleave_complete"
    fi
    
    
    ## end of block file management
    skip="N"
    chklog "initqc_file_management_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> File Management..."

        # copy final files to local input folder
            # in case the process was interrupted after moving /interleave but before reaching completing file management, avoid errors by first checking to see if the /interleave folder is still in the output folder
        if [[ -d ${output}/interleave ]]; then
            for file in ${output}/interleave/*; do
                basename=${file##*/}
                echo ${basename}

                cp ${file} ${input}/${basename}
            done
        fi

        # move contents of output folder to destination folder
        if [[ ! -d ${DESTDIR}/initqc ]]; then
            mkdir ${DESTDIR}/initqc
        fi

        mv ${output}/* ${DESTDIR}/initqc
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> initqc_file_management_complete"
    fi
    
    
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> initqc_complete"
}



### initial sample processing
initproc() {

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
    
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> initproc_complete"
}



### analysis and visualization
mojo() {

    setname=$1
    removelist=$2
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
                cp ${direc}/original_tables/${basename} ${table}
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


# should be done automatically by adding complete list of transcripts to ballgown/output/lists -- see rnaseq_ballgown
#    # BLAST transcriptome for PANTHER reference
#    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Query BLAST database with transcriptome (BLASTX)"
#    
#    $BLASTXAPP -query ${BALLGOWNLOC}/fv_transcriptome.fa \
#        -db ${UNIPROTDIR}/uniprot_agaricales \
#        -out ${UNIPROTDIR}/${setname}_tome_blastx_results.csv \
#        -evalue 1e-3 -num_threads ${NUMCPUS} \
#        -num_alignments 1 -outfmt "6 qseqid stitle sacc evalue pident bitscore length qstart qend sstart send"
#    
#    # extract infromation from UniProt stitle
#    Rscript ${FIXUNIPROT} ${UNIPROTDIR}/${setname}_tome_blastx_results.csv ${UNIPROTDIR}/${setname}_tome_blastx_results_readable.csv ${UNIPROTDIR}


     ## generate heat maps
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Creating Heat Maps"

    if [[ ! -d ${BALLGOWNLOC}/bg_output/heatmaps ]]; then
        mkdir ${BALLGOWNLOC}/bg_output/heatmaps
    fi

    for resultfile in ${BALLGOWNLOC}/bg_output/named/*_expr.csv; do
        basename=${resultfile##*/}
        base=${basename%.csv}
       
        # map all genes
        Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${BALLGOWNLOC}/bg_output/heatmaps/${base}_heatmap_all.pdf "all" "all" "all" "F" ${BALLGOWNLOC}/bg_output/
       
        # map specific genes
        #genelist=("cytochrome" "hydrophobin,psh" "mitogen,mapk,map" "elongation")
        #for ((k=0; k<=${#genelist[@]}-1; k++)); do
               # Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${BALLGOWNLOC}/bg_output/heatmaps/${base}_heatmap_genes_${k}.pdf ${genelist[k]} "all" "all" "F" ${BALLGOWNLOC}/bg_output/
        #done
       
        # map genes from specific organisms
        #orglist=("Flammulina")
        #for ((k=0; k<=${#orglist[@]}-1; k++)); do
                #Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${BALLGOWNLOC}/bg_output/heatmaps/${base}_heatmap_orgs_${k}.pdf "all" "all" ${orglist[k]} "F" ${BALLGOWNLOC}/bg_output/
        #done
       
        # map specific samples
        #samplelist=("nor,cul")
        #for ((k=0; k<=${#samplelist[@]}-1; k++)); do
             #   Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${BALLGOWNLOC}/bg_output/heatmaps/${base}_heatmap_samples_${k}.pdf "all" ${samplelist[k]} "all" "F" ${BALLGOWNLOC}/bg_output/
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
    cd ${PANTHERSCORE%/*}

    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Converting UniProtKB IDs to FASTA..."
    # for each named protein, get the NCBI Protein database sequence and add it to a FASTA file
    python3 ${UNIPROTFASTA} ${BALLGOWNLOC}/panther/tome_ncbi_list.csv ${BALLGOWNLOC}/panther/tome_fasta.fa "" ${UNIPROTDIR}/${setname}_tome_blastx_results_readable.csv
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Running PantherScore..."
    # score FASTA file against PANTHER HMMs to map to PANTHER IDs
    perl ${PANTHERSCORE} -l ${PANTHERLIBDIR} -D B -i ${BALLGOWNLOC}/panther/tome_fasta.fa -o ${BALLGOWNLOC}/panther/output/tome_panther_mapping.txt -n -s

    for deglist in ${BALLGOWNLOC}/bg_output/named/*_expr.csv; do
        base=${deglist##*/}
        basename=${base%.csv}
        echo "Current list: ${basename}"

        # for each named protein, get the NCBI Protein database sequence and add it to a FASTA file
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Converting UniProtKB IDs to FASTA..."
        python3 ${UNIPROTFASTA} ${BALLGOWNLOC}/panther/${basename}_ncbi_list.csv ${BALLGOWNLOC}/panther/${basename}_fasta.fa "" ${deglist}
        # score FASTA file against PANTHER HMMs to map to PANTHER IDs
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Running PantherScore..."
        perl ${PANTHERSCORE} -l ${PANTHERLIBDIR} -D B -i ${BALLGOWNLOC}/panther/${basename}_fasta.fa -o ${BALLGOWNLOC}/panther/output/${basename}_panther_mapping.csv -n -s
            # output is tab-delimited: sequence ID, panther acc, panther family/subfamily, HMM e-value, HMM bitscore, alignment range
        
        # map PANTHER IDs to fpkm tables
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Matching PANTHER IDs to FPKM tables..."
        if [ ! -f ${BALLGOWNLOC}/panther/output/${basename} ]; then
            mkdir ${BALLGOWNLOC}/panther/output/${basename}
        fi
        Rscript ${PANTHERFPKM} ${deglist} ${BALLGOWNLOC}/panther/output/${basename}_panther_mapping.csv ${BALLGOWNLOC}/panther/${basename}_ncbi_list.csv ${PHENODATA} "ROEHL,myc,sti,pil,gil,you,pri,cul,nor" ${basename}_panther_mapping ${BALLGOWNLOC}/panther/output/${basename}
    done
    
    cd ${WRKDIR}
    
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_mojo_complete"
}



#### PIPELINE PROCESS



# initial setup

## load variables
SCRIPTARGS="$@"
if [[ ! -f ${SCRIPTARGS[1]} ]]; then
    usage
    echo "Error: configuration file (rnaseq_pipeline.config.sh) missing or not included as argument in command line!"
    exit 1
else
    source ${SCRIPTARGS[1]}
fi

set -e

if [[ $OUTDIR != "." ]]; then
  mkdir -p $OUTDIR
  cd $OUTDIR
fi

ALIGNLOC=./hisat2
BALLGOWNLOC=${OUTDIR}/ballgown
BALLGOWNLOC1=${BALLGOWNLOC}

LOGFILE=${LOGLOC}/rnaseq_pipeline_run.log

for d in "$TEMPLOC" "$ALIGNLOC" "$BALLGOWNLOC" ; do
 if [ ! -d $d ]; then
    mkdir -p $d
 fi
done

export PATH="$PATH:${BLASTDIR}:${TRIMMOMATICADAPTERS}"




# check for previous log file - if found, ask user whether to start anew and overwrite or to skip completed files
resume="N"
if [[ -f ${LOGFILE} ]]; then
    chklog "PIPELINE-COMPLETE"
    if [[ ${chkresult} == "F" ]]; then
        echo "Previous log file found. You can either resume from the log file or rerun the entire pipeline. Do you want to resume from the log? (Y/N)"
        read resume
    fi
fi
if [[ ${resume} == "Y" ]] || [[ ${resume} == "y" ]] || [[ ${resume} == "Yes" ]] || [[ ${resume} == "YES" ]] || [[ ${resume} == "yes" ]]; then
    # note: log file can be used as record of previous run because we append to the file so the previous record remains intact
    resume="Y"
else
    resume="N"
    # if not resuming, overwrite log file
    echo "rnaseq_pipeline.sh log begins" > $LOGFILE
fi



# create local input and output folders
input=${WRKDIR}/input
output=${WRKDIR}/output
if [[ ! -d ${input} ]]; then
    mkdir ${input}
fi
if [[ ! -d ${output} ]]; then
    mkdir ${output}
fi



echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $SCRIPTARGS

# check that required files are present for block 1
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Checking programs for part 1..."
chkprog 1

# (if skip) check if initial QC was previously completed - if so, skip initial QC
skip="N"
chklog "initqc_complete"
if [[ ${resume} == "Y" && ${chkresult} == "Y" ]]; then
    skip="Y"
fi
if [[ ${skip} == "N" ]]; then

    # do initial QC
    initqc 2>&1 | tee -a $LOGFILE
    
fi



# check that required files are present for block 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Checking programs for part 2..."
chkprog 2
# (if skip) check if initial analysis was previously completed - if so, skip inital analysis
# do initial analysis (if skip, skip completed files)
initproc 2>&1 | tee -a $LOGFILE
# empty input folder
# copy needed files to input folder
# move output files/directories to destination directory



# check that required files are present for block 3
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Checking programs for part 3..."
chkprog 3
# for each data subset...
    # (if skip) check if initial analysis was previously completed - if so, skip inital analysis
    # run analysis pipeline (if skip, skip completed files)
## run rest of pipeline on various sample sets
for ((j=0; j<="${#SETNAMES[@]}"-1; j++ )); do
    mojo ${setnames[j]} ${removelists[j]} 2>&1 | tee -a $LOGFILE
done
# move output files/directories to destination directory

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> PIPELINE-COMPLETE."
