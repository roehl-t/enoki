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
#
#skip="N"
#chklog "completion_text"
#if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
#    skip="Y"
#fi
#if [[ ${skip} == "N" ]]; then
#    do_normal_function
#    echo "completion_text"
#fi
#
#



#### METHODS DECLARATIONS



### convenience method for emptying folders
    # pass folder you want to empty as a parameter when calling this function
    # specified folder will not be deleted, but all its contents will be removed
    # if statement to protect files from deletion outside working directory
emptydir() {
    folder=$1
    if [[ ${WRKDIR} == ${folder}* ]]; then
        rm -rf ${folder}/*
    else
        echo "bad rm location: ${file}"
        exit 1
    fi
}



### convenience method for checking the log
    # pass paramater of the entry you want to check
    # stores result as "T" or "F" in variable chkresult
chklog() {
    set +e
    nmatch=$(grep -c $1 ${LOGFILE})
    chkresult="F"
    if (( nmatch > 0 )); then
        chkresult="T"
    fi
    set -e
}



### check programs required for each function block
    # when used, add an argument for block number
chkprog() {
    
    blockno=$1
    errprog=""
    
    # programs for first block (initqc)
    if [[ ${blockno} == 1 ]]; then

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
    if [[ ${blockno} == 2 ]]; then
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
    
    # programs for third block (listprep)
    if [[ ${blockno} == 3 ]]; then

        if [[ ! -x $STRINGTIE ]]; then
            errprog="stringtie"
        fi
        if [[ ! -x $BLASTNAPP ]]; then
            errprog="BLASTN"
        fi
        if [[ ! -f $REMOVERRNA ]]; then
            errprog="R rRNA removal"
        fi
        if [[ ! -f $TRANSCRIPTCOUNT ]]; then
            errprog="R transcript counter"
        fi
    fi
    
    # programs for fourth block (mojo)
    if [[ ${blockno} == 4 ]]; then

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
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
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
            chklog "${output}/fastqc/raw/${seqfile##*/}_FastQC_done"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                ${FASTQC} -t ${NUMCPUS} -o ${output}/fastqc/raw ${seqfile}
                echo "${output}/fastqc/raw/${seqfile##*/}_FastQC_done"
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> raw_fastqc_complete."
    fi
    
    
    ## concatenate files from multiple runs if needed
    newdatadir=${DATADIR}
    needsconcat="F"
    if (( ${#SEQBATCHES[@]} > 1 )); then
        newdatadir=${output}/concatenate
        needsconcat="T"
    fi
    
    skip="N"
    chklog "concatenation_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" && ${needsconcat} == "T" ]]; then

        # create new folder for concatenation
        if [[ ! -d ${newdatadir} ]]; then
            mkdir ${newdatadir}
        fi

        if (( ${#SEQBATCHES[@]} > 1 )); then

            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Concatenating files by sample (preserving forward/reverse)..."

            for seqfile in ${DATADIR}/*${SEQBATCHES[0]}*.fastq.gz; do
                # remove 0th seq batch pattern to create the base name for the concatenated file
                nameending=${seqfile##*/}
                basename=${nameending/"${SEQBATCHES[0]}"/}

                skip="N"
                chklog "${newdatadir}/${basename}_cat_done"
                if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                    skip="Y"
                fi
                if [[ ${skip} == "N" ]]; then
                    for ((i=0; i<${#SEQBATCHES[@]}; i++)); do
                        # find 0th seq batch pattern and replace with ith seq batch pattern to retrieve next file in the sample, then append to concatenation file
                        nextfile=${seqfile/"${SEQBATCHES[0]}"/"${SEQBATCHES[i]}"}
                        cat ${nextfile} >> ${newdatadir}/${basename}
                    done
                    
                    echo "${newdatadir}/${basename}_cat_done"
                fi
            done
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> concatenation_complete."
        else
            # remove seqbatch ID from file name
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Renaming files..."
            for seqfile in ${DATADIR}/*.fastq.gz; do
                nameending=${seqfile##*/}
                basename=${nameending/"${SEQBATCHES[0]}"/}
                cp ${seqfile} ${newdatadir}/${basename}
            done
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> concatenation_complete"
        fi
    fi
    
    
    ## Adapter and 3' quality trimming using trimmomatic
    skip="N"
    chklog "trimmomatic_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Adapter and 3' quality trimming using Trimmomatic..."

        # add directory for trimmomatic output
        if [[ ! -d ${output}/trimmomatic ]]; then
            mkdir ${output}/trimmomatic
        fi

        for fwdname in ${newdatadir}/*${FWDREV[0]}*.fastq.gz; do
            nameending=${fwdname##*/}
            name=${nameending%.fastq.gz}
            basename=${name/"${FWDREV[0]}"/}
            revname=${fwdname/"${FWDREV[0]}"/"${FWDREV[1]}"}

            skip="N"
            chklog "trimmomatic_${basename}_complete"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                 java -jar ${TRIMMOMATIC} PE -threads ${NUMCPUS} -trimlog ${LOGLOC}/trimmomatic_trimlog.txt ${fwdname} ${revname} -baseout ${output}/trimmomatic/${basename}.fq ILLUMINACLIP:${TRIMMOMATICADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
                 
                 echo "trimmomatic_${basename}_complete"
             fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> trimmomatic_complete"
    fi
    
    ## file management: move concatenated sequences (no longer needed) to destination
    if [[ ! -d ${DESTDIR}/initqc ]]; then
        mkdir ${DESTDIR}/initqc
    fi
    skip="N"
    chklog "completed_concatenation_mv_to_dest"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" && ${needsconcat} == "T" ]]; then
        mv -t ${DESTDIR}/initqc ${output}/concatenate
        echo "completed_concatenation_mv_to_dest"
    fi
    
    
    if [[ ! -d ${DESTDIR}/initqc/trimmomatic ]]; then
        mkdir ${DESTDIR}/initqc/trimmomatic
    fi
    
    ## Remove low-complexity sequences using fqtrim
    skip="N"
    chklog "fqtrim_1p_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
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
            chklog "${readpair}_complete"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
            
                ${FQTRIM} -l 30 -p ${NUMCPUS} -D -o fqtrimmed.fq --outdir ${output}/fqtrim1 -r ${LOGLOC}/fqtrimlog1p.txt ${readpair}
                
                # file management: move trimmomatic paired files (no longer needed) to destination
                mv -t ${DESTDIR}/initqc/trimmomatic ${file}
                mv -t ${DESTDIR}/initqc/trimmomatic ${file/_1P/_2P}
                
                echo "${readpair}_complete"
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> fqtrim_1p_complete"
    fi
    
    
    ## incorporate unpaired sequences, if option is Y
    skip="N"
    chklog "fqtrim_1u_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        if [[ ${USEUNPAIRED} == "Y" ]]; then
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Running fqtrim on unpaired reads..."
        else
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Moving unpaired reads to destination..."
        fi

            # trimmomatic outputs unpaired sequences as ./*_1U.fq and ./*_2U.fq
        for file in ${output}/trimmomatic/*U.fq; do
            skip="N"
            chklog "${output}/fqtrim1/${file##*/}_complete"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then

                if [[ ${USEUNPAIRED} == "Y" ]]; then
                    ${FQTRIM} -l 30 -p ${NUMCPUS} -D -o fqtrimmed.fq --outdir ${output}/fqtrim1 -r ${LOGLOC}/fqtrimlog1u.txt ${file}
                fi
                
                # file management: move trimmomatic unpaired files (no longer needed) to destination
                mv -t ${DESTDIR}/initqc/trimmomatic ${file}

                echo "${output}/fqtrim1/${file##*/}_complete"
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> fqtrim_1u_complete"

    fi
    
    ## file management: all trimmomatic files have been processed, remove folder
    skip="N"
    chklog "trimmomatic_folder_moved"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        for file in ${output}/trimmomatic/*; do
            if [[ -f ${file} ]]; then
                mv -t ${DESTDIR}/initqc/trimmomatic ${file}
            fi
        done
        rmdir ${output}/trimmomatic
        echo "trimmomatic_folder_moved"
    fi
    
    
    ## run FastQC on trimmed reads
    skip="N"
    chklog "trim1_fastqc_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
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
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
            
                ${FASTQC} -t ${NUMCPUS} -o ${output}/fastqc/trim1 ${seqfile}
                echo ${output}/fastqc/trim1/${seqfile##*/}
                
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> trim1_fastqc_complete"
    fi
    
    
    if [[ ! -d ${DESTDIR}/initqc/fqtrim1 ]]; then
        mkdir ${DESTDIR}/initqc/fqtrim1
    fi
    
    ## rerun fqtrim on paired sequences, but ignore pairs
        # to preserve pair order, fqtrim will leave single-nucleotide sequences in the files
        # rerunning fqtrim is required to remove these junk sequences
    skip="N"
    chklog "fqtrim_2p_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
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
            chklog "${output}/fqtrim2/${seqfile##*/}_trimmed2"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then

                ${FQTRIM} -l 30 -p ${NUMCPUS} -D -o 2.fq --outdir ${output}/fqtrim2 -r ${LOGLOC}/fqtrimlog2p.txt ${seqfile}
                
                # file management: paired fqtrim1 reads are no longer needed -- move to destination
                mv -t ${DESTDIR}/initqc/fqtrim1 ${seqfile}
                
                echo "${output}/fqtrim2/${seqfile##*/}_trimmed2"
                
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> fqtrim_2p_complete"
    fi
    
    
    ## run FastQC on retrimmed reads
    skip="N"
    chklog "trim2_fastqc_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
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
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
            
                ${FASTQC} -t ${NUMCPUS} -o ${output}/fastqc/trim2 ${seqfile}
                echo ${output}/fastqc/trim2/${seqfile##*/}
                
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> trim2_fastqc_complete"
    fi
    
    
    if [[ ! -d ${DESTDIR}/initqc/fqtrim2 ]]; then
        mkdir ${DESTDIR}/initqc/fqtrim2
    fi
    
    ## reorder paired sequences
        # paired sequences need to be in the same order in each pair file for the next step
    skip="N"
    chklog "interleave_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Performing sequence pair matching..."

        # create folder for interleave step
        if [[ ! -d ${output}/interleave ]]; then
            mkdir ${output}/interleave
        fi

        echo "matching pairs..."
        # perform interleave
        for file in ${output}/fqtrim2/*_1P.fqtrimmed.2.fq; do
            reverse=$(echo ${file/_1P/_2P} | tr -d '\n')
            basename=${file##*/}
            base=${output}/interleave/${basename%_*}
            
            skip="N"
            chklog "${output}/interleave/${base}_out_unpaired.fastq_interleaved"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then

                python3 ${INTERLEAVE} ${file} ${reverse} ${base}
                # output files are ./*_out_pairs_fwd.fastq ./*_out_pairs_rev.fastq and ./*_out_unpaired.fastq
                
                # file management: fqtrim2 files are no longer needed -- move them to the destination
                mv -t ${DESTDIR}/initqc/fqtrim2 ${file}
                mv -t ${DESTDIR}/initqc/fqtrim2 ${reverse}
                
                echo "${output}/interleave/${base}_out_unpaired.fastq_interleaved"
                
            fi
        done
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> interleave_complete"
    fi
    
    
    ## end-of-block file management
        # trimmomatic folder has already been moved
        # unpaired reads remain in fqtrim1 folder
        # no reads remain in fqtrim2, but the folder still exists
        # fwd, reverse, and unpaired reads remain in fqtrim2
        # quality data remains in fastqc
    skip="N"
    chklog "initqc_file_management_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> File Management..."

        # remove fqtrim2 folder
        for file in ${output}/fqtrim2/*; do
            if [[ -f ${file} ]]; then
                mv -t ${DESTDIR}/initqc/fqtrim2 ${file}
            fi
        done
        if [[ -d ${output}/fqtrim2 ]]; then
            rmdir ${output}/fqtrim2
        fi
        
        # copy all remaining output files to destination
        cp -r -t ${DESTDIR}/initqc ${output}/*
        
        # move needed output files to input
            # paired files
        for file in ${output}/interleave/*_out_pairs_*.fastq; do
            if [[ -f ${file} ]]; then
                mv -t ${input} ${file}
            fi
        done
        if [[ ${USEUNPAIRED} == "Y" ]]; then
            # unpaired fqtrim1 files
            for file in ${output}/fqtrim1/*U.fqtrimmed.fq; do
                if [[ -f ${file} ]]; then
                    mv -t ${input} ${file}
                fi
            done
            # unpaired interleaved files
            for file in ${output}/interleave/*_out_unpaired.fastq; do
                if [[ -f ${file} ]]; then
                    mv -t ${input} ${file}
                fi
            done
        fi
        
        # empty ouptut folder
        emptydir ${output}

        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> initqc_file_management_complete"
    fi
    
    
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> initqc_complete"
}



### initial sample processing
initproc() {

    temploc=${output}/temp
    if [[ ! -d ${temploc} ]]; then
        mkdir ${temploc}
    fi
    indexloc=${output}/genome_index
    if [[ ! -d ${indexloc} ]]; then
        mkdir ${indexloc}
    fi
    alignloc=${output}/aligned
    if [[ ! -d ${alignloc} ]]; then
        mkdir ${alignloc}
    fi
    unmappeddir=${output}/unmapped
    if [[ ! -d ${unmappeddir} ]]; then
        mkdir ${unmappeddir}
    fi
    
    # create genome index
    skip="N"
    chklog "genome_index_completed"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Indexing the genome..."
        
        ${HISAT2}-build -q -p ${NUMCPUS} --seed 12345 ${GENOME} ${indexloc}/index
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> genome_index_completed"
    fi
    
    
    # the final qality controlled sequences should all be in the local input folder
    # samples should have the same prefix and may have the suffixes: _1U.fqtrimmed.fq _2U.fqtrimmed.fq _out_pairs_fwd.fastq _out_pairs_rev.fastq _out_unpaired.fastq
    reads1=(${input}/*_out_pairs_fwd.fastq)
    reads2=("${reads1[@]/_fwd./_rev.}")
    if [[ ${USEUNPAIRED} == "Y" ]]; then
        unpaired1=("${reads1[@]/_pairs_fwd/_unpaired}")
        unpaired2=(${input}/*_1U.fqtrimmed.fq)
        unpaired3=("${unpaired2[@]/_1U/_2U}")
    fi
    
    for ((i=0; i<=${#reads1[@]}-1; i++ )); do
        sample="${reads1[$i]##*/}"
        sample="${sample%%_out_*}"
        echo "Sample name: " ${sample}

        # in case of interruption, skip finished files
        skip="N"
        chklog "initproc_${sample}_complete"
        if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
            skip="Y"
        fi
        if [[ ${skip} == "N" ]]; then
        
            unmappedloc=${unmappeddir}/${sample}_unmapped.txt

            unpairedBoth=${unpaired1[$i]},${unpaired2[$i]},${unpaired3[$i]}
            stime=`date +"%Y-%m-%d %H:%M:%S"`
            echo "[$stime] Processing sample: $sample"
            echo "[$stime]    * Alignment of reads to genome (HISAT2)"
            
            if [[ ${USEUNPAIRED} == "Y" ]]; then
                $HISAT2 -p $NUMCPUS --dta -x ${indexloc}/index \
                 -1 ${reads1[$i]} \
                 -2 ${reads2[$i]} \
                 -S ${temploc}/${sample}.sam 2>${alignloc}/${sample}.alnstats \
                 -U ${unpairedBoth} \
                 --seed 12345 --un ${unmappedloc} \
                 --un-conc ${unmappedloc}
                 #defaults kept for strandedness, alignment, and scoring options
                 
             else
                $HISAT2 -p $NUMCPUS --dta -x ${indexloc}/index \
                 -1 ${reads1[$i]} \
                 -2 ${reads2[$i]} \
                 -S ${temploc}/${sample}.sam 2>${alignloc}/${sample}.alnstats \
                 --seed 12345 --un ${unmappedloc} \
                 --un-conc ${unmappedloc}
                 #defaults kept for strandedness, alignment, and scoring options
             fi

            echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Alignments conversion (SAMTools)"
            if [[ "$newsamtools" ]]; then
             $SAMTOOLS view -S -b ${temploc}/${sample}.sam | \
              $SAMTOOLS sort -@ $NUMCPUS -o ${alignloc}/${sample}.bam -
            else
             $SAMTOOLS view -S -b ${temploc}/${sample}.sam | \
              $SAMTOOLS sort -@ $NUMCPUS - ${alignloc}/${sample}
            fi
            $SAMTOOLS index ${alignloc}/${sample}.bam

            echo "..removing intermediate files"
            rm ${temploc}/${sample}.sam

            echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
            $STRINGTIE -p $NUMCPUS -o ${alignloc}/${sample}.gtf \
             -l ${sample} ${alignloc}/${sample}.bam \
             -c 1.5
             # -m kept at default of 0.01
             
             echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> initproc_${sample}_complete"
        else
           echo "    ${sample} logged as completed, skipping..."
        fi
    done
    
    
    ## end-of-block file management
    skip="N"
    chklog "initproc_file_management_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> File Management..."
        
        # empty local input folder (a copy should have been saved elsewhere earlier)
        emptydir ${input}
        
        # empty and delete temp folder
        emptydir ${temploc}
        rmdir ${temploc}
        
        # copy final files to local input folder
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Copying necessary files to input folder"
        for file in ${alignloc}/*.gtf; do
            if [[ -f ${file} ]]; then
                basename=${file##*/}
                cp ${file} ${input}/${basename}
            fi
        done
        for file in ${alignloc}/*.bam; do
            if [[ -f ${file} ]]; then
                basename=${file##*/}
                cp ${file} ${input}/${basename}
            fi
        done

        # move contents of output folder to destination folder
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Moving files to destination folder"
        if [[ ! -d ${DESTDIR}/initproc ]]; then
            mkdir ${DESTDIR}/initproc
        fi

        mv -t ${DESTDIR}/initproc ${output}/*
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> initproc_file_management_complete"
    fi
    
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> initproc_complete"
}



### estimate abundances and remove rRNA -- used in blocks 3 and 4
estabund() {
    thisname=$1
    thisout=$2

    # create merge list
    skip="N"
    chklog "${thisname}_transcripts_merged"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Merge transcripts (StringTie)"
        
        # create file to hold list of files to merge
        touch ${thisout}/mergelist.txt
        ls -1 ${input}/*.gtf > ${thisout}/mergelist.txt

        $STRINGTIE --merge -p $NUMCPUS \
            -o ${thisout}/stringtie_merged.gtf ${thisout}/mergelist.txt
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> ${thisname}_transcripts_merged"
    fi

    ## estimate transcript abundance
    skip="N"
    chklog "${thisname}_abundance_estimation_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Estimate abundance for each sample (StringTie)"
        
        samples=${input}/*.bam
        for ((i=0; i<=${#samples[@]}-1; i++ )); do
            sample="${samples[$i]%%.*}"
            sample="${sample##*/}"
            
            skip="N"
            chklog "${thisname}_${sample}_abundance_estimated"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
            
                if [ ! -d ${thisout}/abund ]; then
                    mkdir -p ${thisout}/abund
                fi
                if [ ! -d ${thisout}/abund/${sample} ]; then
                    mkdir -p ${thisout}/abund/${sample}
                fi
                
                if [[ ${thisname} == "listprep" ]]; then
                    set +e # do not exit on errors from the next command
                fi
                
                $STRINGTIE -e -B -p $NUMCPUS -G ${thisout}/stringtie_merged.gtf \
                -o ${thisout}/abund/${sample}/${sample}.gtf ${input}/${sample}.bam
                
                if [[ ${thisname} == "listprep" ]]; then
                    # handle cases where StringTie returns an error because the file is too small
                    status=$?

                    set -e # return to default of exiting on any errors

                    if (( status > 0 )); then
                        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> StringTie error: sample ${sample} will be removed from analysis"

                        REMOVEALWAYS="${REMOVEALWAYS} ${sample}"

                    fi
                fi

                echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${thisname}_${sample}_abundance_estimated"
            fi
        done
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${thisname}_abundance_estimation_complete"
    fi

    # create transcript FASTA file
    skip="N"
    chklog "${thisname}_transcriptome_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Generate fasta transcript file from .gtf"
        
        $GFFREAD -w ${thisout}/transcriptome.fa -g ${GENOME} ${thisout}/stringtie_merged.gtf
        
        echo  [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${thisname}_transcriptome_complete"
    fi

    ### BLAST against ribosomal database

    # make rRNA database
    
    skip="N"
    chklog "rRNA_DB_created"
    if [[ ${chkresult} == "T" && -d ${databases}/rrnadb ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Create rRNA database (makeblastdb)"
        
        if [[ ! -d ${databases}/rrnadb ]]; then
            mkdir ${databases}/rrnadb
        fi
    
        $BLASTDIR/makeblastdb -in ${RRNAFILE} -out ${databases}/rrnadb -dbtype nucl
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> rRNA_DB_created"
    fi

    # BLAST rRNA database
    skip="N"
    chklog "${thisname}_rrna_blast_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Query rRNA database with transcroptome (BLASTN)"

        $BLASTNAPP -query ${thisout}/transcriptome.fa \
            -db ${databases}/rrnadb \
            -out ${thisout}/rrna_blastn_results.csv \
            -evalue 1e-6 -num_threads ${NUMCPUS} \
            -num_alignments 1 -outfmt "6 qseqid stitle sacc evalue pident bitscore length qstart qend sstart send"
            
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${thisname}_rrna_blast_complete"
    fi

    ## remove matched transcripts from ctab files
    skip="N"
    chklog "${thisname}_rRNA_removed"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Removing rRNA..."

        # save a copy of the original tables in a new folder
        skip="N"
        chklog "${thisname}_ctabs_copied"
        if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
            skip="Y"
        fi
        if [[ ${skip} == "N" ]]; then
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Save a copy of original expression tables"

            if [[ -d ${thisout}/rrna_free ]]; then
                # we want the contents of /abund to be copied to /rrna_free
                # if /rrna_free exists, cp will instead copy /abund to /rrna_free/abund
                # to simplify the process, remove /rrna_free, if it exists
                emptydir ${thisout}/rrna_free
                rm ${thisout}/rrna_free
            fi
            
            cp -R ${thisout}/abund ${thisout}/rrna_free
            
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${thisname}_ctabs_copied"
        fi

        # use R script to remove rRNA and overwrite the .ctab files
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Remove rRNA from expression tables (R)"

        Rscript ${REMOVERRNA} ${thisout}/blastn_results.csv ${thisout}/rrna_free ${LOGDIR} ${thisname}
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${thisname}_rRNA_removed"
    fi
    
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${thisname}_estabund_complete"
}



### automatically generate removelists based on user-specified cutoffs
listprep() {
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Automatically generating remove lists..."
    
    # create temporary folder for these processes
    temp=${output}/temp
    if [[ ! -d ${temp} ]]; then
        mkdir ${temp}
    fi
    
    ## estimate abundance and remove rRNA
    skip="N"
    chklog "listprep_estabund_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        estabund "listprep" ${temp}
    fi
    
    ## for each auto remove list, add samples to that list that are below the cutoff and add any REMOVEALWAYS samples
    skip="N"
    chklog "removelists_updated"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Checking removelists..."
        
        # run through transcript counter once to generate count file
        Rscript ${TRANSCRIPTCOUNTER} ${input} ${SETNAMES[i]} ${LOGLOC} ${temp} 0

        removecounter=0
        newremovelists=()
        for ((i=0; i<${#SETNAMES[@]}; i++)); do
            # check if setname begins with "auto-"
            if [[ ${SETNAMES[i]} == "auto-all" ]]; then
                currentremovelist="${REMOVEALWAYS}"
                
            # if auto...
            elif [[ ${SETNAMES[i]} == auto-* ]]; then
                # extract cutoff from setname
                cutoff=${SETNAMES[i]}
                cutoff=${cutoff/"auto-"/}
                
                Rscript ${TRANSCRIPTCOUNTER} ${temp}/transcript_counts.csv ${SETNAMES[i]} ${LOGLOC} ${temp} ${cutoff}
                
                # extract list from removelist.txt and store as currentremovelist
                currentremovelist=(cat ${temp}/removelist.txt) ####################################################### does this work?
                
            else
            # if not auto...
                # extract next removelist and store as currentremovelist
                currentremovelist="${REMOVELISTS[removecounter]}"
                removecounter=${removecounter}+1
                
            fi
            # add currentremovelist to newremovelists
            echo ${currentremovelist}
            newremovelists+=("${currentremovelist}")
        done
        
        REMOVELISTS=${newremovelists}
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> removelists_updated"
    fi
    
    # end-of-block file management: move transcript_counts.csv to destination and remove /temp
    if [[ -f ${temp}/transcript_counts.csv ]]; then
        mv ${temp}/transcript_counts.csv ${DESTDIR}/transcript_counts.csv
    fi
    if [[ -d ${temp} ]]; then
        emptydir ${temp}
        rm ${temp}
    fi
    
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> listprep_complete"
}


### analysis and visualization
mojo() {

    setname=$1
    removelist=$2
    bg_cov=$3
    bg_adjvars=$4
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Beginning analysis of sample set ${setname}"
    
    analysis=${DESTDIR}/analysis
    destination=${DESTDIR}/analysis/${setname}
    mojoinput=${input}/${setname}
    direcs=("${analysis}" "${destination}" "${mojoinput}")
    
    # file management strategy
        # merge files, est abundance, create transcriptome, remove rRNA
            # do in output/abundances
            # store abundances in output/abundances/abund
            # store removed rRNA abundances in output/abundances/rrna_free
            # move tome to destination/products
            # move count file to destination/calculations
            # move abund to destination/calculations/raw_abundances
            # move rrna_free to destdir/calculations/rrna_removed
    abundances=${output}/abundances
    productdest=${destination}/products
    calcdest=${destination}/calculations
    abunddest=${destination}/calculations/raw_abundances
    rrnadest=${destination}/calculations/rrna_removed
    direcs+=("${abundances}" "${productdest}" "${calcdest}" "${abunddest}" "${rrnadest}")
        # busco test
            # must be done in output/abundances
            # copy output/abundances/busco_output_enoki*/short_summary*.txt to destdir/products/short_summary*.txt
            # move output/abundances/busco_* to destination/calculations/busco/busco_*
    buscodest=${destination}/calculations/busco
    direcs+=("${buscodest}")
    # copy rrna_removed directories to mojoinput/abund
    # copy tome to mojoinput
    # move files from output as described above
    # move busco files from output as described above
    mojoinabund=${mojoinput}/abund
    direcs+=("${mojoinabund}")
        # ballgown
            # do in output/ballgown
            # move PCA plots to destdir/products/plots
            # move remainder to destination/calculations/ballgown
    plotsdest=${destination}/products/plots
    ballgownout=${output}/ballgown
    ballgowndest=${destination}/calculations/ballgown
    direcs+=("${prodplot}" "${ballgownout}" "${ballgowndest}")
        # proteome
            # do in output/proteome
            # move lists to destdir/products/proteome
            # move rest to destdir/calculations/proteome
    protout=${output}/proteome
    protproducts=${destination}/products/proteome
    protdest=${destination}/calculations/proteome
    direcs+=("${protout}" "${protproducts}" "${protdest}")
        # heat maps/barplots
            # do in output/plots
            # move heat maps to destdir/products/plots
    plotsout=${output}/plots
    direcs+=("${plotsout}")
        # panther
            # do in output/panther
            # move lists to destdir/products/go
            # move rest to destdir/calculations/panther
    pantherout=${output}/panther
    pantherdest=${destination}/calculations/panther
    godest=${destination}/products/go
    direcs+=("${pantherout}" "${pantherdest}" "${godest}")
    # move files from output as described above
    # empty and delete mojoinput
    
    # file location for anything I forgot to add to this strategy
    otherdest=${destination}/other
    direcs+=("${otherdest}")
    
    for direc in ${direcs[@]}; do
        if [[ ! -d ${direc} ]]; then
            mkdir ${direc}
        fi
    done
    
    
    ### estimate abundances and remove rRNA
    skip="N"
    chklog "${setname}_estabund_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
    
        estabund ${setname} ${abundances}
        
    fi
    

    ### busco completeness test
    skip="N"
    chklog "${setname}_busco_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Testing Busco completeness..."

        busconame="busco_output_enoki"

        if [[ -d ${BALLGOWNLOC}/${busconame} ]]; then
            busconame="busco_output_enoki_"`date +"%Y-%m-%d_%H-%M-%S"`
        fi

        docker run -u $(id -u) -v ${abundances}:/busco_wd/ ezlabgva/busco:v4.1.4_cv1 busco \
            --mode transcriptome \
            --in /busco_wd/transcriptome.fa \
            --out ${busconame} \
            --lineage_dataset ${BUSCODATASET}
        # note: -v [specify complete file path where you want busco to live]:/busco_wd/
        # note: the working directory (specified in -v) must be /busco_wd/
        # note: all file paths (except in -v) must be relative to the image root directory. In this case, -v maps the image root (/busco_wd/) to ${abundances}/ . As a result, the transcriptome file must be referenced as /busco_wd/transcriptome.fa because it is found in ${abundances}/transcriptome.fa
        # note: --out file must have a name that has never been used for a busco run. Ex: busco was previously run on the same machine with "--out busco_output" so this run must use a new name, "--out busco_output_enoki"
        # note: if script must be run more than once, delete the output directories before running OR make sure to run the busconame if statement above
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_busco_complete"
    fi
    
    
    ## file management
    skip="N"
    chklog "${setname}_mojo_files_1_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Performing file management..."

        # transcriptome
        if [[ -f ${abundances}/transcriptome.fa ]]; then
            cp ${abundances}/transcriptome.fa ${mojoinput}/transcriptome.fa
            mv ${abundances}/transcriptome.fa ${destination}/products/transcriptome.fa
        fi
        # BUSCO results
        buscoresult=${output}/abundances/${busconame}/short_summary*.txt
        buscoresultname=${buscoresult##*/}
        if [[ -f ${buscoresult} ]]; then
            cp ${buscoresult} ${productdest}/BUSCO_${buscoresultname}
        fi
        # BUSCO calculations
        for dir in ${output}/abundances/busco_*; do
            mv -t ${buscodest} ${dir} 
        done
        # raw abundances
        if [[ -d ${abundances}/abund ]]; then
            mv ${abundances}/abund ${calcdest}/raw_abundances
        fi
        # rRNA removed abundances
        if [[ -d ${abundances}/rrna_free ]]; then
            mv ${abundances}/rrna_free ${calcdest}/rrna_removed
        fi
        # anything else that might be left in output/abundances
        for file in ${abundances}/*; do
            mv -t ${otherdest} file
        done
        # delete output/abundances
        if [[ -d ${abundances} ]]; then
            emptydir ${abundances}
            rm ${abundances}
        fi
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> ${setname}_mojo_files_1_complete"
    fi


    ### use ballgown to create differential expression tables
    skip="N"
    chklog "${setname}_ballgown_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Generate the DE tables (Ballgown)"
        
        if [[ ! -d ${ballgownout}/PCA ]]; then
            mkdir ${ballgownout}/PCA
        fi
        if [[ ! -d ${ballgownout}/lists ]]; then
            mkdir ${ballgownout}/lists
        fi

        Rscript ${BALLGOWN} ${PHENODATA} ${mojoinabund} ${setname} ${bg_cov} ${bg_adjvars} ${PCAPAIRS} ${ballgownout} ${LOGLOC}
    
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_ballgown_complete"
    fi


    ## BLAST DEGs
    skip="N"
    chklog "${setname}_gene_naming_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then

        # create BLAST database
        skip="N"
        chklog "protein_db_created"
        if [[ ${chkresult} == "T" ]]; then
            skip="Y"
        fi
        if [[ ${skip} == "N" ]]; then
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Create protein database (makeblastdb)"
            if [[ ! -d ${databases}/protein ]]; then
                mkdir ${dataabses}/protein
            fi
            
            $BLASTDIR/makeblastdb -in ${UNIPROTFILE} -out ${databases}/protein -dbtype prot
            
            echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> protein_db_created"
        fi

        for deglist in ${ballgownout}/lists/*.txt; do
            currentDEGList=${deglist##*/}
            currentDEGs=${currentDEGList%.txt}
            
            skip="N"
            chklog "${setname}_${currentDEGs}_gene_naming_complete"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                echo "Current DEG list: ${currentDEGs}"

                skip="N"
                chklog "${setname}_${currentDEGs}_matchgenes_complete"
                if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                    skip="Y"
                fi
                if [[ ${skip} == "N" ]]; then
                    # biopython: match MSTRNG numbers from DEGs and .fa sequences
                    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Match MSTRNGs to genes for DEGs (Biopython)"
                    filename=${ballgownout}/lists/${currentDEGs}.fa
                    if [[ -f ${filename} ]]; then
                        rm ${filename}
                    fi

                    # 3 arguments:
                    # output file name (fasta format)
                    # fasta transcriptome with MSTRNG/sequence pairs
                    # .txt list of DEGs separated by line
                    python3 ${MATCHGENES} ${filename} ${mojoinput}/transcriptome.fa ${deglist}
                    
                    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_${currentDEGs}_matchgenes_complete"
                fi

                # search BLAST database
                skip="N"
                chklog "${setname}_${currentDEGs}_blastx_complete"
                if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                    skip="Y"
                fi
                if [[ ${skip} == "N" ]]; then
                    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Query BLAST database for gene names (BLASTX)"

                    $BLASTXAPP -query ${filename} \
                        -db ${databases}/protein \
                        -out ${protout}/${setname}_${currentDEGs}_blastx_results.csv \
                        -evalue 1e-3 -num_threads ${NUMCPUS} \
                        -num_alignments 1 -outfmt "6 qseqid stitle sacc evalue pident bitscore length qstart qend sstart send"
                    
                    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_${currentDEGs}_blastx_complete"
                fi

                # extract infromation from UniProt stitle
                Rscript ${FIXUNIPROT} ${protout}/${setname}_${currentDEGs}_blastx_results.csv ${protout}/${setname}_${currentDEGs}_blastx_results_readable.csv ${LOGLOC}
                
                if [[ ! -d ${protout}/named ]]; then
                    mkdir ${protout}/named
                fi

                # match protein names to output files
                for bgresults in ${ballgownout}/${currentDEGs}*csv; do
                    resultname=${bgresults##*/}
                    Rscript ${NAMEGENES} ${protout}/${setname}_${currentDEGs}_blastx_results_readable.csv ${bgresults} ${protout}/named/${resultname} ${LOGLOC}
                done
                
                echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_${currentDEGs}_gene_naming_complete"
            fi
        done
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_gene_naming_complete"
    fi


######################################################################
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


     ## generate heat maps or bar plots
    skip="N"
    chklog "${setname}_plots_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Creating Heat Maps and/or Bar Plots"

        for resultfile in ${protout}/named/*_expr.csv; do
            basename=${resultfile##*/}
            base=${basename%.csv}

            
            skip="N"
            chklog "${setname}_all_heatmap_complete"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                # map all genes
                Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${plotsout}/${base}_heatmap_all.pdf "all" "all" "all" "F" ${LOGLOC}
                
                echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_all_heatmap_complete"
            fi

            # some sets may not include the features of interest and could therefore return an error, but keep running anyway
            set +e

            skip="N"
            chklog "${setname}_gene_plots_complete"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                # map specific genes
                for ((k=0; k<=${#MAPGENELISTS[@]}-1; k++)); do
                       Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${plotsout}/${base}_heatmap_genes_${k}.pdf ${MAPGENELISTS[k]} "all" "all" "F" ${LOGLOC}
                done
                echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_gene_plots_complete"
            fi

            skip="N"
            chklog "${setname}_organism_plots_complete"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                # map genes from specific organisms
                for ((k=0; k<=${#MAPORGLISTS[@]}-1; k++)); do
                       Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${plotsout}/${base}_heatmap_orgs_${k}.pdf "all" "all" ${MAPORGLISTS[k]} "F" ${LOGLOC}
                done
                echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_organism_plots_complete"
            fi

            skip="N"
            chklog "${setname}_sample_plots_complete"
            if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
                skip="Y"
            fi
            if [[ ${skip} == "N" ]]; then
                # map specific samples
                for ((k=0; k<=${#MAPSAMPLELISTS[@]}-1; k++)); do
                       Rscript ${HEATMAP} ${resultfile} ${PHENODATA} ${plotsout}/${base}_heatmap_samples_${k}.pdf "all" ${MAPSAMPLELISTS[k]} "all" "F" ${LOGLOC}
                done
                echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_sample_plots_complete"
            fi

            # reset default behavior of exiting on errors
            set -e

        done
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_plots_complete"
    fi
    
    
    ## file management
    skip="N"
    chklog "${setname}_mojo_files_2_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Performing file management..."
        
        # copy needed files to new input folder
        mojogoin=${input}/mojogo
        if [[ ! -d ${mojogoin} ]]; then
            mkdir ${mojogoin}
        fi
        mv ${mojoinput}/transcriptome.fa ${mojogoin}/transcriptome.fa
        cp -t ${mojogoin} ${protout}/named/*.csv
        
        # move output files to destinations
        mv -t ${plotsdest} ${ballgownout}/PCA/*
        rm ${ballgownout}/PCA
        mv -t ${ballgowndest} ${ballgownout}/*
        mv -t ${protproducts} ${protout}/*_readable.csv
        mv -t ${protdest} ${protout}/*
        mv -t ${plotsdest} ${plotsout}/*
        
        # in case anything was forgotten, move it to the "other" folder
        mv -t ${otherdest} ${output}/*.*
        
        # remove mojoinput folder
        emptydir ${mojoinput}
        rm ${mojoinput}
        
        # remove output folders
        emptydir ${output}
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> ${setname}_mojo_files_2_complete"
    fi
    
    
    ## GO analysis
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Generating tables for GO analysis..."

    # for each DEG list and reference list...
    
    # PantherScore won't work unless the working directory is the directory that contains PantherScore and the program's lib folder
    cd ${PANTHERSCORE%/*}

    #######################################################
    # shouldn't need to do because the full transcriptome is included in deglists
    #echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Converting UniProtKB IDs to FASTA..."
    ## for each named protein, get the NCBI Protein database sequence and add it to a FASTA file
    #python3 ${UNIPROTFASTA} ${pantherout}/tome_ncbi_list.csv ${pantherout}/tome_fasta.fa ${NCBIAPI} ${mojogoin}/${setname}_tome_blastx_results_readable.csv
    #
    #echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Running PantherScore..."
    ## score FASTA file against PANTHER HMMs to map to PANTHER IDs
    #perl ${PANTHERSCORE} -l ${PANTHERLIBDIR} -D B -i ${pantherout}/tome_fasta.fa -o ${pantherout}/tome_panther_mapping.txt -n -s

    for deglist in ${mojogoin}/*_expr.csv; do
        base=${deglist##*/}
        basename=${base%.csv}
        echo "Current list: ${basename}"

        # for each named protein, get the NCBI Protein database sequence and add it to a FASTA file
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Converting UniProtKB IDs to FASTA..."
        python3 ${UNIPROTFASTA} ${pantherout}/${basename}_ncbi_list.csv ${pantherout}/${basename}_fasta.fa "" ${deglist}
        
        # score FASTA file against PANTHER HMMs to map to PANTHER IDs
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Running PantherScore..."
        perl ${PANTHERSCORE} -l ${PANTHERLIBDIR} -D B -i ${pantherout}/${basename}_fasta.fa -o ${pantherout}/${basename}_panther_mapping.csv -n -s
            # output is tab-delimited: sequence ID, panther acc, panther family/subfamily, HMM e-value, HMM bitscore, alignment range
        
        # map PANTHER IDs to fpkm tables
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Matching PANTHER IDs to FPKM tables..."
        if [[ ! -d ${pantherout}/final/${basename} ]]; then
            mkdir ${pantherout}/final/${basename}
        fi
        Rscript ${PANTHERFPKM} ${deglist} ${pantherout}/${basename}_panther_mapping.csv ${pantherout}/${basename}_ncbi_list.csv ${PHENODATA} "ROEHL,myc,sti,pil,gil,you,pri,cul,nor" ${basename}_panther_mapping ${pantherout}/final/${basename}
    done
    
    cd ${WRKDIR}
    
    
    # end-of-block file management
    skip="N"
    chklog "${setname}_mojo_files_3_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "##> Performing file management..."
        
        # move PANTHER files
        mv -t ${godest} ${pantherout}/final/*
        rm ${pantherout}/final
        mv -t ${pantherdest} ${pantherout}
        rm ${pantherout}
        
        # remove mojogoin
        emptydir ${mojogoin}
        rm ${mojogoin}
        
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> ${setname}_mojo_files_3_complete"
    fi
    
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> mojo_${setname}_complete"
}



#### PIPELINE PROCESS



# initial setup

## load variables
SCRIPTARGS="$@"
cd ${WRKDIR}

if [[ ! -f ${SCRIPTARGS[0]} ]]; then
    usage
    echo "Error: configuration file (rnaseq_pipeline.config.sh) missing or not included as argument in command line!"
    exit 1
else
    source ${SCRIPTARGS[0]}
fi

set -e

export PATH="$PATH:${BLASTDIR}:${TRIMMOMATICADAPTERS}"




# check for previous log file - if found, ask user whether to start anew and overwrite or to skip completed files
LOGFILE=${LOGLOC}/rnaseq_pipeline.log
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



# create local input, output, and database folders
# but first, check to see whether those folders exist -- if they do, prompt user whether files should be overwritten
input=${WRKDIR}/input
output=${WRKDIR}/output
databases=${WRKDIR}/databases
checkfolders=("${input}" "${output}" "${databases}")
for chkfldr in ${checkfolders[@]}; do
    if [[ -d ${chkfldr} ]]; then
        if [[ ${resume} == "N" ]]; then
            fldrname=${chkfldr##*/}
            echo "A folder named ${fldrname} already exists in the working location. If you continue, its files will be deleted. Continue anyway? (Y/N)"
            read continueanyway
            if [[ ${continueanyway} == "Y" ]] || [[ ${continueanyway} == "y" ]] || [[ ${continueanyway} == "Yes" ]] || [[ ${continueanyway} == "YES" ]] || [[ ${continueanyway} == "yes" ]]; then
                echo "Deleting contents and continuing..."
                emptydir ${chkfldr}
            else
                echo "Please choose a new working directory or save the files you need in another location. When ready, restart this pipeline."
                exit 1
            fi
        fi
    else
        mkdir ${chkfldr}
    fi
done



echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $SCRIPTARGS

# (if skip) check if initial QC was previously completed - if so, skip initial QC
skip="N"
chklog "initqc_complete"
if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
    skip="Y"
fi
if [[ ${skip} == "N" ]]; then

    # check that required files are present for block 1
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Checking programs for part 1..."
    chkprog 1

    # do initial QC
    initqc 2>&1 | tee -a $LOGFILE
    
    chklog "initqc_complete"
    if [[ ${chkresult} == "F" ]]; then
        exit 1
    fi
    
fi



# (if skip) check if initial analysis was previously completed - if so, skip inital analysis
skip="N"
chklog "initproc_complete"
if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
    skip="Y"
fi
if [[ ${skip} == "N" ]]; then

    # check that required files are present for block 2
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Checking programs for part 2..."
    chkprog 2

    # do initial analysis
    initproc 2>&1 | tee -a $LOGFILE
    
    chklog "initproc_complete"
    if [[ ${chkresult} == "F" ]]; then
        exit 1
    fi
    
fi



# (if skip) check if initial analysis was previously completed - if so, skip inital analysis
skip="N"
chklog "listprep_complete"
if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
    skip="Y"
fi
if [[ ${skip} == "N" ]]; then
    # check that required files are present for block 3
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Checking programs for part 3..."
    chkprog 3

    # prepare removelists
    listprep 2>&1 | tee -a $LOGFILE
    
    chklog "listprep_complete"
    if [[ ${chkresult} == "F" ]]; then
        exit 1
    fi
    
fi



# check that required files are present for block 4
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Checking programs for part 4..."
chkprog 4
# for each data subset...
    # (if skip) check if initial analysis was previously completed - if so, skip inital analysis
    # run analysis pipeline (if skip, skip completed files)
## run rest of pipeline on various sample sets
for ((index=0; index<${#SETNAMES[@]}; index++ )); do
    skip="N"
    chklog "mojo_${SETNAMES[index]}_complete"
    if [[ ${resume} == "Y" && ${chkresult} == "T" ]]; then
        skip="Y"
    fi
    if [[ ${skip} == "N" ]]; then
        
        # perform main analyses
        mojo ${SETNAMES[index]} ${REMOVELISTS[index]} ${COVARIATES[index]} ${ADJVARSETS[index]} 2>&1 | tee -a $LOGFILE
        
        chklog "mojo_${SETNAMES[index]}_complete"
        if [[ ${chkresult} == "F" ]]; then
            exit 1
        fi
        
    fi
done

# end-of-pipeline file management

# move forgotten files to destination/other
mv -t ${otherdest} ${output}/*

# remove local input, output, and database folders
rm ${output}
emptydir ${input}
rm ${input}
emptydir ${databases}
rm ${databases}

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> PIPELINE-COMPLETE."
