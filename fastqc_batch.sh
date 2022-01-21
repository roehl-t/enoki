#!/bin/sh

#  fastqc_batch.sh
#  
#
#  Created by Todd Osmundson and Thomas Roehl on 10/6/21.
#  
  
  # set up environment
  data_folder=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data
  data_folder_new=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/new_data
  data_folder_arr=(${data_folder} ${data_folder_new})
  
  fastqc_app=/Applications/FastQC.app/Contents/MacOS/fastqc
  fqtrim_app=/Applications/fqtrim-0.9.7/fqtrim
  HISAT2_dir=/Library/Apple/usr/bin:/Applications/hisat2-2.2.1
  trimmomatic_app=/Applications/Trimmomatic-0.38/trimmomatic-0.38.jar
  trimmomatic_adapter_dir=/Applications/Trimmomatic-0.38/adapters
  interleave_pairs=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_concat/fqtrim_output/interleave_pairs.py
  
  NCores=6
  
  export PATH="$PATH:${HISAT2_dir}:${trimmomatic_adapter_dir}"

cd ${data_folder}
  
  mkdir FastQC_reports

# generate FastQC reports from raw data
for seqfile in ./*.fastq.gz; do
${fastqc_app} -t ${NCores} -o ./FastQC_reports ${seqfile}
done

# Adapter and 3' quality trimming using trimmomatic
mkdir ./trimmed
cd trimmed
cp ${trimmomatic_adapter_dir}/TruSeq3-PE-2.fa ./TruSeq3-PE-2.fa
echo "Adapter and 3' quality trimming using Trimmomatic..."
for filename in ${data_folder}/*_R1_001.fastq.gz; do
     inBase=${filename##*/}
     outBase=${inBase%_R1_001.fastq.gz}
     echo $inBase
     echo $outBase

   java -jar ${trimmomatic_app} PE -threads ${NCores} -trimlog trimmomatic_trimlog.txt -basein $filename -baseout ${outBase}.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
done

# Remove low-complexity sequences using fqtrim
mkdir ./fqtrim_output

echo "Removing low-complexity sequences and poly-A/poly-T sequences using fqtrim..."
for file in ${data_folder}/trimmed/*_1P.fq
do
readpair=$(echo $file,${file/_1P/_2P} | tr -d '\n')
${fqtrim_app} -l 30 -p ${NCores} -D -o fqtrimmed.fq --outdir ${data_folder}/trimmed/fqtrim_output -r fqtrimlog.txt ${readpair}
done

# concatenate sequence runs 1 and 2
cd ${data_folder}
mkdir ./data_concat
for seqfile in /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/trimmed/fqtrim_output do
basename=${seqfile##*/}
cat $seqfile /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/new_data/trimmed/fqtrim_output/$basename > ./data_concat/$basename
done

# Generate quality reports using FastQC
cd ./data_concat
mkdir ./FastQC_reports
echo "Conducting FastQC quality check of processed reads..."

for seqfile in ./*.fqtrimmed.fq; do
${fastqc_app} -t ${NCores} -o ./FastQC_reports ${seqfile}
done

#
# incorporate unpaired sequences-----------------

# Remove low-complexity sequences using fqtrim
for i in ${!data_folder_arr[@]}; do
cd ${data_folder_arr[$i]}/trimmed
mkdir ./fqtrim_output_unpaired

echo "Removing low-complexity sequences and poly-A/poly-T sequences using fqtrim..."
for file in ${data_folder_arr[$i]}/trimmed/*U.fq
do
${fqtrim_app} -l 30 -p ${NCores} -D -o fqtrimmed.fq --outdir ${data_folder_arr[$i]}/trimmed/fqtrim_output_unpaired -r fqtrimlog.txt ${file}
done
done

# concatenate sequence runs 1 and 2
cd ${data_folder}
mkdir ./data_concat_unpaired
for seqfile in /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/trimmed/fqtrim_output_unpaired/*U.fqtrimmed.fq; do
basename=${seqfile##*/}
cat $seqfile /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/new_data/trimmed/fqtrim_output_unpaired/$basename > ./data_concat_unpaired/$basename
done

# Generate quality reports using FastQC
cd ./data_concat_unpaired
mkdir ./FastQC_reports
echo "Conducting FastQC quality check of processed reads..."

for seqfile in ./*.fqtrimmed.fq; do
${fastqc_app} -t ${NCores} -o ./FastQC_reports ${seqfile}
done

#------------------------------------------------
#

# rerun fqtrim on concatenated paired sequences (ignore pairs with fqtrim)
cd /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_concat/
mkdir ./fqtrim_output

echo "Removing low-complexity sequences and poly-A/poly-T sequences using fqtrim..."

for seqfile in /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_concat/*.fqtrimmed.fq; do
${fqtrim_app} -l 30 -p ${NCores} -D -o 2.fq --outdir /Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_concat/fqtrim_output -r fqtrimlog.txt ${seqfile}
done

# Generate quality reports using FastQC on reruns
cd ./fqtrim_output
mkdir ./FastQC_reports
echo "Conducting FastQC quality check of processed reads..."

for seqfile in ./*.fqtrimmed.2.fq; do
${fastqc_app} -t ${NCores} -o ./FastQC_reports ${seqfile}
done

# reorder paired sequences
cd ${data_folder}/data_concat/fqtrim_output/
echo "Performing sequence pair matching..."

for file in ./*_1P.fqtrimmed.2.fq; do
reverse=$(echo ${file/_1P/_2P} | tr -d '\n')
basename=${file##*/}
base=${basename%_*}
python3 ${interleave_pairs} ${file} ${reverse} ${base}
done

mkdir ./interleaved
mv ./*_out_*.fastq ./interleaved

# housekeeping: collect files in one place
cd ${data_folder}
mkdir ./data_qc_done

cp ${data_folder}/data_concat/fqtrim_output/interleaved/*_out_*.fastq ./data_qc_done
cp ${data_folder}/data_concat_unpaired/*.fqtrimmed.fq ./data_qc_done
