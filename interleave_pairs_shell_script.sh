#!/bin/sh

#  interleave_pairs_shell_script.sh
#  
#
#  Created by Thomas Roehl on 11/9/21.
#  

  
# set up environment
data_folder=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/data/data_concat/fqtrim_output/
interleave_pairs=/Volumes/RAID_5_data_array/Todd/Thomas_Roehl_RNASeq/scripts/interleave_pairs.py

# reorder paired sequences using interleave_pairs.py
cd ${data_folder}
echo "Performing sequence pair matching..."

# [make sure to change file match criteria to fit your naming patterns]
for file in ./*_1P.fqtrimmed.2.fq; do
reverse=$(echo ${file/_1P/_2P} | tr -d '\n')
basename=${file##*/}
base=${basename%_*}
echo ${base}
python3 ${interleave_pairs} ${file} ${reverse} ${base}
done

# move processed files to a new folder
# [if your original files have "_out_" in their file names, change the match criteria below]
mkdir ./interleaved
mv ./*_out_*.fastq ./interleaved
