#!/usr/bin/env bash


bam_path="/mnt/datavault/SPP2002/analysis/exp28/bam"
scriptpath="bin"
input_annotation=""
output_rawreads=""
output_annotation=""
output_mapped=""
output_length=""

# Handling input
while getopts "h?i:o:b:s:r:m:l:" opt; do
    case "$opt" in
    h|\?)
        exit 0
        ;;
    b)  bam_path=$OPTARG
        ;;
    s)  scriptpath=$OPTARG
        ;;
    i)  input_annotation=$OPTARG
        ;;
    r)  output_rawreads=$OPTARG
        ;;
    m)  output_mapped=$OPTARG
        ;;
    l)  output_length=$OPTARG
        ;;
    o)  output_annotation=$OPTARG
        ;;
    esac
done

bam_list=()
for file in $bam_path/*.bam; do
    bam_list+=("$file")
done

echo "${bam_list[@]}"
echo "$output_annotation"
echo "$input_annotation"

python3 $scriptpath/total_mapped_reads.py -b "${bam_list[@]}" -m "$output_mapped" -l "$output_length"
python3 $scriptpath/call_featurecounts.py -b "${bam_list[@]}" -s 1 --with_O -o "$output_rawreads" -t 8 -a "$input_annotation"
python3 $scriptpath/map_reads_to_annotation.py -i "$output_rawreads" -a "$input_annotation" -o "$output_annotation"
