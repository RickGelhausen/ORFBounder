#!/usr/bin/env bash

# bash TTS_analysis/tts_finder_analysis_TIS.sh -p /mnt/datavault/SPP2002/analysis/campy_tt_analysis -s /mnt/datavault/SPP2002/analysis/TTS_analysis -a /mnt/datavault/SPP2002/analysis/haloferax_TIS/annotation/annotation.gff -g /mnt/datavault/SPP2002/analysis/haloferax_TIS/genomes/genome.fa -e haloferax_TIS -m fiveprimetracks -m threeprimetracks -n raw -o 0 -o 31 -c RIBO-A-1_TIS-A-1 -c RIBO-A-2_TIS-A-1 -b /mnt/datavault/SPP2002/analysis/haloferax_TIS/bam/ -t /mnt/datavault/SPP2002/analysis/campy_tt_analysis/tmp/ -r 5

# path="/mnt/datavault/SPP2002/analysis/campy_tt_analysis"
# scriptpath="/mnt/datavault/SPP2002/analysis/TTS_analysis"
# annotationpath="/mnt/datavault/SPP2002/analysis/haloferax_TIS/annotation/annotation.gff"
# genomepath="/mnt/datavault/SPP2002/analysis/haloferax_TIS/genomes/genome.fa"
# experiments=("haloferax_TIS")
# normalizations=("raw" "min" "mil")
# mappings=("threeprimetracks" "fiveprimetracks")
# offsets=("31" "0")
# bamfolder="/mnt/datavault/SPP2002/analysis/haloferax_TIS/bam/"
# tmpfolder="/mnt/datavault/SPP2002/analysis/campy_tt_analysis/tmp/"

read_count_threshold=5
# Handling input
while getopts "h?p:s:a:g:e:m:n:o:c:b:t:r:" opt; do
    case "$opt" in
    h|\?)
        exit 0
        ;;
    p)  path=$OPTARG
        ;;
    s)  scriptpath=$OPTARG
        ;;
    a)  annotationpath=$OPTARG
        ;;
    g)  genomepath=$OPTARG
        ;;
    e)  experiments+=("$OPTARG")
        ;;
    n)  normalizations+=("$OPTARG")
        ;;
    m)  mappings+=("$OPTARG")
        ;;
    o)  offsets+=("$OPTARG")
        ;;
    c)  contrasts+=("$OPTARG")
        ;;
    b)  bamfolder=$OPTARG
        ;;
    t)  tmpfolder=$OPTARG
        ;;
    r)  readcountthreshold=$OPTARG
        ;;
    esac
done

mkdir -p $tmpfolder

echo "----------------------------------------------------"
for experiment in ${experiments[*]}; do
    for m_i in ${!mappings[@]}; do
        for norm in ${normalizations[*]}; do
            wigpath="$path/$experiment/${mappings[m_i]}/$norm"
            respath="$path/TTS_finder_results/TIS/$experiment/${mappings[m_i]}/${offsets[m_i]}/$norm"
            mkdir -p $respath
            mkdir -p $wigpath

            echo $wigpath
            echo $respath

            prefix_list=()
            for file in $wigpath/*.wig; do
                file_base="$(basename "$file")"
                IFS="." read -r -a prefix_arr <<< "$file_base"

                prefix=${prefix_arr[0]}
                if [[ $prefix != *"RNATIS-"* && $prefix != *"RNA-"* ]]; then
                  #echo "$prefix"
                  prefix_list+=("$prefix")
                fi
            done
            echo "$experiment ${mappings[m_i]} $norm"
            uniq_prefix=($(printf "%s\n" "${prefix_list[@]}" | sort -u | tr '\n' ' '))

            for sample in ${uniq_prefix[@]}; do
                python3 $scriptpath/TTS_finder.py --fwd_file=$path/$experiment/${mappings[m_i]}/$norm/$sample.$norm.forward.wig --rev_file=$path/$experiment/${mappings[m_i]}/$norm/$sample.$norm.reverse.wig \
                                                  --annotation_file=$annotationpath --genome_file=$genomepath -o=$respath/$sample.$norm.csv --target_site=TIS --p_offset=${offsets[m_i]} \
                                                  --output_gff=$respath/$sample.$norm.gff --codon_interval_out=$respath/$sample.${norm}_codons.gff -c $readcountthreshold
            done

            infiles=()
            for entry in "$respath"/*.csv
            do
                infiles+=($entry)
            done

            python3 $scriptpath/merge_TTS.py -t ${infiles[@]} --contrasts ${contrasts[@]} -x $respath/${experiment}_${norm}_overview.xlsx
            python3 $scriptpath/xlsx_to_gff.py -i $respath/${experiment}_${norm}_overview.xlsx -o $tmpfolder/${experiment}_${norm}_overview.gff
            bash $scriptpath/get_readcount_gff.sh -i $tmpfolder/${experiment}_${norm}_overview.gff -r $tmpfolder/${experiment}_${norm}_overview_readcounts.raw -o $tmpfolder/${experiment}_${norm}_overview_readcounts.gff -b $bamfolder -m $tmpfolder/${experiment}_${norm}_overview_total_mapped.txt -l $tmpfolder/${experiment}_${norm}_overview_lengths.txt
            python3 $scriptpath/calculate_expression.py -i $respath/${experiment}_${norm}_overview.xlsx -m $tmpfolder/${experiment}_${norm}_overview_total_mapped.txt -r $tmpfolder/${experiment}_${norm}_overview_readcounts.gff -o $respath/${experiment}_${norm}_overview_final.xlsx

            python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview_final.xlsx --size 25 --direction both --target_site TTS -a $annotationpath -o $respath/${experiment}_${norm}_overview_final_filtered_both_ends.xlsx
            python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview_final.xlsx --size 25 --direction up --target_site TTS -a $annotationpath -o $respath/${experiment}_${norm}_overview_final_filtered_upstream.xlsx
            python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview_final.xlsx --size 25 --direction down --target_site TTS -a $annotationpath -o $respath/${experiment}_${norm}_overview_final_filtered_downstream.xlsx

        done
    done
done
