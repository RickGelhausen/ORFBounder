#!/usr/bin/env bash

# path="/mnt/datavault/SPP2002/analysis/campy_tt_analysis"
# scriptpath="/mnt/datavault/SPP2002/analysis/TTS_analysis"
# annotationpath="/mnt/datavault/SPP2002/analysis/exp28/annotation/"
# genomepath="/mnt/datavault/SPP2002/analysis/exp28/genomes/"
# experiments=("exp28" "exp28_no_tRNA")
# normalizations=("raw" "min" "mil")
# mappings=("fiveprime56" "fiveprime56" "fiveprime56" "fiveprime56" "fiveprime56-57" "fiveprime56-57" "fiveprimetracks" "fiveprimetracks" "threeprime56" "threeprime56" "threeprime29-31" "threeprime30" "threeprime30" "threeprime30" "threeprime30" "threeprimetracks")
# offsets=("-17" "-30" "-42" "-43" "-42" "-43" "-42" "-43" "13" "40" "13" "6" "13" "22" "40" "13")
# contrasts=("RIBO-A-1_TIS-A-1" "RIBO-A-1_TIS-A-3" "RIBO-A-1_TIS-A-5" "TIS-A-1_TIS-A-3" "TIS-A-1_TIS-A-5")
# bamfolder="/mnt/datavault/SPP2002/analysis/exp28/bam/"
# tmpfolder="/mnt/datavault/SPP2002/analysis/campy_tt_analysis/tmp/"

# Handling input
while getopts "h?p:s:a:g:e:m:n:o:c:b:t:" opt; do
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
    e)  experiments+=($OPTARG)
        ;;
    n)  normalizations+=($OPTARG)
        ;;
    m)  mappings=+=($OPTARG)
        ;;
    o)  offsets+=($OPTARG)
        ;;
    c)  contrasts+=($OPTARG)
        ;;
    b)  bamfolder=$OPTARG
        ;;
    t)  tmpfolder=$OPTARG
        ;;
    esac
done

mkdir -p $tmpfolder

echo "----------------------------------------------------"
for experiment in ${experiments[*]}; do
    for m_i in ${!mappings[@]}; do
        for norm in ${normalizations[*]}; do
            wigpath="$path/$experiment/${mappings[m_i]}/$norm"
            respath="$path/TTS_finder_results/TTS/$experiment/${mappings[m_i]}/${offsets[m_i]}/$norm"
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
                                                  --annotation_file=$annotationpath/annotation.gff --genome_file=$genomepath/genome.fa -o=$respath/$sample.$norm.csv --target_site=TTS --p_offset=${offsets[m_i]} \
                                                  --output_gff=$respath/$sample.$norm.gff --codon_interval_out=$respath/$sample.${norm}_codons.gff
            done

            infiles=()
            for entry in "$respath"/*.csv
            do
                infiles+=($entry)
            done
            python3 $scriptpath/merge_TTS.py -t ${infiles[@]} --contrasts ${contrasts[@]} -x $respath/${experiment}_${norm}_intermediate.xlsx

            python3 $scriptpath/detect_longest_potential_ORFs.py -i=$respath/${experiment}_${norm}_intermediate.xlsx -g=$genomepath/genome.fa -o=$respath/${experiment}_${norm}_overview.xlsx
            python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction both --target_site TTS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_both_ends.xlsx
            python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction up --target_site TTS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_upstream.xlsx
            python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction down --target_site TTS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_downstream.xlsx

            python3 $scriptpath/xlsx_to_gff.py -i $respath/${experiment}_${norm}_overview.xlsx -o $tmpfolder/${experiment}_${norm}_overview.gff
            python3 $scriptpath/xlsx_to_gff.py -i $respath/${experiment}_${norm}_overview_filtered_both_ends.xlsx -o $tmpfolder/${experiment}_${norm}_overview_filtered_both_ends.gff
            python3 $scriptpath/xlsx_to_gff.py -i $respath/${experiment}_${norm}_overview_filtered_upstream.xlsx -o $tmpfolder/${experiment}_${norm}_overview_filtered_upstream.gff
            python3 $scriptpath/xlsx_to_gff.py -i $respath/${experiment}_${norm}_overview_filtered_downstream.xlsx -o $tmpfolder/${experiment}_${norm}_overview_filtered_downstream.gff

            bash $scriptpath/get_readcount_gff.sh -i $tmpfolder/${experiment}_${norm}_overview.gff -r $tmpfolder/${experiment}_${norm}_overview_readcounts.raw -o $tmpfolder/${experiment}_${norm}_overview_readcounts.gff -b $bamfolder -m $tmpfolder/${experiment}_${norm}_overview_total_mapped.txt -l $tmpfolder/${experiment}_${norm}_overview_lengths.txt
            bash $scriptpath/get_readcount_gff.sh -i $tmpfolder/${experiment}_${norm}_overview_filtered_both_ends.gff -r $tmpfolder/${experiment}_${norm}_overview_filtered_both_ends_readcounts.raw -o $tmpfolder/${experiment}_${norm}_overview_filtered_both_ends_readcounts.gff -b $bamfolder -m $tmpfolder/${experiment}_${norm}_overview_filtered_both_ends_total_mapped.txt -l $tmpfolder/${experiment}_${norm}_overview_filtered_both_ends_lengths.txt
            bash $scriptpath/get_readcount_gff.sh -i $tmpfolder/${experiment}_${norm}_overview_filtered_upstream.gff -r $tmpfolder/${experiment}_${norm}_overview_filtered_upstream_readcounts.raw -o $tmpfolder/${experiment}_${norm}_overview_filtered_upstream_readcounts.gff -b $bamfolder -m $tmpfolder/${experiment}_${norm}_overview_filtered_upstream_total_mapped.txt -l $tmpfolder/${experiment}_${norm}_overview_filtered_upstream_lengths.txt
            bash $scriptpath/get_readcount_gff.sh -i $tmpfolder/${experiment}_${norm}_overview_filtered_downstream.gff -r $tmpfolder/${experiment}_${norm}_overview_filtered_downstream_readcounts.raw -o $tmpfolder/${experiment}_${norm}_overview_filtered_downstream_readcounts.gff -b $bamfolder -m $tmpfolder/${experiment}_${norm}_overview_filtered_downstream_total_mapped.txt -l $tmpfolder/${experiment}_${norm}_overview_filtered_downstream_lengths.txt

            python3 $scriptpath/calculate_expression.py -i $respath/${experiment}_${norm}_overview.xlsx -m $tmpfolder/${experiment}_${norm}_overview_total_mapped.txt -r $tmpfolder/${experiment}_${norm}_overview_readcounts.gff -o $respath/${experiment}_${norm}_overview_final.xlsx
            python3 $scriptpath/calculate_expression.py -i $respath/${experiment}_${norm}_overview_filtered_both_ends.xlsx -m $tmpfolder/${experiment}_${norm}_overview_filtered_both_ends_total_mapped.txt -r $tmpfolder/${experiment}_${norm}_overview_filtered_both_ends_readcounts.gff -o $respath/${experiment}_${norm}_overview_filtered_both_ends_final.xlsx
            python3 $scriptpath/calculate_expression.py -i $respath/${experiment}_${norm}_overview_filtered_upstream.xlsx -m $tmpfolder/${experiment}_${norm}_overview_filtered_upstream_total_mapped.txt -r $tmpfolder/${experiment}_${norm}_overview_filtered_upstream_readcounts.gff -o $respath/${experiment}_${norm}_overview_filtered_upstream_final.xlsx
            python3 $scriptpath/calculate_expression.py -i $respath/${experiment}_${norm}_overview_filtered_downstream.xlsx -m $tmpfolder/${experiment}_${norm}_overview_filtered_downstream_total_mapped.txt -r $tmpfolder/${experiment}_${norm}_overview_filtered_downstream_readcounts.gff -o $respath/${experiment}_${norm}_overview_filtered_downstream_final.xlsx
        done
    done
done
