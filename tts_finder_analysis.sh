#!/usr/bin/env bash

path="/mnt/datavault/SPP2002/analysis/campy_tt_analysis"
scriptpath="/mnt/datavault/SPP2002/analysis/TTS_analysis"
annotationpath="/mnt/datavault/SPP2002/analysis/exp3/annotation/"
genomepath="/mnt/datavault/SPP2002/analysis/exp3/genomes/"
targetsite=("TTS")
experiments=("exp28", "exp28_no_tRNA")
mappings=("fiveprime56" "fiveprime56" "fiveprime56" "fiveprime56" "fiveprime56-57" "fiveprime56-57" "fiveprimetracks" "fiveprimetracks" "threeprime56" "threeprime56" "threeprime29-31" "threeprime30" "threeprime30" "threeprime30" "threeprime30" "threeprimetracks")
#mappings=("fiveprimetracks" "fiveprimetracks" "threeprimetracks")

normalizations=("raw" "min" "mil")

offsets=("-17" "-30" "-42" "-43" "-42" "-43" "-42" "-43" "13" "40" "13" "6" "13" "22" "40" "13")
#offsets=("-42" "-43" "13")
contrasts=("RIBO-A-1_TIS-A-1" "RIBO-A-1_TIS-A-3" "RIBO-A-1_TIS-A-5" "TIS-A-1_TIS-A-3" "TIS-A-1_TIS-A-5")

echo "----------------------------------------------------"
for experiment in ${experiments[*]}; do
    for method in ${targetsite[*]}; do
        for m_i in ${!mappings[@]}; do
            for norm in ${normalizations[*]}; do
                wigpath="$path/$experiment/${mappings[m_i]}/$norm"
                respath="$path/TTS_finder_results/$method/$experiment/${mappings[m_i]}/${offsets[m_i]}/$norm"
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
                                                      --annotation_file=$annotationpath/annotation.gff --genome_file=$genomepath/genome.fa -o=$respath/$sample.$norm.csv --target_site=$method --p_offset=${offsets[m_i]} \
                                                      --output_gff=$respath/$sample.$norm.gff --codon_interval_out=$respath/$sample.${norm}_codons.gff
                done

                infiles=()
                for entry in "$respath"/*.csv
                do
                    infiles+=($entry)
                done
                python3 $scriptpath/merge_TTS.py -t ${infiles[@]} --contrasts ${contrasts[@]} -x $respath/${experiment}_${norm}_overview.xlsx
                python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction both --target_site TTS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_both_ends.xlsx
                python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction up --target_site TTS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_upstream.xlsx
                python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction down --target_site TTS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_downstream.xlsx

            done
        done
    done
done



path="/mnt/datavault/SPP2002/analysis/campy_tt_analysis"
scriptpath="/mnt/datavault/SPP2002/analysis/scripts"
annotationpath="/mnt/datavault/SPP2002/analysis/exp3/annotation/"
genomepath="/mnt/datavault/SPP2002/analysis/exp3/genomes/"
targetsite=("TIS")
experiments=("exp28_no_tRNA")
mappings=("threeprimetracks" "threeprime30-32" "threeprime32" "fiveprimetracks" "fiveprime30-32" "fiveprime32")
normalizations=("raw" "min" "mil")
offsets=("16" "16" "16" "-16" "-16" "-16")


echo "----------------------------------------------------"
for experiment in ${experiments[*]}; do
    for method in ${targetsite[*]}; do
        for m_i in ${!mappings[@]}; do
            for norm in ${normalizations[*]}; do
                wigpath="$path/$experiment/${mappings[m_i]}/$norm"
                respath="$path/TTS_finder_results/$method/$experiment/${mappings[m_i]}/${offsets[m_i]}/$norm"
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
                                                      --annotation_file=$annotationpath/annotation.gff --genome_file=$genomepath/genome.fa -o=$respath/$sample.$norm.csv --target_site=$method --p_offset=${offsets[m_i]} \
                                                      --output_gff=$respath/$sample.$norm.gff --codon_interval_out=$respath/$sample.${norm}_codons.gff
                done

                infiles=()
                for entry in "$respath"/*.csv
                do
                    infiles+=($entry)
                done
                python3 $scriptpath/merge_TTS.py -t ${infiles[@]} --contrasts RIBO-A-1_TIS-A-1 RIBO-A-1_TIS-A-3 RIBO-A-1_TIS-A-5 TIS-A-1_TIS-A-3 TIS-A-1_TIS-A-5 -x $respath/${experiment}_${norm}_overview.xlsx
                python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction both --target_site TIS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_both_ends.xlsx
                python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction up --target_site TIS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_upstream.xlsx
                python3 $scriptpath/post_filter_tts.py -i $respath/${experiment}_${norm}_overview.xlsx --size 25 --direction down --target_site TIS -a $annotationpath/annotation.gff -o $respath/${experiment}_${norm}_overview_filtered_downstream.xlsx
            done
        done
    done
done
