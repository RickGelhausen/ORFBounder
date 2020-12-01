#!/usr/bin/env bash

scriptpath="bin"

read_count_threshold=5
# Handling input
while getopts "h?p:s:a:g:e:m:n:o:c:b:t:r:x:y:w:f" opt; do
    case "$opt" in
    h|\?)
        exit 0
        ;;
    p)  path=$OPTARG
        ;;
    w)  coveragepath=$OPTARG
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
    x)  start_codons+=("$OPTARG")
        ;;
    y)  stop_codons+=("$OPTARG")
        ;;
    f)  postfiltering=true
        ;;
    esac
done

if [ ${#start_codons[@]} -eq 0 ]; then
    start_codons=("ATG" "GTG" "TTG")
fi

if [ ${#stop_codons[@]} -eq 0 ]; then
    stop_codons=("TAG" "TAA" "TGA")
fi

mkdir -p $tmpfolder

echo "----------------------------------------------------"
for experiment in ${experiments[*]}; do
    for m_i in ${!mappings[@]}; do
        for norm in ${normalizations[*]}; do
            wigpath="$coveragepath/$experiment/${mappings[m_i]}/$norm"
            respath="$path/ORFBounder_results/TIS/$experiment/${mappings[m_i]}/${offsets[m_i]}/$norm"
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
                python3 $scriptpath/ORFBounder.py --fwd_file $coveragepath/$experiment/${mappings[m_i]}/$norm/$sample.$norm.forward.wig --rev_file $coveragepath/$experiment/${mappings[m_i]}/$norm/$sample.$norm.reverse.wig \
                                                  --annotation_file $annotationpath --genome_file $genomepath -o $respath/$sample.$norm.csv --target_site TIS --offset ${offsets[m_i]} \
                                                  --output_gff $respath/$sample.$norm.gff --codon_interval_out $respath/$sample.${norm}_codons.gff -c $readcountthreshold \
                                                  --start_codons ${start_codons[@]} --stop_codons ${stop_codons[@]}
            done

            infiles=()
            for entry in "$respath"/*.csv
            do
                infiles+=($entry)
            done

            python3 $scriptpath/merge_TTS.py -t ${infiles[@]} --contrasts ${contrasts[@]} -x $respath/${experiment}_${norm}_overview.xlsx
            python3 $scriptpath/xlsx_to_gff.py -i $respath/${experiment}_${norm}_overview.xlsx -o $tmpfolder/${experiment}_${norm}_overview.gff
            bash $scriptpath/get_readcount_gff.sh -s $scriptpath -i $tmpfolder/${experiment}_${norm}_overview.gff -r $tmpfolder/${experiment}_${norm}_overview_readcounts.raw -o $tmpfolder/${experiment}_${norm}_overview_readcounts.gff -b $bamfolder -m $tmpfolder/${experiment}_${norm}_overview_total_mapped.txt -l $tmpfolder/${experiment}_${norm}_overview_lengths.txt
            python3 $scriptpath/calculate_expression.py -i $respath/${experiment}_${norm}_overview.xlsx -m $tmpfolder/${experiment}_${norm}_overview_total_mapped.txt -r $tmpfolder/${experiment}_${norm}_overview_readcounts.gff -o $respath/${experiment}_${norm}_overview_final.xlsx

            if [ "$postfiltering" = true ] ; then
                python3 $scriptpath/post_filtering.py -i $respath/${experiment}_${norm}_overview_final.xlsx --size 25 --direction both --target_site TIS -a $annotationpath -o $respath/${experiment}_${norm}_overview_final_filtered_both_ends.xlsx
                python3 $scriptpath/post_filtering.py -i $respath/${experiment}_${norm}_overview_final.xlsx --size 25 --direction up --target_site TIS -a $annotationpath -o $respath/${experiment}_${norm}_overview_final_filtered_upstream.xlsx
                python3 $scriptpath/post_filtering.py -i $respath/${experiment}_${norm}_overview_final.xlsx --size 25 --direction down --target_site TIS -a $annotationpath -o $respath/${experiment}_${norm}_overview_final_filtered_downstream.xlsx
            fi
        done
    done
done
