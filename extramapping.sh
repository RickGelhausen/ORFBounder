
#samtools
#wigToBigWig
# Handling input
while getopts "h?b:r:o:s:g:" opt; do
    case "$opt" in
    h|\?)
        exit 0
        ;;
    b)  bamfolder=$OPTARG
        ;;
    r)  readlengths+=("$OPTARG")
        ;;
    o)  outputfolder=$OPTARG
        ;;
    s)  hriboscriptpath=$OPTARG
        ;;
    g)  genomesizes=$OPTARG
        ;;
    esac
done

normalization=("mil" "min" "raw")
mappings=("threeprimetracks" "fiveprimetracks")

bam_list=()
for file in $bamfolder/*.bam; do
    bam_list+=("$file")
done


for length in "${readlengths[@]}"; do
    mkdir -p $outputfolder/bam${length}

    boundaries=(${length//-/ })
    if [ ${#boundaries[@]} -eq 2 ]; then
        for bam in "${bam_list[@]}"; do
            bam_base="$(basename "$bam")"
            samtools index $bam
            samtools view -h $bam | awk -v left=${boundaries[0]} -v right=${boundaries[1]} 'length($10) >= left  && length($10) <= right || $1 ~ /^@/' | samtools view -bS - > $outputfolder/bam$length/noheader_$bam_base
            samtools view -H $bam > $outputfolder/bam$length/$bam_base.header
        	  samtools reheader $outputfolder/bam$length/$bam_base.header $outputfolder/bam$length/noheader_$bam_base > $outputfolder/bam$length/$bam_base
        	  samtools index $outputfolder/bam$length/$bam_base
            rm $outputfolder/bam$length/$bam_base.header $outputfolder/bam$length/noheader_$bam_base
        done
    else
        for bam in "${bam_list[@]}"; do
            bam_base="$(basename "$bam")"
            samtools index $bam
            samtools view -h $bam | awk -v bound=${boundaries[0]} 'length($10) == bound || $1 ~ /^@/' | samtools view -bS - > $outputfolder/bam$length/noheader_$bam_base
            samtools view -H $bam > $outputfolder/bam$length/$bam_base.header
            samtools reheader $outputfolder/bam$length/$bam_base.header $outputfolder/bam$length/noheader_$bam_base > $outputfolder/bam$length/$bam_base
            samtools index $outputfolder/bam$length/$bam_base
            rm $outputfolder/bam$length/$bam_base.header $outputfolder/bam$length/noheader_$bam_base
        done
    fi

    new_bam_list=()
    for file in $outputfolder/bam$length/*.bam; do
        new_bam_list+=("$file")
    done

    $hriboscriptpath/total_mapped_reads.py -b ${new_bam_list[@]} -m $outputfolder/bam$length/readcounts -l $outputfolder/bam$length/lengths

    mkdir -p $outputfolder/bam$length/wig/threeprimetracks; mkdir -p $outputfolder/bam$length/wig/threeprimetracks/raw; mkdir -p $outputfolder/bam$length/wig/threeprimetracks/mil; mkdir -p $outputfolder/bam$length/wig/threeprimetracks/min;
    for bam in "${new_bam_list[@]}"; do
            y=${bam%.bam}
            $hriboscriptpath/mapping.py --mapping_style last_base_only --bam_path $bam --wiggle_file_path $outputfolder/bam$length/wig/threeprimetracks/ --no_of_aligned_reads_file_path $outputfolder/bam$length/readcounts --library_name ${y##*/}
    done

    mkdir -p $outputfolder/bam$length/wig/fiveprimetracks; mkdir -p $outputfolder/bam$length/wig/fiveprimetracks/raw; mkdir -p $outputfolder/bam$length/wig/fiveprimetracks/mil; mkdir -p $outputfolder/bam$length/wig/fiveprimetracks/min;
    for bam in "${new_bam_list[@]}"; do
            y=${bam%.bam}
            $hriboscriptpath/mapping.py --mapping_style first_base_only --bam_path $bam --wiggle_file_path $outputfolder/bam$length/wig/fiveprimetracks/ --no_of_aligned_reads_file_path $outputfolder/bam$length/readcounts --library_name ${y##*/}
    done


    for norm in ${normalization[@]}; do
        for mapping in ${mappings[@]}; do
            for entry in $outputfolder/bam$length/wig/$mapping/$norm/*.wig; do
                y=${entry%.wig}
                wigToBigWig $entry $genomesizes $outputfolder/bam$length/wig/$mapping/$norm/${y##*/}.bw
            done
        done
    done
done
