# TTS_analysis

Detection of potential start/stop codons based on Translation Initiation Site (TIS) or Translatation Termination Site (TTS) peaks, using a similar concept than the [RETscript for TIS](https://www.sciencedirect.com/science/article/pii/S1097276519301078).

These scripts were created to be used with the metagene-profiling and coverage (.wig) files created by the [HRIBO workflow](https://github.com/RickGelhausen/HRIBO) [[1]](#1).
Nevertheless, the TTS/TIS_finder scripts can also be used with any other coverage files and metagene-profiling tools.

# Requirements
## Required packages
The scripts used in the analysis exclusively require python3.
Packages required are:
* pandas
* numpy
* pysam
* biopython
* subread
* interlap
* xlrd
* xlsxwriter

All required packages are easily retrievable via conda.

```
conda create -n "TTS_finder" -c bioconda -c conda-forge pysam numpy pandas biopython subread interlap xlrd xlsxwriter
conda activate TTS_finder
```

## Required files
If you used the HRIBO workflow, you will have all files required to run the analysis. 
It is important to note that you can also run the analysis partially, by manually calling the individual scripts provided in this repository (e.g. if you do not require expression values, you do not require bam files).

* `wig files:` The wig files for the desired mappings/normalizations. For easy usage, these should be in the HRIBO notation. `path/experiment/mapping/normalization/|method|-|condition|-|replicate|.normalization.forward.wig`. (e.g `/path_to_user/exp1/threeprimetracks/min/TIS-A-1.min.forward.wig`)
Wig files must be split into two individual files, one for each strand. (forward, reverse) A folder can contain wig files of the according mapping and normalization for multiple samples. The scripts will be run on all files and bundled into one result file.
If you do not have .wig files from `HRIBO`, either create an according folder structure or write your own script tailored to your data, using the scripts provided in this repository. Explanation for each script are provided in the [scripts section](#Scripts).

* `bam files`: `HRIBO` provides `.bam` files containing all read counts. These should be named using the `|method|-|condition|-|replicate|.bam` naming scheme. `|method|` is either `RIBO`, `RNA`, `TIS`, `RNATIS`, `TTS`, `RNATTS` . `|condition|` can be any string and `|replicate|` can be any integer.

* `genome file`: a genome file in fasta format for the analysed organism.
* `annotation file`: an annotation file in .gff3 format for the analysed organism. (Tested using annotation files from NCBI)


# Analysis 
To run the scripts, coverage files in `.wig format` are required. We generated the coverage files using `HRIBO` [[1]](#1). `HRIBO` generates coverage files with centered, threeprime, fiveprime and global mappings. Additionally, it provides raw and normalized coverage files for each of the methods (raw, min, mil). 
Then metagene-profiling is performed on all coverage files. 
From the resulting metagene-profiling results, we determined the best p-site-offsets for interesting mapping / normalization combinations. These were then analysed with the `TTS_finder` scripts.
`HRIBO` additionally provides `.bam files` containing all read lengths. If you want to investigate specific read lengths, an example script for filtering .bam files for different read-lengths is provided in this repository.

The analysis is done in multiple steps:

1. First, all potential stop (or start) codons are collected for the TTS (TIS) analysis. The stop codon position is then expanded into an interval of 7nt. These intervals are collected and saved in an InterLap object. The p/a-site offsets are added/substracted from these intervals, in order to ensure that the correct regions are investigated. (The codon intervals are also written to a .gff file. These can be loaded and investigated in a genome-browser.)

2. Next, the requested .wig files are read and if a position passes a given read_count_threshold (default 5), all codon intervals overlapping with the given position are retrieved and their peak_height is incremented by the read count of the position detected. This matches the coverage peaks with given stop(start) codons. 

3. Then, we iterate over all potential codons that have a peak attributed to them. For each stop(codon) the next in-frame start(stop) codon is searched and formed into an ORF prediction. These ORFs are collected and written into a `.gff` and a `.csv` (table) file.  Additionally, .gff files for each gene_type are generated for easier investigation in a genome_browser.

4. If the script was used on different RIBO-seq, RNA-seq and TIS or TTS samples, all result tables are bundled into one big excel file, by combining ORFs that have been predicted for multiple samples into one row, providing the peak_height information for all involved samples. The excel file contains a lot of additional information for each predicted ORF (e.g. gene_type, start, stop, strand, locus_tag, codon_count, peak_height, 15nt upstream of the start, nucleotide sequence, amino acid sequence, etc...). Contrasts can be given in form of a list of file prefixes (e.g RIBO-A-1_TIS_A-1, RIBO-A-2_TIS_A-1). This will add additional columns with log2foldchange for the given prefix combinations.

5. (TTS_only) For the TTS predictions, it is hard to find the best start codon matching the predicted stop codon, as multiple start codons can be present in-frame upstream of the predicted stop codon. The method used in the TTS_finder script, finds the shortest possible ORF, by choosing the first in-frame start-codon. In this step, the longest possible ORF is added to the results for a given predicted stop. First, the first in-frame stop codon upstream of the current predicted stop-codon is searched, then the first start-codon downstream of the upstream stop-codon is chosen. This ensures that the detected start-codon is the furthest possible in-frame start-codon that ensures that no additional in-frame stop-codon is between the predicted stop-codon and the attributed start-codon. This provides us with the longest possible ORF.
For TIS predictions this is not necessary, as we start from the predicted start-codon and look for the first in-frame stop-codon.

6. Next, the the excel files are filtered and split into three different files. The idea is to better distinguish un-annotated ORFs that were predicted. To do this, we ensure that no annotated stop is within 25nt of a predicted stop, otherwise it is filtered out. This is done once for the upstream direction, the downstream direction and both, resulting in 3 different files. This is just an additional method that might help easier manual investigation of potentially new ORFs. The original table is also valid, and contains all information.

7. In a final step, expression information is added to all tables for all predicted ORF intervals. This includes both read per kilobase million (RPKM) values and translational efficiency (TE) values. To do this, the read counts are collected using subread-featureCounts. These readcounts are then used in order to calculate both the RPKM and the TE for every sample.

 :warning: **IMPORTANT:** The scripts are written to be compatible with the HRIBO workflow, all samples must be in the form `|method|-|condition|-|replicate|`. `|method|` is either `RIBO`, `RNA`, `TIS`, `RNATIS`, `TTS`, `RNATTS`. `|condition|` can be any string and `|replicate|` can be any integer.

The chosen thresholds, offsets and coverage mappings can change for each organism, therefore it is advised to investigate the data first to ensure that the right parameters are chosen. :warning:


# Running the analysis scripts
The analysis is made up of multiple python3 and bash scripts. If all data is collected as described in the [required files section](#Required-files), running the script will be straight-forward.
Simply run either `tts_finder_analysis_TTS.sh` or `tts_finder_analysis_TIS.sh` depending on the site that is to be analysed.

The following commandline arguments are required:
| Name                | Command Line Argument | Description                                                                                                           |
|---------------------|-----------------------|-----------------------------------------------------------------------------------------------------------------------|
| path                | -p                    | Path where the result folder will be placed.                                                                          |
| coveragepath        | -p                    | Path to the folder containing the experiments that will be analysed                                                   |
| scriptpath          | -s                    | Path to the TTS_analysis folder, in which all scripts necessary for computation are located.                          |
| annotationpath      | -a                    | The annotation file for the organism that is analysed (`.gff3` format)                                                |
| genomepath          | -g                    | The genome file for the organism that is analysed (`.fasta` format)                                                   |
| experiments         | -e                    | The experiment to be analysed (if more than one, use this option multiple times e.g `-e exp1 -e exp2 ...`)            |
| mappings            | -m                    | The mappings to be used for analysis (e.g threeprime, fiveprime, etc...) ( any subfolder under experiment). If more than one, use this option multiple times (e.g. `-m threeprime -m threeprime30 -m fiveprime ...`)       |    
| offsets             | -o                    | The p-site offset used for each of the mappings. If you use multiple mappings, ensure that you use the same amount of offsets. (e.g. `-m threeprime -m fiveprime`, `-o 8 -o 10` means that the 3' files have an offset of 8 and the 5' files have an offset of 10).                          |
| normalizations      | -n                    | The normalizations to be used (e.g. raw, mil, min) (any subfolder under mapping).                                     |
| contrasts           | -c                    | The contrasts used for the experiment. If you want log2FC for certain peak-heights in a table you can use this option to indicate which samples should be compared (e.g. RIBO-A-1_TIS-A-1)                                                                                                       |
| bamfolder           | -b                    | The path to the bamfiles, ensure that each bamfile has an according index file. If not, use `samtools index |bamfile|` for all files missing the index. If installed you can also use `parallel`[[2]](#2),  `parallel  samtools index ::: *.bam` to run it on all bam files.              |
| tmpfolder           | -t                    | The folder where temporary files will be dumped.                                                                      |
| readcountthreshold  | -r                    | The readcount threshold used in the analysis, if a position has less than this amount of reads it is ignored (>0)     |
| start_codons        | -x                    | A list of start_codons. (Default: ATG, GTG, TTG) (e.g `-x ATG -x GTG -x TTG ...`)                                                  |
| stop_codons         | -y                    | A list of stop_codons. (Default: TAG, TAA, TGA) (e.g `-y TAG -y TAA -y TGA ...`)                                                      |


If you have your own data, you can run the scripts individually, each of them is described in the [scripts section](#Scripts) below.

# Scripts
This section contains short descriptions of each of the scripts (in execution order) and the commandline parameters.
* **TTS_finder.py:** is the main script which uses annotation, genome and wig files to detect potential ORFs using TIS or TTS read coverage peaks. 

| Name                 | Command Line Argument | Description                                                                                                           |
|----------------------|-----------------------|-----------------------------------------------------------------------------------------------------------------------|
| fwd_file             | -f                    | The forward wig file used for the analysis.                                                                           |
| rev_file             | -r                    | The reverse wig file used for the analysis, complementary to the forward wig file.                                    |
| annotation_file      | -a                    | The annotation file for the organism that is analysed (`.gff3` format)                                                |
| genome_file          | -g                    | The genome file for the organism that is analysed (`.fasta` format)                                                   |
| start_codons         | --start_codons        | A space-seperated list of start_codons. (Default: ATG, GTG, TTG)                                                      |
| stop_codons          | --stop_codons         | A space-seperated list of stop_codons. (Default: TAG, TAA, TGA)                                                       |    
| p_offset             | --p_offset            | The p-site offset to be used for the current wig files.                                                               |
| output_gff           | --output_gff          | The output folder for the .gff files for genome browser inspection of the result.                                     |
| target_site          | --target_site         | The site you are interested in (TIS / TTS)                                                                            |
| read_count_threshold | -c                    | The readcount threshold used in the analysis, if a position has less than this amount of reads it is ignored (>0)     |
| codon_interval_out   | -codon_interval_out   | The output .gff file for the codon intervals which are used to test overlap with a potential codon.                   |
| output_file          | -o                    | The output .csv file for further processing in the other included scripts.                                            |

* **merge_TTS.py:** merges the .csv files resulting from the `TTS_finder.py` script, for different samples (RIBO-A-1, RIBO-A-2, TIS-A-1, etc...). The merged results are infused with additional information including log2 fold-changes for the different desired contrasts, nucleotide and amino-acid sequences. The resulting information is then written to an excel output file (.xlsx).

| Name                 | Command Line Argument | Description                                                                                                   |
|----------------------|-----------------------|---------------------------------------------------------------------------------------------------------------|
| tables               | -t/--tables           | A list of .csv tables resulting from `TTS_finder.py` that will be merged                                      |
| contrasts            | --contrasts           | The contrasts used for the experiment. If you want log2FC for certain peak-heights in the table you can use this option to indicate which samples should be compared (e.g. RIBO-A-1_TIS-A-1)|
| xlsx                 | -x                    | The output excel file                                                                                         |

* **detect_longest_potential_ORFs.py:** (optional) is a script only useful for TTS predictions, it detects the longest ORF in addition to the shortest ORF. This means that for a given predicted stop, the next upstream start codon is chosen (shortest ORF) and farthest upstream start codon is chosen, that still ensures that no additional in-frame stop-codon is between the start codon and the predicted stop codon.

| Name                 | Command Line Argument | Description                                                                                                   |
|----------------------|-----------------------|---------------------------------------------------------------------------------------------------------------|
| input_xlsx           | -i                    | An input xlsx file resulting from the `merge_TTS.py` script                                                   |
| genome_file          | -g                    | The genome file for the organism that is analysed (`.fasta` format)                                           |
| output_xlsx          | -o                    | The output excel file                                                                                         |
| start_codons         | --start_codons        | A space-seperated list of start_codons. (Default: ATG, GTG, TTG)                                              |
| stop_codons          | --stop_codons         | A space-seperated list of stop_codons. (Default: TAG, TAA, TGA)                                               |            

* **xlsx_to_gff.py:** creates simple .gff files from an .xlsx table to be used to run featureCounts. This will add read counts for each entry in the .gff file for each sample.

| Name                 | Command Line Argument | Description                                                                                                   |
|----------------------|-----------------------|---------------------------------------------------------------------------------------------------------------|
| input_xlsx           | -i                    | An input excel file                                                                                           |
| output_gff           | -o                    | The output annotation file (.gff)                                                                             |


* **get_readcount_gff.sh:** creates the .gff files containing readcounts for all samples. This bash script uses the additional scripts `call_featurrecounts.py`, `total_mapped_reads.py` and `map_read_to_annotation.py`. 

| Name                 | Command Line Argument | Description                                                                                                   |
|----------------------|-----------------------|---------------------------------------------------------------------------------------------------------------|
| bam_path             | -b                    | Path to the bam file folder. Bam files for the samples used in the analysis.                                  |
| input_annotation     | -i                    | The input annotation created by `xlsx_to_gff.py`                                                              |
| output_rawreads      | -r                    | The output file for the raw files (temporary file)                                                            |
| output_mapped        | -m                    | The output file for the total mapped reads for each plasmid/chromosome.                                       | 
| output_length        | -l                    | The output file for the average read lengths for each plasmid/chromosome. (currently not required.)           |   
| output_annotation    | -o                    | The output annotation containing the read counts for each entry.                                              |

* **calculate_expression.py:** calculates expression values (TE and RPKM) for each sample used and adds them to the final table. For TTS predictions, it uses both the short and the long ORF.

| Name                 | Command Line Argument | Description                                                                                                   |
|----------------------|-----------------------|---------------------------------------------------------------------------------------------------------------|
| input_xlsx           | -i                    | An input excel file that will be expended by adding translational efficiency and RPKM values.                 |
| mapped_reads         | -m                    | A file containing the total number of mapped reads for each sample and chromosome/plasmid, created by `get_readcount_gff.sh` |
| read_counts          | -r                    | A file containing all read counts for each ORF and each sample, created by `get_readcount_gff.sh`             |
| output_xlsx          | -o                    | The filtered output excel file                                                                                |  

* **post_filter_tts.py:** does an additional post-filtering step on the final .xlsx file, ensuring that no annotated start/stop codons are close to the predicted ORF boundaries. This might help narrowing down results to find good candidates for closer experimental inspection.

| Name                 | Command Line Argument | Description                                                                                                   |
|----------------------|-----------------------|---------------------------------------------------------------------------------------------------------------|
| input_xlsx           | -i                    | An input xlsx file to be filtered                                                                             |
| annotation_gff       | -a                    | The annotation file for the organism that is analysed (`.gff3` format)                                        |
| target_site          | --target_site         | The site that is currently analysed (TIS / TTS)                                                               |
| direction            | --direction           | The direction from the (start/stop) that will be filtered. (up/down/both)                                     |
| size                 | --size                | The size of the interval used for filtering in the current direction (if both is selected the interval spans over 2xsize | 
| output_xlsx          | -o                    | The filtered output excel file                                                                                |  

## Additional scripts
* **excel_utils.py:** is a python library script used for the `calculate_expression.py` script.
* **call_featurecounts.py:** calls `subread featureCount` to count the number of reads overlapping with a certain entry of an annotation file for a given sample.

| Name                 | Command Line Argument | Description                                                                                                           |
|----------------------|-----------------------|-----------------------------------------------------------------------------------------------------------------------|
| bamfiles             | -b                    | Bam files of the samples used in the analysis, used for read-counting.                                                |
| strandness           | -s                    | Whether strandedness should be considered (Default: 1 stranded)                                                       |
| with_O               | --with_O              | Activate `subread featureCount` command O. (Assign reads to all overlapping meta-features.)                           |
| with_M               | --with_M              | Activate `subread featureCount` command M. (Multi-mapping reads will be counted)                                      |
| fraction             | --fraction            | Assign fractional counts to features.                                                                                 |
| annotation           | -a                    | The annotation file that will be processed with `subread featureCount`.                                               |
| threads              | -t                    | Number of threads used by `subread featureCount`.                                                                     |
| output               | -o                    | The output annotation file with the read counts.                                                                      |
| for_diff_expr        | --for_diff_expr       | not required                                                                                                          |

* **map_reads_to_annotation.py:** map the read-counts calculated by `subread featureCount` back to the original annotation file.

| Name                 | Command Line Argument | Description                                                                                                           |
|----------------------|-----------------------|-----------------------------------------------------------------------------------------------------------------------|
| input                | -i                    | The raw file output by `subread featureCount`.                                                                        |
| annotation           | -a                    | Original annotation used for `subread featureCount`.                                                                  |
| output               | -o                    | The output gff file with read counts.                                                                                 |

* **total_mapped_reads.py:** counts the total number of mapped reads for a list of bam files. Additionally, determine the average read length.

| Name                 | Command Line Argument | Description                                                                                                           |
|----------------------|-----------------------|-----------------------------------------------------------------------------------------------------------------------|
| bamfiles             | -b                    | Bam files of the samples used in the analysis, used for counting the total number of mapped reads.                    |
| out_mapped           | -m                    | The output file with total number of mapped reads.                                                                    |
| out_length           | -l                    | The output file with the average read length.                                                                         |

## Extramapping script

If you want to analyse specific read-lengts, you will require both filtered wig and bam files.
In facilitate the creation of these files, we added the bash script `extramapping.sh`, which automatically generates the required files.

This is a script that is seperate from the analysis and requires different dependencies.
All required packages are easily retrievable via conda.

```
conda create -n "extramapping" -c bioconda -c conda-forge pandas samtools ucsc-wigtobigwig
conda activate extramapping
```

In addition, it also requires some scripts from `HRIBO` in order to do the mapping. If you do not have an `HRIBO` installation, you can simply download the `mapping.py` script from [GitHub repository](https://github.com/RickGelhausen/HRIBO/blob/master/scripts/mapping.py). (You will also need the `total_mapped_reads.py` script, but this is present in the TTS_analysis folder too. Ensure that both scripts are in the same folder.)

| Name                 | Command Line Argument | Description                                                                                                           |
|----------------------|-----------------------|-----------------------------------------------------------------------------------------------------------------------|
| bamfolder            | -b                    | Bam files of the samples used in the analysis, used for counting the total number of mapped reads.                    |
| readlengths          | -r                    | A list of read length for which you want new bam/wig files (intervals are allowed) (e.g. `-r 32 -r 29-31` etc... an interval will create one bam files containing only reads from 29 to 31 in length.                                                                                               |
| outputfolder         | -o                    | The output folder containing all .bam and .wig files.                                                                 |
| hriboscriptpath      | -s                    | The path to the `HRIBO scripts` folder. Must contain `total_mapped_reads.py` and  `mapping.py`.                       |
| genomesizes          | -g                    | A file with the genome sizes. This can be created using samtools: `samtools faidx genome.fa; cut -f1,2 {input[0]} > genomes/sizes.genome;`                                                                                                                                                 |


## References
<a id="1">[1]</a> 
Gelhausen, R. (2020).
HRIBO - High-throughput analysis of bacterial ribosome profiling data 
([BioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.27.046219v1))

<a id="2">[2]</a> 
Tange, O. (2011).
[GNU Parallel](http://www.gnu.org/software/parallel/) - The Command-Line Power Tool
[http://dx.doi.org/10.5281/zenodo.16303](http://dx.doi.org/10.5281/zenodo.16303)
