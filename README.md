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
If you do not have .wig files from `HRIBO`, either create an according folder structure or write your own script tailored to your data, using the scripts provided in this repository. Explanation for each script are provided in the [scripts section](#Scripts).

* `bam files`: `HRIBO` provides `.bam` files containing all read counts. These should be named using the `|method|-|condition|-|replicate|.bam` naming scheme. `|method|` is either `RIBO`, `RNA`, `TIS`, `RNATIS`. `TTS` and `RNATTS` will be supported soon, until then we suggest labeling `TTS` files `TIS`. `|condition|` can be any string and `|replicate|` can be any integer.

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

3. Then, we iterate over all potential codons that have a peak attributed to them. For each stop(codon) the next in-frame start(stop) codon is searched and formed into an ORF prediction. These ORFs are collected and written into a `.gff` and a `.xlsx` (excel table) file. The excel file contains a lot of additional information for each predicted ORF (e.g. gene_type, start, stop, strand, locus_tag, codon_count, peak_height, 15nt upstream of the start, nucleotide sequence, amino acid sequence, etc...). Additionally, .gff files for each gene_type are generated for easier investigation in a genome_browser.

4. If the script was used on different RIBO-seq, RNA-seq and TIS or TTS samples, all result tables are bundled into one big excel file, by combining ORF predicted for multiple samples into one row, providing the peak_height information for all involved samples. Contrasts can be given in form of a list of file prefixes (e.g RIBO-A-1_TIS_A-1, RIBO-A-2_TIS_A-1). This will add additional columns with log2foldchange for the given prefix combinations.

5. (TTS_only) For the TTS predictions, it is hard to find the best start codon matching the predicted stop codon, as multiple start codons can be present in-frame upstream of the predicted stop codon. The method used in the TTS_finder script, finds the shortest possible ORF, by choosing the first in-frame start-codon. In this step, the longest possible ORF is added to the results for a given predicted stop. First, the first in-frame stop codon upstream of the current predicted stop-codon is searched, then the first start-codon downstream of the upstream stop-codon is chosen. This ensures that the detected start-codon is the furthest possible in-frame start-codon that ensures that no additional in-frame stop-codon is between the predicted stop-codon and the attributed start-codon. This provides us with the longest possible ORF.
For TIS predictions this is not necessary, as we start from the predicted start-codon and look for the first in-frame stop-codon.

6. Next, the the excel files are filtered and split into three different files. The idea is to better distinguish un-annotated ORFs that were predicted. To do this, we ensure that no annotated stop is within 25nt of a predicted stop, otherwise it is filtered out. This is done once for the upstream direction, the downstream direction and both, resulting in 3 different files. This is just an additional method that might help easier manual investigation of potentially new ORFs. The original table is also valid, and contains all information.

7. In a final step, expression information is added to all tables for all predicted ORF intervals. This includes both read per kilobase million (RPKM) values and translational efficiency (TE) values. To do this, the read counts are collected using subread-featureCounts. These readcounts are then used in order to calculate both the RPKM and the TE for every sample.

 :warning: **IMPORTANT:** The scripts are written to be compatible with the HRIBO workflow, all samples must be in the form `|method|-|condition|-|replicate|`. `|method|` is either `RIBO`, `RNA`, `TIS`, `RNATIS`. `TTS` and `RNATTS` will be supported soon, until then we suggest labeling `TTS` files `TIS`. `|condition|` can be any string and `|replicate|` can be any integer.

The chosen thresholds, offsets and coverage mappings can change for each organism, therefore it is advised to investigate the data first to ensure that the right parameters are chosen. :warning:


# Running the analysis scripts
The analysis is made up of multiple python3 and bash scripts. If all data is collected as described in the [required files section](#Required-files), running the script will be straight-forward.
Simply run either `tts_finder_analysis_TTS.sh` or `tts_finder_analysis_TIS.sh` depending on the site that is to be analysed.

The following commandline arguments are required:
| Name                | Argument | Description                                                                                                           |
|---------------------|----------|-----------------------------------------------------------------------------------------------------------------------|
| path                | -p       | Path to the project folder, where the wig file folder is located and the result folder will be placed.                |
| scriptpath          | -s       | Path to the TTS_analysis folder, in which all scripts necessary for computation are located.                          |
| annotationpath      | -a       | The annotation file for the organism that is analysed (`.gff3` format)                                                |
| genomepath          | -g       | The genome file for the organism that is analysed (`.fasta` format)                                                   |
| experiments         | -e       | The experiment to be analysed (if more than one, use this option multiple times e.g `-e exp1 -e exp2 ...`)            |
| mappings            | -m       | The mappings to be used for analysis (e.g threeprime, fiveprime, etc...) ( any subfolder under experiment). If more than one, use this option multiple times (e.g. `-m threeprime -m threeprime30 -m fiveprime ...`)       |    
| offsets             | -o       | The p-site offset used for each of the mappings. If you use multiple mappings, ensure that you use the same amount of offsets. (e.g. `-m threeprime -m fiveprime`, `-o 8 -o 10` means that the 3' files have an offset of 8 and the 5' files have an offset of 10).                          |
| normalizations      | -n       | The normalizations to be used (e.g. raw, mil, min) (any subfolder under mapping).                                     |
| contrasts           | -c       | The contrasts used for the experiment. If you want log2FC for certain peak-heights in a table you can use this option to indicate which samples should be compared (e.g. RIBO-A-1_TIS-A-1)                                                                                                       |
| bamfolder           | -b       | The path to the bamfiles, ensure that each bamfile has an according index file. If not, use `samtools index |bamfile|` for all files missing the index. If installed you can also use `parallel`[[2]](#2),  `parallel  samtools index ::: *.bam` to run it on all bam files.              |
| tmpfolder           | -t       | The folder where temporary files will be dumped.                                                                      |
| readcountthreshold  | -r       | The readcount threshold used in the analysis, if a position has less than this amount of reads it is ignored (>0)     |



If you have your own data, you can run the scripts individually, each of them is described in the [scripts section](#Scripts) below.



# Scripts

## Sub-scripts
* **TTS_finder.py:**
* **calculate_expression.py:**
* **call_featurecounts.py:**
* **detect_longest_potential_ORFs.py:**
* **excel_utils.py:**
* **get_readcount_gff.sh:**
* **map_reads_to_annotation.py:**
* **merge_TTS.py:**
* **post_filter_tts.py:**
* **

## References
<a id="1">[1]</a> 
Gelhausen, R. (2020).
HRIBO - High-throughput analysis of bacterial ribosome profiling data 
([BioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.27.046219v1))

<a id="2">[2]</a> 
Tange, O. (2011).
[GNU Parallel](http://www.gnu.org/software/parallel/) - The Command-Line Power Tool
[DOI](http://dx.doi.org/10.5281/zenodo.16303)
