# TTS_analysis

Detection of potential start/stop codons based on Translation Initiation Site (TIS) or Translatation Termination Site (TTS) peaks, using a similar concept than the [RETscript for TIS](https://www.sciencedirect.com/science/article/pii/S1097276519301078).

These scripts were created to be used with the metagene-profiling and coverage (.wig) files created by the [HRIBO workflow](https://github.com/RickGelhausen/HRIBO) [[1]](#1).
Nevertheless, the TTS/TIS_finder scripts can also be used with any other coverage files and metagene-profiling tools.

# Requirements

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

# Running the analysis scripts

# Analysis 
To run the scripts, coverage files in .wig format are required. We generated the coverage files using HRIBO [[1]](#1). HRIBO generates coverage files with centered, threeprime, fiveprime and global mappings. Additionally, it provides raw and normalized coverage files for each of the methods (raw, min, mil). 
Then metagene-profiling is performed on all coverage files. 
From the resulting metagene-profiling results, we determined the best p-site-offsets for interesting mapping / normalization combinations. These were then analysed with the TTS_finder scripts.
HRIBO additionally provides .bam files containing all read lengths. If you want to investigate specific read lengths, an example script for filtering .bam files for different read-lengths is provided in this repository.

The analysis is done in multiple steps:

1. First, all potential stop (or start) codons are collected for the TTS (TIS) analysis. The stop codon position is then expanded into an interval of 7nt. These intervals are collected and saved in an InterLap object. The p/a-site offsets are added/substracted from these intervals, in order to ensure that the correct regions are investigated. (The codon intervals are also written to a .gff file. These can be loaded and investigated in a genome-browser.)

2. Next, the requested .wig files are read and if a position passes a given read_count_threshold (default 5), all codon intervals overlapping with the given position are retrieved and their peak_height is incremented by the read count of the position detected. This matches the coverage peaks with given stop(start) codons. 

3. Then, we iterate over all potential codons that have a peak attributed to them. For each stop(codon) the next in-frame start(stop) codon is searched and formed into an ORF prediction. These ORFs are collected and written into a .gff and a .xlsx (excel table) file. The excel file contains a lot of additional information for each predicted ORF (e.g. gene_type, start, stop, strand, locus_tag, codon_count, peak_height, 15nt upstream of the start, nucleotide sequence, amino acid sequence, etc...). Additionally, .gff files for each gene_type are generated for easier investigation in a genome_browser.

4. If the script was used on different RIBO-seq, RNA-seq and TIS or TTS samples, all result tables are bundled into one big excel file, by combining ORF predicted for multiple samples into one row, providing the peak_height information for all involved samples. Contrasts can be given in form of a list of file prefixes (e.g RIBO-A-1_TIS_A-1, RIBO-A-2_TIS_A-1). This will add additional columns with log2foldchange for the given prefix combinations.



**IMPORTANT:** The scripts are written to be compatible with the HRIBO workflow, all samples must be in the form <method>-<condition>-<replicate>. <method> is either RIBO, RNA, TIS, RNATIS. TTS and RNATTS will be supported soon, until then we suggest labeling TTS files TIS. <condition> can be any string and <replicate> any number.

# Scripts

## Main script
* **tts_finder_analysis_TIS.sh:**
* **tts_finder_analysis_TTS.sh:**

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

