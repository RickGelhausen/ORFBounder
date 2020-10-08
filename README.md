# TTS_analysis

Detection of potential start/stop codons based on TIS or TTS peaks, using a similar concept than the [RETscript for TIS](https://www.sciencedirect.com/science/article/pii/S1097276519301078).


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

# Scripts

## Main script
* **tts_finder_analysis.sh:**

## Sub-scripts
* **TTS_finder.py:** 
