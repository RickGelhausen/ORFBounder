# ORFBounder

Detection of potential start/stop codons based on Translation Initiation Site (TIS) or Translatation Termination Site (TTS) peaks, using a similar concept than the [RETscript for TIS](https://www.sciencedirect.com/science/article/pii/S1097276519301078).

These scripts were created to be used with the metagene-profiling and alignment files created by the [HRIBO workflow](https://github.com/RickGelhausen/HRIBO) [[1]](#1).
Nevertheless, `ORFBounder` can also be used with other standard alignment files and metagene-profiling tools.

- [ORFBounder](#orfbounder)
  - [What is ORFBounder?](#what-is-orfbounder)
- [Requirements](#requirements)
  - [Required packages](#required-packages)
  - [Required files](#required-files)
- [Analysis](#analysis)
- [Running ORFBounder](#running-orfbounder)
    - [Offset JSON](#offset-json)
    - [Config Spreadsheet](#config-spreadsheet)
      - [Required parameters](#required-parameters)
      - [Optional parameters](#optional-parameters)
  - [Running call_ORFBounder.py](#running-call_orfbounderpy)
- [Output files](#output-files)
- [Additional information](#additional-information)
  - [Normalization methods](#normalization-methods)
  - [Mapping methods](#mapping-methods)
- [References](#references)


## What is ORFBounder?

ORFBounder is a tool that aids in the detection of potential OpenReadingFrames using TIS, TTS and RIBOseq data, based on previously determined offsets.

ORFBounder can be run using only TIS or only TTS data. In these cases, the start codons with a strong TIS peak or the stop codons with a strong TTS peak will be detected and an ORF will be created using the next inframe stop for TIS. For TTS either the closest inframe start codon or the furthest inframe start codon (that does not overlap with another inframe stop codon) will be used to create an ORF.

ORFBounder can also be used with a combination of TIS and TTS data. This will run the original analysis for both TIS and TTS individually. Additionally an experimental approach that combines the start codons detected by TIS and the stop codons detected by TTS to form high confidence ORFs.

Ribo-seq data is used as a control and to calculate fold-changes in order to determine ORFs with high confidence.

# Requirements
## Required packages
The scripts used in the analysis exclusively require python3 (>3.6).


| Tool       | Tested version |
|------------|----------------|
| biopython  | 1.70           |
| interlap   | 0.2.7          |
| numpy      | 1.12.1         |
| pandas     | 0.24.2         |
| pysam      | 0.16.0.1       |
| xlrd       | 2.0.1          |
| xlsxwriter | 3.0.2          |

All required packages are easily retrievable via conda or by using the provided `environment.yml`:

```
conda env create -f environment.yml
conda activate orfbounder_env
```

---
:warning: This tool was developed and tested on a linux system. It should not contain linux specific commands, but it was never tested on Windows or iOS.

---

## Required files

| File              | Description                                    |
|-------------------|------------------------------------------------|
| `genome.fa`     | a genome file in fasta format for the analysed organism. |
| `annotation.gff` | an annotation file in .gff3 format for the analysed organism. (Tested using annotation files from NCBI)
| `alignment.(sam\|bam)` | alignment files in `.sam` or `.bam` format.

---
:warning: **IMPORTANT:** The scripts are written to be compatible with the [HRIBO workflow](https://github.com/RickGelhausen/HRIBO), all samples (`.sam|.bam`) must be in the form `|method|-|condition|-|replicate|.(sam|bam)`.
* `|method|` is either `RIBO`, `RNA`, `TIS`, `RNATIS`, `TTS`, `RNATTS`.
* `|condition|` can be any string, avoid using special characters. (e.g. A, B, C, pH4, xyz123, ...)
* `|replicate|` can be any integer. (e.g. 1, 2, 3, 4, ...)

---

# Analysis

For details on the analysis itself please refer to our publication (in progress). :construction:

Analyzing data using ORFBounder is done in two steps:

1. First, using the metagene profiling output, the best mapping and offset combinations are determined before running the scripts. This step is currently still manual work and requires the user to determine the best mapping/offset and eventually read lengths to be used by `ORFBounder`.

2. Running `call_ORFBounder.py`, which allows you to run single experiments or multiple experiments at the same time.

# Running ORFBounder

We recommend using the `call_ORFBounder.py` script. It allows running ORFBounder on multiple experiments and with multiple parameterizations at the same time. It requires both an `offset JSON` and an `config spreadsheet`. Using these files the amount of input parameters used is reduced and they make it easier to reproduce the results later.

:construction:Templates for both files can be found in the [templates]() folder.:construction:

### Offset JSON

In order to work, ORFBounder requires a set of offsets for each input file (TIS, TTS or RIBO). These offsets can be determined using any metagene-profiling tool.

ORFBounder allows different offsets for different samples and per read-length. To allow this, the offsets should be given in JSON format.

:construction:
A template file with the following example can be found in the [templates folder].
:construction:


```
{
	"TIS-A-1": {
		"33": 12,
		"30": 6,
		"default": 10
  	},
	"RIBO-A-1": {
		"28": 5,
		"30": 6,
		"default": 10
  	},
	"TIS-B-4": {
		"24": 6,
		"30": 6,
		"40": -10,
		"default": 10
  	},
	"default": {
		"31": 10,
		"32": 10,
		"default": 10
  	}
}
```
This file allows us to set specific offsets for each file and each read length within the file.

---
:warning: Please note that in the offset file the prefixes of the alignment files shold be used in format: |method|-|condition|-|replicate|

:warning: If a value is not specified in the offset JSON file. The default value will be used. Make sure that default values are set.
Thus the minimal JSON file would be:

```
{
	"default": {
		"default": 10
  	}
}
```
In this case, the offset for all files and all read lengths is set to 10.

:bulb: Negative offsets are supported and are common for certain organisms and mapping methods.

:no_entry: Multiple offsets for the same read-length (within the same file) are currently not supported.

---

### Config Spreadsheet

The config spreadsheet is a tab-seperated table that contains all the input parameters for ORFBounder.

- Columns in the config file describe the different parameters. Some are required and some are optional.
- Rows in the config file describe different independant experiments.


#### Required parameters


| Column               | Description                                    |
|----------------------|------------------------------------------------|
| experiment_name      | The name given to the experiment (e.g exp1, run001, ...). We suggest not using special characters. There should not be an issue with most of them though. |
| annotation_file_path | The path to the gff3 format annotation file.   |
| genome_file_path     | The path to the fasta format genome file.      |
| TIS_folder_path      | The path to the folder containing the TIS alignment files in .sam\|.bam format. Can be the same as TTS or RIBO. |
| TTS_folder_path      | The path to the folder containing the TIS alignment files in .sam\|.bam format. Can be the same as TIS or RIBO. |
| RIBO_folder_path     | The path to the folder containing the TIS alignment files in .sam\|.bam format. Can be the same as TIS or TTS. |
| normalization_method | The normalization method to be used, multiple methods can be given and result in multiple ORFBounder runs. (e.g. raw,mil or min or mil,min,raw ...). Detailed information can be found in the [Normalization]() section. :construction:|
| mapping_method       | The mapping method to be used. This determines which part of reads are used in the analysis. For TIS and TTS analysis five or threeprime usually perform better, because they result in sharper peaks. (e.g. fiveprime or threeprime or centered or global). Detailed information can be found in the [Mapping]() section. :construction: |
| offset_file_path     | The path to the custom offset file. Detailed explanations can be found in the [Offset JSON](#offset-json) section. |

---
:bulb: `TIS_folder_path`, `TTS_folder_path` and `RIBO_folder_path` (even `alignment_folder_path`) can point to the same directory. Due to the |method|-|condition|-|replicate| naming scheme, ORFBounder will collect and match the correct samples and accumulate all results in one file.
This is only split into different parameters, should the different file types be stored in different folders.

:bulb: You do not have to run ORFBounder seperately for each replicate. All replicates present in the respective folder paths (`TIS_folder_path`, `TTS_folder_path` and `RIBO_folder_path`) will be run automatically.

---

#### Optional parameters


| Column                | Default value    | Options | Description                    |
|-----------------------|------------------|---------|--------------------------------|
| read_lengths          | -1 (all lengths) | 24,25,26,28<br> or <br> 24-26,28        | Read lengths that are used in the analysis. These can be given as intervals and/or single values. |
| start_codons          | ATG,TTG,GTG      | any list of codons | Start codons used in the analysis. |
| stop_codons           | TAG,TAA,TGA      | any list of codons | Stop codons used in the analysis. |
| alignment_folder_path | No TE/RPKM       | path to alignment files <br> or <br> nothing | Path to a folder containing RIBO and/or RNA alignment files. These will be used to calculate RPKM and TE values for every detected ORF. If no path is given this analysis step is skipped. |
| min_peak_height       | 5                | any positive integer | The minimum height of a TIS/TTS or RIBO peak required to be considered a valid start/stop codon. |
| peak_height_operator  | max              | max \| sum        | The peaks are searched in a 5nt interval around the start/stop codon. This specifies whether the maximum value or the sum of values within the interval are used. |
| tts_start_selection   | furthest_inframe | furthest_inframe <br> or <br> next_inframe | When constructing an ORF based on a stop codon detected using TTS data, either use the next or the furthest inframe start codon, without overlapping with another stop codon. |
| log_fold_contrasts    | No Fold Change   | List of contrasts <br> RIBO-A-1_TIS-A-1, RIBO-A-2_TIS-A-2, ...)  | Fold changes for any samples present in the experiment can be requested using this parameter. |
| max_ORF_length        | 150              | any positive integer         | When both TIS and TTS data is used, the data is combined to form ORFs up to a length of max_ORF_length. |
| rpkm_read_usage       | all              | all \| specific     | When calculating RPKM and TE values, eiter all reads can be used or the specifc read lengths given by the read-lengths parameter. |
| gff_output_mode       | combined         | combined \| split      | Output a single combined gff file or a gff file for each orf type. |

---
:warning: Even though these parameters are optional, the column headers are not. ORFBounder will tell you which headers are missing. The values are optional and default values will be used.

:warning: Ensure that the file is tab-seperated, you can modify the template file with any spreadsheet viewer.

---

## Running call_ORFBounder.py

After creating both input files, running ORFBounder is easy.

Simply run the `call_ORFBounder.py` by specifing the config spreadsheet to be used and the path to the output folder.

```
python3 call_ORFBounder.py -c <path/to/config_spreadsheet> -o <path/to/output/folder>
```


# Output files



# Additional information
## Normalization methods
## Mapping methods

# References
<a id="1">[1]</a>
Gelhausen, R. (2020).
HRIBO - High-throughput analysis of bacterial ribosome profiling data
([BioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.27.046219v1))
