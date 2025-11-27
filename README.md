# ORFBounder

Detection of potential start/stop codons based on Translation Initiation Site (TIS) or Translation Termination Site (TTS) peaks, using a similar concept to the [RETscript for TIS](https://www.sciencedirect.com/science/article/pii/S1097276519301078).

These scripts were created to be used with the metagene-profiling and alignment files created by the [HRIBO workflow](https://github.com/RickGelhausen/HRIBO) [[1]](#references). Nevertheless, `ORFBounder` can also be used with other standard alignment files and metagene-profiling tools.

---

## Table of Contents

- [What is ORFBounder?](#what-is-orfbounder)
- [Requirements](#requirements)
  - [Required Packages](#required-packages)
  - [Required Files](#required-files)
- [Analysis Overview](#analysis-overview)
- [Running ORFBounder](#running-orfbounder)
  - [Step 1: Offset JSON](#step-1-offset-json)
  - [Step 2: Config Spreadsheet](#step-2-config-spreadsheet)
  - [Step 3: Running call_ORFBounder.py](#step-3-running-call_orfbounderpy)
- [Output Files](#output-files)
  - [Results Directory Structure](#results-directory-structure)
  - [Output Tables](#output-tables)
- [Toy Example](#toy-example)
- [References](#references)
- [Citation](#citation)
- [License](#license)

---

## What is ORFBounder?

ORFBounder is a tool that aids in the detection of potential Open Reading Frames using TIS, TTS and Ribo-seq data, based on previously determined offsets.

**ORFBounder can be run in three modes:**

1. **TIS-only mode**: Detects start codons with strong TIS peaks and creates ORFs using the next in-frame stop codon
2. **TTS-only mode**: Detects stop codons with strong TTS peaks and creates ORFs using either the closest or furthest in-frame start codon (that does not overlap with another in-frame stop codon)
3. **Combined TIS/TTS mode**: Runs both TIS and TTS analyses individually, then combines the results in an experimental approach that pairs start codons detected by TIS with stop codons detected by TTS to form high-confidence ORFs

Ribo-seq data is used as a control and to calculate fold-changes in order to determine ORFs with high confidence.

---

## Requirements

### Required Packages

The scripts used in the analysis exclusively require Python 3 (>3.6).

| Package | Minimum Version |
|------------|----------------|
| pysam | 0.23.3 |
| pandas | 2.3.3 |
| numpy | 2.3.4 |
| biopython | 1.8.6 |
| interlap | 0.2.7 |
| openpyxl | 3.1.5 |
| xlsxwriter | 3.2.9 |
| xlrd | 2.0.2 |
| pytest | 8.4.2 |
| pytest-mock | 3.15.1 |
| pyyaml | 6.0.3 |

**Installation:**

We recommend using [uv](https://docs.astral.sh/uv/) for dependency management. The required packages will be automatically downloaded when you run ORFBounder:

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Run ORFBounder (dependencies will be automatically managed)
uv run call_ORFBounder.py -c config.tsv -r output
```

> **⚠️ Platform Note:** This tool was developed and tested on a Linux system. It should not contain Linux-specific commands, but it was never tested on Windows or macOS.

### Required Files

| File | Description |
|-------------------|------------------------------------------------|
| `genome.fa` | A genome file in FASTA format for the analyzed organism |
| `annotation.gff` | An annotation file in GFF3 format for the analyzed organism (tested using annotation files from NCBI) |
| `alignment.(sam\|bam)` | Alignment files in `.sam` or `.bam` format |

> **⚠️ IMPORTANT - File Naming Convention:**
>
> The scripts are written to be compatible with the [HRIBO workflow](https://github.com/RickGelhausen/HRIBO). All samples (`.sam|.bam`) must be in the form:
>
> **`<method>-<condition>-<replicate>.(sam|bam)`**
>
> - `<method>` is either `RIBO`, `RNA`, `TIS`, `RNATIS`, `TTS`, or `RNATTS`
> - `<condition>` can be any string (avoid using special characters, e.g., A, B, C, pH4, xyz123)
> - `<replicate>` can be any integer (e.g., 1, 2, 3, 4)

---

## Analysis Overview

For details on the analysis itself, please refer to our publication *(in progress)*.

Analyzing data using ORFBounder is done in three steps:

1. **Determine offsets**: Using metagene profiling output, determine the best mapping and offset combinations. This step is currently manual and requires the user to determine the best mapping/offset and read lengths to be used by `ORFBounder`.
We suggest using the metagene workflow described in `

2. **Prepare configuration files**: Create an offset JSON file and a config spreadsheet with experimental parameters.

3. **Run ORFBounder**: Execute `call_ORFBounder.py`, which allows you to run single or multiple experiments with different parameterizations simultaneously.

---

## Running ORFBounder

We recommend using the `call_ORFBounder.py` script, as it allows running ORFBounder on multiple experiments and with multiple parameterizations at the same time.

It requires two input files:
- An **offset JSON** file
- A **config spreadsheet**

Using these files reduces the number of input parameters and makes it easier to reproduce results later.

> 🚧 Templates for both files can be found in the [templates](templates/) folder.

### Step 1: Offset JSON

ORFBounder requires a set of offsets for each input file (TIS, TTS, or RIBO). These offsets can be determined using any metagene-profiling tool.

ORFBounder allows different offsets for different samples and per read length. To enable this, offsets should be provided in JSON format.

**Example offset JSON:**

```json
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

This file allows setting specific offsets for each file and each read length within the file.

> **⚠️ Important Notes:**
>
> - In the offset file, use the prefixes of the alignment files in format: `<method>-<condition>-<replicate>`
> - If a value is not specified in the offset JSON file, the default value will be used. **Make sure that default values are set.**
> - The minimal JSON file would be:
>   ```json
>   {
>     "default": {
>       "default": 10
>     }
>   }
>   ```
>   In this case, the offset for all files and all read lengths is set to 10.

> **💡 Tips:**
>
> - Negative offsets are supported and are common for certain organisms and mapping methods
> - Multiple offsets for the same read length (within the same file) are currently not supported

---

### Step 2: Config Spreadsheet

The config spreadsheet is a tab-separated table that contains all input parameters for ORFBounder.

- **Columns** describe different parameters (some required, some optional)
- **Rows** describe different independent experiments

#### Required Parameters

| Column | Description |
|----------------------|------------------------------------------------|
| `experiment_name` | The name given to the experiment (e.g., exp1, run001). We suggest not using special characters. |
| `annotation_file_path` | The path to the GFF3 format annotation file |
| `genome_file_path` | The path to the FASTA format genome file |
| `TIS_folder_path` | The path to the folder containing the TIS alignment files in `.sam\|.bam` format. Can be the same as TTS or RIBO. |
| `TTS_folder_path` | The path to the folder containing the TTS alignment files in `.sam\|.bam` format. Can be the same as TIS or RIBO. |
| `RIBO_folder_path` | The path to the folder containing the RIBO alignment files in `.sam\|.bam` format. Can be the same as TIS or TTS. |
| `normalization_method` | The normalization method(s) to be used. Multiple methods can be given and result in multiple ORFBounder runs (e.g., `raw,mil` or `min` or `mil,min,raw`). See [Normalization Methods](#normalization-methods) for details. |
| `mapping_method` | The mapping method to be used. Determines which part of reads are used in the analysis. For TIS and TTS analysis, `fiveprime` or `threeprime` usually perform better because they result in sharper peaks (e.g., `fiveprime`, `threeprime`, `centered`, or `global`). See [Mapping Methods](#mapping-methods) for details. |
| `offset_file_path` | The path to the custom offset JSON file. See [Offset JSON](#step-1-offset-json) section for details. |

> **💡 Tips:**
>
> - `TIS_folder_path`, `TTS_folder_path`, and `RIBO_folder_path` (even `alignment_folder_path`) can point to the same directory. Due to the `<method>-<condition>-<replicate>` naming scheme, ORFBounder will collect and match the correct samples and accumulate all results in one file. These are only split into different parameters should the different file types be stored in different folders.
> - You do not have to run ORFBounder separately for each replicate. All replicates present in the respective folder paths will be run automatically.

#### Optional Parameters

| Column | Default Value | Options | Description |
|-----------------------|------------------|---------|--------------------------------|
| `read_lengths` | -1 (all lengths) | `24,25,26,28` or `24-26,28` | Read lengths used in the analysis. These can be given as intervals and/or single values. |
| `start_codons` | `ATG,TTG,GTG` | Any list of codons | Start codons used in the analysis |
| `stop_codons` | `TAG,TAA,TGA` | Any list of codons | Stop codons used in the analysis |
| `alignment_folder_path` | No TE/RPKM | Path to alignment files or nothing | Path to a folder containing RIBO and/or RNA alignment files. These will be used to calculate RPKM and TE values for every detected ORF. If no path is given, this analysis step is skipped. |
| `min_peak_height` | 5 | Any positive integer | The minimum height of a TIS/TTS or RIBO peak required to be considered a valid start/stop codon |
| `peak_height_operator` | `max` | `max` \| `sum` | Peaks are searched in a 5nt interval around the start/stop codon. This specifies whether the maximum value or the sum of values within the interval are used. |
| `tts_start_selection` | `furthest_inframe` | `furthest_inframe` \| `next_inframe` | When constructing an ORF based on a stop codon detected using TTS data, either use the next or the furthest in-frame start codon, without overlapping with another stop codon. |
| `log_fold_contrasts` | No Fold Change | List of contrasts (e.g., `RIBO-A-1_TIS-A-1, RIBO-A-2_TIS-A-2`) | Fold changes for any samples present in the experiment can be requested using this parameter. |
| `max_ORF_length` | 150 | Any positive integer | When both TIS and TTS data is used, the data is combined to form ORFs up to a length of `max_ORF_length`. |
| `rpkm_read_usage` | `all` | `all` \| `specific` | When calculating RPKM and TE values, either all reads can be used or the specific read lengths given by the `read_lengths` parameter. |
| `gff_output_mode` | `combined` | `combined` \| `split` | Output a single combined GFF file or a GFF file for each ORF type. |

> **⚠️ Important Notes:**
>
> - Even though these parameters are optional, the column headers are not. ORFBounder will tell you which headers are missing.
> - The values are optional and default values will be used if not specified.
> - Ensure that the file is tab-separated. You can modify the template file with any spreadsheet viewer.

---

### Step 3: Running call_ORFBounder.py

After creating both input files, running ORFBounder is straightforward:

```bash
uv run call_ORFBounder.py -c  -r
```

**Parameters:**
- `-c <config_spreadsheet>`: Path to the configuration spreadsheet
- `-r <output_folder>`: Path to the output directory where results will be saved

---

## Output Files

### Results Directory Structure

For each experiment, ORFBounder creates a structured output directory containing the following:

```
<output_folder>/
└── <experiment_name>/
    └── <mapping_method>_
		├── <normalization_method>/
		│   ├── <experiment>_final.xlsx
		│   └── <experiment>_final.gff
		├── coverage_files/
		│   └── *.wig (or similar coverage files)
		├── table_per_condition/
		│   └── <condition>-<replicate>.xlsx
		└── gff_per_condition/
			├── <condition>-<replicate>_<orf_type>.gff (if gff_output_mode=split)
		  └── <condition>-<replicate>.gff (if gff_output=combined)
```

#### Directory Contents

**1. Mapping/Normalization Folders** (`<mapping_method>_<normalization_method>/`)

For each combination of mapping method (e.g., `threeprime`, `fiveprime`) and normalization method (e.g., `raw`, `mil`, `min`), a separate folder is created containing:

- **`<experiment>_final.xlsx`**: The primary output file containing all detected ORFs with comprehensive annotations and statistics
- **`<experiment>_final.gff`**: Single file with all ORF types.

**2. Coverage Folder** (`coverage_files/`)

Contains coverage files adjusted to the experimental settings (offsets, read lengths, mapping method).

**3. Result Tables Folder** (`table_per_sample/`)

Contains tab-separated tables with detected ORFs for each sample:
- **`<condition>-<replicate>.xlsx`**: Separate tables for each individual condition

**4. GFF Files Folder** (`gff_per_sample/`)

Contains GFF3 format files for genome browser visualization per sample:
- **`<condition>-<replicate>_<orf_type>.gff`**: Separate files for each ORF type (if `gff_output_mode=split`)

---

### Output Tables

The main output file (`final_results.xlsx`) contains the following columns:

#### Static Columns (Always Present)

| Column | Description |
|--------|-------------|
| `Type` | The type/category of the detected ORF |
| `Identifier` | Unique identifier for the ORF |
| `Genome` | Reference genome/chromosome name |
| `Start` | Start position of the ORF (genomic coordinate) |
| `Stop` | Stop position of the ORF (genomic coordinate) |
| `Strand` | Strand orientation (+ or -) |
| `Locus_tag` | Associated locus tag from annotation (if applicable) |
| `Codon_count` | Number of codons in the ORF |
| `Evidence` | Type of evidence supporting the ORF (TIS, TTS, or combined) |
| `Start_codon` | The start codon sequence |
| `Stop_codon` | The stop codon sequence |
| `15nt_window` | 15 nucleotide window around the start codon |
| `Nucleotide_Seq` | Complete nucleotide sequence of the ORF |
| `Amino_Acid_Seq` | Translated amino acid sequence |
| `5'-distance` | Distance to the nearest upstream feature |
| `3'-distance` | Distance to the nearest downstream feature |

#### Dynamic Columns (Depend on Experimental Setup)

The following column types are generated for each sample in your experiment (e.g., `RIBO-WT-1`, `RIBO-WT-2`, `TIS-WT-1`, `TIS-WT-2`):

| Column Pattern | Description |
|----------------|-------------|
| `<sample>_peak_height` | Peak height at the start/stop codon for each sample |
| `<sample1>_<sample2>_log2FC` | Log2 fold change between specified sample pairs |
| `<sample>_relative_density` | Relative read density across the ORF |
| `<sample>_rpkm` | Reads Per Kilobase per Million mapped reads |
| `<sample>_TE` | Translation Efficiency (if RNA-seq data is provided) |

---

## Toy example

We provide a toy example to show how to run the workflow.

The toy example can be generated using:

```
uv run helpers/toy_example_generation.py
```

This will create toy versions of each input file required for `ORFBounder`.

You can run the toy example using:

```
uv run call_orfbounder.py -c toy_example/config/config.tsv -r toy_example/results
```

---

## References

<a id="1">[1]</a> Gelhausen, R. (2020). HRIBO - High-throughput analysis of bacterial ribosome profiling data. [BioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.27.046219v1)

---

## Citation

🚧 *Publication in progress*

## License

GPL-3.0 license


