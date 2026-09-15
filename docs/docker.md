# Run ORFBounder with Docker

The Docker image contains ORFBounder and its Python dependencies. Use it to run
an analysis on files on your computer and save results back to your computer.

## Build the image

Install and start Docker, download this repository, and open a terminal in the
repository folder. Build a local image:

```bash
docker build -t orfbounder:2.0.0 .
```

The first build downloads the base image and dependencies, so it needs an
internet connection. Once built, the image can analyze your local inputs offline.
The name `orfbounder:2.0.0` in this guide refers to the image you built locally.

View the available options:

```bash
docker run --rm orfbounder:2.0.0
docker run --rm orfbounder:2.0.0 orfbounder --help
docker run --rm orfbounder:2.0.0 orfbounder-batch --help
```

Use `orfbounder` for one sample or matched set of assays. Use `orfbounder-batch`
for experiments described by a TSV configuration table.

## Run one sample

The following Linux shell example assumes your current directory contains:

```text
inputs/
├── genome.fa
├── annotation.gff
├── offsets.json
└── alignments/
    └── TIS-WT-1.bam
```

Supply your own reference, annotations, alignments and calibrated offsets.
See [preparing inputs](../README.md#prepare-your-own-inputs) and
[setting offsets](../README.md#set-offsets-and-read-lengths). This example uses
`fiveprime` mapping; change it to match your calibration.

```bash
mkdir -p results
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$(pwd)/inputs,target=/data,readonly" \
  --mount "type=bind,source=$(pwd)/results,target=/results" \
  orfbounder:2.0.0 orfbounder \
  --alignment_file_tis /data/alignments/TIS-WT-1.bam \
  --annotation_file /data/annotation.gff \
  --genome_file /data/genome.fa \
  --offset_json /data/offsets.json \
  --mapping_method fiveprime \
  --output_path /results \
  --output_basename WT-1
```

The first mount makes your `inputs` folder available as `/data` inside the
container, with read-only access. The second mount saves `/results` into your
local `results` folder. Paths passed to ORFBounder must use these **container
paths**. Files outside the mounted folders are unavailable to the analysis.
The `--user` option uses your current user and group so the result files belong
to you. Create the results folder before running the command and ensure you
can write to it.

Open **`results/WT-1.report.html`** on your computer after the command finishes.
The container is removed by `--rm`; the files in your local results folder
remain. For another analysis, choose a new results folder or output basename.
ORFBounder protects existing result files from being overwritten.

For TTS input, use `--alignment_file_tts /data/alignments/TTS-WT-1.bam` with the
appropriate mapping method and offsets. To include a matched TIS/TTS pair or
RIBO measurements, add their alignment options to the same command. Other
analysis options work as shown in [the main usage guide](../README.md#run-one-sample).

## Run multiple samples

Put your configuration table at `inputs/experiments.tsv` and your alignments
in `inputs/alignments`. Create the table using
[the TSV template](../templates/template_config_spreadsheet.tsv). Every input
path in it must refer to a mounted location inside the container. For example:

| TSV column | Example value |
| --- | --- |
| `experiment_name` | `WT_tis` |
| `annotation_file_path` | `/data/annotation.gff` |
| `genome_file_path` | `/data/genome.fa` |
| `offset_file_path` | `/data/offsets.json` |
| `TIS_folder_path` | `/data/alignments` |
| `mapping_method` | `fiveprime` |
| `normalization_method` | `raw` |

Save these as columns in a tab-separated file, with one row per experiment.
Use absolute container paths, including for optional read-length, expression
and mapped-count inputs. A host path such as `/home/me/inputs/genome.fa` will
not resolve inside this container; the mounted path is `/data/genome.fa`.

Validate the configuration first:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$(pwd)/inputs,target=/data,readonly" \
  orfbounder:2.0.0 orfbounder-batch \
  -c /data/experiments.tsv --validate-only
```

Then run the analysis:

```bash
mkdir -p batch-results
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$(pwd)/inputs,target=/data,readonly" \
  --mount "type=bind,source=$(pwd)/batch-results,target=/results" \
  orfbounder:2.0.0 orfbounder-batch \
  -c /data/experiments.tsv -r /results
```

For the example settings, open
**`batch-results/WT_tis/fiveprime/raw/report.html`**. Each experiment, mapping
method and normalization method has its own result folder.

## Use the results

The HTML report opens in a web browser and works offline. Search or filter
candidates, expand their evidence, then follow the links to spreadsheets,
nucleotide and protein sequences, GFF annotations and WIG coverage tracks.
Keep the whole result folder together so those links continue to work.

For a single sample, the spreadsheet is `table_per_sample/WT-1.xlsx` inside
the results folder. For the batch example it is
`batch-results/WT_tis/fiveprime/raw/WT_tis_final.xlsx`. A successful single run
also writes `WT-1.complete.json`; a batch result folder contains `complete.json`.
If the command fails, read the terminal error before using any output.

See [reading and using results](results.md) for interpreting ORF boundaries,
comparing sample evidence, choosing candidates and loading genome-browser files.
