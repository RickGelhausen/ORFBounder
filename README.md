# ORFBounder

ORFBounder finds candidate bacterial open reading frames (ORFs) from ribosome
profiling peaks at translation initiation sites (TIS) and termination sites
(TTS). Use its HTML report to explore candidates, compare samples and open
tables, sequences and genome-browser files.

- **TIS data:** a start-codon peak proposes an ORF ending at the next in-frame stop.
- **TTS data:** a stop-codon peak proposes an ORF beginning at an upstream in-frame start.
- **TIS and TTS together:** predictions with identical coordinates share a row.
  Either assay can contribute a candidate; both assays are not required.
- **RIBO and RNA data:** optional measurements describe peak ratios, expression
  and translation efficiency for the candidates.

Use [Docker](docs/docker.md) or install with Python below, then follow
[the results guide](docs/results.md) to inspect an ORF.

## Install

### Docker

With Docker installed and running, open a terminal in this repository and run:

```bash
docker build -t orfbounder:2.0.0 .
docker run --rm orfbounder:2.0.0
```

The second command displays the help. The image includes both `orfbounder` and
`orfbounder-batch`. Follow [the Docker guide](docs/docker.md) to mount your inputs,
run an analysis and open the results.

### Python

Use Python **3.11 or later**. Download this repository, open a terminal in its
folder, and run:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
orfbounder --version
```

Activate the environment with `source .venv/bin/activate` again when opening a
new terminal. Linux is the tested platform. If you use conda, run
`conda env create -f environment.yml` followed by `conda activate orfbounder`
instead. With [uv](https://docs.astral.sh/uv/), use `uv sync --locked`
and prefix subsequent commands with `uv run`.

The commands are `orfbounder` for one sample or matched set of assays and
`orfbounder-batch` for multiple samples.

## Try the supplied example

From the repository folder, with the environment active:

```bash
python -m helpers.toy_example_generation --output-dir demo
orfbounder-batch -c demo/config/config.tsv --validate-only
orfbounder-batch -c demo/config/config.tsv -r demo/results
```

Choose a new folder name if `demo` already contains files. The example generates
a small reference genome and two TIS/RIBO replicates, so no external data are
needed.

Open **`demo/results/toy_example/fiveprime/raw/report.html`** in a web browser.
The report works offline. You should find three annotated ORFs on `NC_000913.3`,
at **100–162, 300–362 and 500–562**, each 63 nucleotides long, with TIS calls in
both replicates. The spreadsheet is `toy_example_final.xlsx` in the same folder.
The example has local statistics switched off; missing q-values are expected.

To practice using the output, search for `100-162`, inspect its sample evidence,
then load the result GFF and WIG tracks alongside `demo/data/genome.fa` in your
genome browser. See [reading and using results](docs/results.md) for the next steps.

## Prepare your own inputs

You need a reference genome in **FASTA**, CDS annotations in **GFF3**, and at
least one **TIS or TTS SAM/BAM** alignment. All must use the same assembly,
contig identifiers and reference lengths. SAM and unindexed BAM are accepted;
your genome browser may require BAM indexes. FASTA sequence identifiers and
alignment-header reference names must be nonempty and unique.

Multipart CDS rows are joined by their GFF3 `ID`. Annotation labels can be on
the CDS itself or inherited through its `Parent` chain (for example,
CDS → mRNA → gene); use normal GFF3 percent encoding for reserved characters in
IDs and Parent values. Each attribute name may occur only once on a feature row;
write multiple parents as one comma-separated `Parent=parent1,parent2` value.

Supply offsets calibrated for your assay, organism, read lengths and selected
read end. Calibration is usually based on metagene profiles: an incorrect offset
can assign a peak to the wrong codon. [The offset instructions below](#set-offsets-and-read-lengths)
explain the file format and sign convention.

For batch runs and matched expression measurements, name files like this:

```text
alignments/
├── TIS-WT-1.bam
├── TTS-WT-1.bam
├── RIBO-WT-1.bam
├── RNA-WT-1.bam
├── TIS-WT-2.bam
└── RIBO-WT-2.bam
```

Names follow `METHOD-condition-replicate.bam` or `.sam`. Conditions contain
letters and digits; replicate numbers are positive integers. A suffix such as
`TIS-WT-1_aligned.bam` is accepted; JSON sample names are still `TIS-WT-1`.
Keep only one alignment per sample in each input folder. Matching uses the exact
condition and replicate, so replicate 1 and replicate 10 are separate samples.
Direct calling files may use custom labels when no folder-based expression
matching is needed. If all calling files in one direct command use the standard
pattern, ORFBounder verifies their assay prefixes and requires the same condition
and replicate before creating output.

### Circular references

Circularity is opt-in and is never guessed from a contig name. For every
circular FASTA sequence, declare a whole-reference `region` in the input GFF3;
including the matching sequence-region directive is recommended:

```gff3
##sequence-region plasmid 1 5000
plasmid	RefSeq	region	1	5000	.	.	.	ID=plasmid-region;Is_circular=true
```

The region must span exactly `1..reference_length`, matching the FASTA length
and, when the sequence occurs in an alignment, its SAM/BAM header length. Input
CDS coordinates remain physical and within that range.
Represent an annotated CDS crossing the origin as multiple physical GFF3 rows
with the same `ID`; do not use `start > end` or a coordinate above the reference
length. Undeclared sequences remain linear, so their searches and windows do
not wrap.

ORF tables linearize a crossing candidate: `Start` remains in `1..L`, while
`Stop` can be greater than `L`. `Reference_length` and `Is_circular` make that
representation explicit. Candidate GFF3 files convert it back to physical
multipart CDS rows, and WIG tracks always use physical coordinates. See
[circular coordinates in the results guide](docs/results.md#circular-reference-coordinates).

## Run one sample

For a TIS sample, using your own file paths and calibrated mapping method:

```bash
orfbounder \
  --alignment_file_tis alignments/TIS-WT-1.bam \
  --annotation_file annotation.gff --genome_file genome.fa \
  --offset_json offsets.json --mapping_method fiveprime \
  --output_path results/tis-WT-1 --output_basename WT-1
```

For a TTS sample:

```bash
orfbounder \
  --alignment_file_tts alignments/TTS-WT-1.bam \
  --annotation_file annotation.gff --genome_file genome.fa \
  --offset_json tts-offsets.json --mapping_method threeprime \
  --tts_start_selection furthest_inframe \
  --output_path results/tts-WT-1 --output_basename WT-1
```

`furthest_inframe` selects the furthest upstream allowed start before the next
in-frame stop. Use `next_inframe` to select the closest upstream start instead.
In a TTS-only result, the start is inferred from sequence. It has no measured
TIS support unless TIS data are also supplied and call that boundary.

For a matched TIS/TTS pair, with optional RIBO measurements:

```bash
orfbounder \
  --alignment_file_tis alignments/TIS-WT-1.bam \
  --alignment_file_tts alignments/TTS-WT-1.bam \
  --alignment_file_ribo alignments/RIBO-WT-1.bam \
  --annotation_file annotation.gff --genome_file genome.fa \
  --offset_json offsets.json --mapping_method fiveprime \
  --statistics local \
  --output_path results/paired-WT-1 --output_basename WT-1
```

Each run uses one mapping method for every supplied assay, with separate offsets
allowed per sample. Choose a method for which all supplied assays are calibrated.
TTS offsets need termination-specific calibration; initiation advice alone does
not provide that calibration.

Open `WT-1.report.html` in the selected output folder after the command finishes.
Use `orfbounder --help` for all options.

## Set offsets and read lengths

An offset JSON maps sample names to read-length offsets. For example, a
fiveprime calibration might produce:

```json
{
  "TIS-WT-1": {"28": -12, "30": -13},
  "RIBO-WT-1": {"28": -12, "30": -13}
}
```

These numbers illustrate the format; use your own calibrated values. An offset
of **-12 with fiveprime** moves the endpoint 12 bases downstream in transcript
direction. An offset of **+12 with threeprime** moves it 12 bases upstream.
This applies on both strands. When using offsets from another tool, check its
sign convention before preparing this file.

A `"default"` length supplies a fallback within a sample's map. A top-level
`"default"` map is used when the sample itself is absent. A sample-specific map
does not inherit missing lengths from the top-level map. Leave out fallbacks
when uncalibrated samples or lengths should cause an error.
Write length keys as ordinary positive decimal integers without leading zeros
(`"30"`, not `"030"`), so they match the read's query length exactly.

Use a separate read-length JSON to select calibrated lengths:

```json
{"TIS-WT-1": "28,30", "RIBO-WT-1": "28,30"}
```

Pass it with `--read_length_json read_lengths.json`. Lengths can be integers or
comma-separated values and ranges, such as `"28-30,32"`; omitting this file uses
all lengths. Expanded selections are limited to 100,000 lengths per sample to
catch accidental enormous ranges. Read length is query length, including
soft-clipped bases.
When the file is supplied, every TIS, TTS and RIBO calling sample must have an
exact key or the JSON must contain a top-level `"default"`. Expression samples
also need a key when `rpkm_read_usage` is `specific`. Validation rejects missing
entries so a misspelled sample name cannot silently disable length filtering.

## Run several samples with a configuration table

Copy [the TSV template](templates/template_config_spreadsheet.tsv), fill in your
paths and save as **tab-separated text**. Each row is one experiment. Optional
columns can be omitted; empty optional cells use defaults.

| Required column | What to enter |
| --- | --- |
| `experiment_name` | Unique name, such as `WT_tis`; start with a letter or digit and use only letters, digits, `_`, `-`, `.` |
| `annotation_file_path` | GFF3 annotation file |
| `genome_file_path` | FASTA genome file |
| `offset_file_path` | Calibrated offset JSON file |
| `TIS_folder_path` or `TTS_folder_path` | At least one alignment folder; supply both for paired assays |

TIS, TTS and RIBO columns can point to the same folder; filenames select the
assay. **Relative input paths start from the directory where you run the
command**, not the TSV's directory. Absolute paths avoid that ambiguity.

```bash
orfbounder-batch -c experiments.tsv --validate-only
orfbounder-batch -c experiments.tsv -r results
```

Validation checks configuration, sample matching, offset coverage by sample and
reference headers before analysis. A missing read-length-specific offset can
still appear only when the analysis encounters a record of that length.

The report for each experiment is at
`results/<experiment>/<mapping>/<normalization>/report.html`. Replicates are
merged by ORF coordinates while preserving their separate measurements.

### Choose analysis settings

| TSV column | Default | How it affects your analysis |
| --- | --- | --- |
| `RIBO_folder_path` | blank | Include matched RIBO peak measurements |
| `mapping_method` | `threeprime` | `fiveprime`, `threeprime`, `centered`, `global`; match your calibration |
| `read_length_json` | blank | File selecting read lengths; blank accepts all lengths |
| `min_peak_height` | `5` | A position must be **strictly greater** than this value to contribute to a called peak |
| `peak_height_operator` | `max` | Highest passing position in a codon window; `sum` adds all passing positions |
| `tts_start_selection` | `furthest_inframe` | Furthest or closest (`next_inframe`) upstream start |
| `genetic_code` | `11` | NCBI translation-table number for the organism |
| `start_codons`, `stop_codons` | See below | Comma-separated DNA triplets for boundary searching |
| `normalization_method` | `raw` | `raw`, `mil` or `min`, explained below |
| `normalization_scope` | `contig` | Scale using reads on the same contig; `library` uses all contigs |
| `mapped_counts_file_path` | blank | Required for `min`; tab-separated sample, contig and count, without a header |
| `alignment_folder_path` | blank | Matched assay/RNA alignments for RPKM and TE |
| `rpkm_read_usage` | `all` | All accepted lengths for expression; `specific` uses selected lengths |
| `gff_output_mode` | `combined` | One GFF per sample; `split` also writes category-group files (ordinary and complex exact annotation matches share the annotated file) |
| `statistics` | `none` | `local` adds local enrichment tests with endpoint mapping |
| `background_width` | `50` | Bases in each statistical background flank |
| `fdr` | `0.05` | Adjusted p-value threshold for local support |
| `fdr_method` | `by` | `by` or `bh` correction; see [statistics](docs/statistics.md) |
| `min_mapq` | `0` | Exclude alignment records below this MAPQ |
| `duplicates` | `include` | Include or exclude records carrying the SAM duplicate flag |
| `multimappers` | `exclude` | Exclude or include records with an `NH` tag greater than 1 |
| `matched_comparison_file_path` | blank | Optional comparison TSV for exploratory fixed-margin, same-assay support; requires `statistics=local` |
| `matched_fdr` | `0.05` | Adjusted p-value threshold used by the exploratory matched-support rule |
| `matched_fdr_method` | `by` | `by` or `bh` correction for each matched comparison |

Comma-separated mapping or normalization methods create separate runs, for
example `raw,mil`. Choose settings before inspecting statistical results;
q-values do not account for selecting a favorable parameter run.

### Genetic code and boundary codons

The default is **NCBI table 11**, with starts `ATG,GTG,TTG` and stops
`TAG,TAA,TGA`. Set `--genetic-code` in a direct run or `genetic_code` in the TSV
for another code. Other codes use their full table-defined start and stop sets
when these lists are omitted. Choose the code appropriate to your organism using
[NCBI's genetic-code tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

Explicit lists restrict boundary searching to a subset of that code's allowed
start/stop codons. Direct commands use spaces, such as
`--start_codons ATG GTG TTG`; TSV cells use commas. Tables with context-dependent
stop/sense codons (27, 28 and 31) are unsupported. The selected code also determines
amino-acid translation; see [using sequences](docs/results.md#use-nucleotide-and-protein-sequences).

### Mapping and normalization

`fiveprime` and `threeprime` put one count at each selected read endpoint after
applying its offset. They support local enrichment testing. `global` counts each
aligned base. `centered` trims 11 query bases from each end and shares one count
among the remaining aligned bases; reads without a remaining aligned base
contribute no coverage.

`raw` keeps counts unchanged. `mil` expresses counts per million accepted
reads. With `contig` scope, each contig uses its own accepted-read total; with
`library`, the denominator is the total across all contigs. Both strands
contribute to these totals. The scope also sets the RPKM denominator, including
when peak counts are `raw`. Use the same scope across samples you compare.

`min` scales to the smallest provided accepted-read count. With `contig` scope,
the target is calculated separately for each contig. With `library` scope, the
rows are first summed by sample across contigs and the smallest sample total is
used as the library-wide target. To prepare the count file:

```bash
python -m helpers.recover_total_counts -c experiments.tsv -r mapped_counts
```

Set `mapped_counts_file_path` to the generated `<experiment>_mapped_reads.tsv`.
Counts use the same alignment policy and selected read lengths as the analysis.
Each input row must be unique and positive, and the file must contain every
calling sample in that experiment; a missing or misspelled calling label is
rejected so it cannot change the minimum silently. Additional samples may be
included intentionally as normalization references. With contig scope, every analyzed
contig needs a count; with library scope, every analyzed sample needs one or
more rows.

For direct runs, use `--normalization_method mil --normalization-scope library`
for whole-library scaling. Normalization changes the meaning of the peak-height
threshold; it does not change the raw counts used for statistics.

### Include expression measurements

Set `alignment_folder_path` in a batch TSV, or pass
`--alignment_file_path alignments` in a direct run, to calculate RPKM. Include
matched assay files and any RNA libraries in that folder. Translation efficiency
(TE) uses the same condition and replicate in these pairs:

| Assay | Matched RNA name |
| --- | --- |
| `RIBO-WT-1` | `RNA-WT-1` |
| `TIS-WT-1` | `RNATIS-WT-1` |
| `TTS-WT-1` | `RNATTS-WT-1` |

Batch expression counting uses all lengths by default. Direct runs use the
selected lengths unless you add `--all_reads_rpkm`. Expression length selection
is separate from peak-calling length selection. See
[expression and ratios](docs/results.md#interpret-expression-and-ratios)
before comparing values.

## Add statistical support

Add `--statistics local` to a direct command or enter `local` in the TSV's
`statistics` column. Results then include raw endpoint/background counts,
width-normalized local fold enrichment, p-values, adjusted p-values (q-values)
and local support flags for each assay.
The report and table still include calls without statistical support.

Use the effect size and q-values to prioritize inspection alongside called
peaks, annotation, coverage and replicate consistency. The test describes local endpoint
enrichment; it does not compare treatments with controls or establish that an
ORF is translated. Read [the statistical guide](docs/statistics.md) and
[the candidate inspection workflow](docs/results.md#choose-candidates-to-inspect).

For a batch, ORFBounder can also calculate **exploratory fixed-margin matched
support** between two conditions within one TIS or TTS assay. Copy
[the comparison template](templates/template_matched_comparisons.tsv), or create
a tab-separated file such as:

```text
comparison	assay	numerator_condition	denominator_condition
treated_vs_control	TIS	treated	control
```

The corresponding files must be named, for example, `TIS-treated-1.bam` and
`TIS-control-1.bam`, with at least two identical replicate IDs in both
conditions. Each ID must represent a real, predeclared experimental pair or
block, and those blocks must be independent; arbitrary ordinal matching and
technical lane splits are not biological replication. Set `statistics=local`
and put this file in `matched_comparison_file_path`.

The directional fixed-margin test asks whether the numerator has greater
peak-versus-local-background endpoint odds while retaining each block's local
read depth. A candidate needs at least two informative pairs; support requires
its q-value to meet the configured threshold, a common odds ratio above 1 and a
strict majority of informative pairs favoring the numerator. Local and matched
effect directions at the neutral boundary are evaluated from exact count
relationships rather than rounded display ratios. This is read-depth evidence
under a conditional endpoint-count model: it does not estimate
biological dispersion, test a population-level condition effect, compare TIS
against RIBO, test overall abundance or combine TIS and TTS evidence.
The aggregate candidate table and a companion `_strata.tsv` table retain the
summary and per-pair counts, respectively. Join them by `candidate_id`. See
[matched-condition inference](docs/statistics.md#matched-condition-inference).

## Use the installed utilities

The Python package and Docker image include the following utility modules. Run
them with the same Python environment as ORFBounder:

```bash
python -m helpers.toy_example_generation --output-dir demo
python -m helpers.recover_total_counts -c experiments.tsv -r mapped_counts
python -m lib.merging -t first.csv second.xlsx -o merged.xlsx
python -m helpers.final_to_gff -i merged.xlsx -o colored.gff
```

`lib.merging` accepts ORFBounder Excel files and tab-separated `.csv`/`.tsv`
tables, then writes a merged `.xlsx`, tab-separated `.csv` and `.gff`. Sample
identities come from table column names; repeated measurements must be missing
or identical, otherwise the merge stops instead of choosing one by input order.
Local and matched statistical evidence is retained. Integer audit fields are
parsed exactly; malformed count representations are rejected, and spreadsheet
counts above `2^53` are stored as decimal text to avoid silent rounding.
Nullable support flags accept only Boolean/blank values (including `1`/`0`),
so TSV and XLSX inputs cannot silently disagree about statistical support.
`helpers.merge_tables` remains an installed compatibility entry point, but
`python -m lib.merging` is the supported interface. `helpers.final_to_gff`
creates a GFF3 from an existing result table and colors features by their
TIS/TTS call evidence. Its importer preserves literal identity labels such as
`NA` and leading-zero reference names, rejects duplicate headers, and checks a
coordinate-form identifier against its coordinate columns. These utilities
publish complete files without replacing an existing target. Choose new output
paths with no symbolic-link components; derived files are not written into a
completed result folder in place.

## Troubleshooting

| Problem | What to do |
| --- | --- |
| Command not found | Activate your environment; use `orfbounder-batch` for TSV runs |
| Missing/unknown configuration columns | Save tab-separated text using the current template; the read-length column is `read_length_json` |
| No samples found | Check folders and the exact `METHOD-condition-replicate` filename pattern |
| Duplicate sample | Keep one SAM or BAM per sample; move alternative alignments out of the input folder |
| Reference mismatch | Use the same assembly and contig identifiers in FASTA, GFF and alignments |
| Missing offset for a sample or length | Add that sample to the offset JSON (or an intentional top-level default); for a missing length, select calibrated lengths or add it to the sample's offset map |
| Zero ORFs | Check accepted counts in run JSON, offset signs, selected lengths, genetic code and threshold; height 5 does not pass threshold 5 |
| No q-values | Enable `statistics=local` with `fiveprime` or `threeprime`; see [missing values](docs/results.md#understand-missing-values) |
| Missing RPKM or TE | Supply the expression folder and exactly matched assay/RNA names |
| Existing results error | Select a new output directory to preserve the earlier analysis |
| Interrupted or failed analysis | Read the terminal error, correct the input and rerun; use only outputs marked complete |

ORFBounder retains primary mapped alignments passing the SAM QC flag. Secondary,
supplementary, unmapped and QC-failed records are always excluded. Configure
MAPQ, duplicate-flag and `NH > 1` filtering with `min_mapq`, `duplicates` and
`multimappers`; without an `NH` tag, a multimapper cannot be identified by that
rule. A present `NH` tag must be a positive integer; malformed values stop the
run rather than being treated as unique alignments. Coverage respects CIGAR gaps,
and every primary mapped alignment must provide a CIGAR. When a valid SAM record
omits `SEQ`, CIGAR supplies the query length needed for read-length filters and
offset selection. Expression counts an alignment once when any of its blocks
overlaps the ORF on the same strand. Paired mates are separate alignment records;
ORFBounder does not infer fragments. The run JSON records the policy and counts
for every exclusion reason.

Interpret candidates alongside your controls and biological replicates.

## References and license

The TIS peak-calling approach was inspired by
[RETscript](https://doi.org/10.1093/nar/gkaa304).
License: [GPL-3.0](LICENSE).
