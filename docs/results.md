# Read and use your results

Start with the HTML report, then use the spreadsheet and genome-browser files
to inspect the candidates relevant to your question. ORFBounder proposes ORF
boundaries from peak evidence and reference sequence. A candidate's presence
in the results, annotation category and statistical support describe different
aspects of that proposal.

## Open the report and check the run

For a batch experiment, open:

```text
results/<experiment>/<mapping>/<normalization>/report.html
```

For a direct command with `--output_basename WT-1`, open
`<output_path>/WT-1.report.html`. Both reports work in a web browser without an
internet connection. Keep the report in its output folder so links to tables,
FASTA, GFF and coverage files continue to work.

The report contains every candidate and exported field, so its file size grows
with the result table and sequence lengths. Report creation writes that payload
incrementally, and opening the page does not create a second full-table search
index. If a very large report still exceeds your browser's practical limits,
use the linked TSV/Excel table for bulk filtering and return to the report for
shortlisted candidates.

Check the run settings before comparing candidates: reference genome, mapping
end, offsets, selected read lengths, normalization, genetic code and whether
local statistics were enabled. A completed batch result has `complete.json`;
a direct result has `<basename>.complete.json`. These files are written only
after the run's outputs have been saved. An empty but complete result means
that no ORFs passed the selection rules.

The `<basename>.run.json` records the sample's settings, accepted read counts,
alignment filtering policy and exclusion diagnostics, input paths, software
versions and ORF count. Its `input_files` entries include
SHA-256 fingerprints and file sizes, which help identify the exact input files
used. A fingerprint identifies file contents; it does not assess data quality.
For a batch, `config.resolved.tsv` also records the configuration with defaults
filled in. Retain these files with the analysis.

## Choose candidates to inspect

A practical first pass through the report is:

1. Use **Search candidates** to find a locus tag or coordinate, or use
   **Category** to focus on a class such as `Unannotated`.
2. Use **Called support** to select TIS, TTS or both, and **Called sample** to
   focus on a particular sample. Selecting both means the same ORF has calls
   from both assay types; check the details to see which replicates contributed.
3. If you enabled statistics, use **Local statistics** to find candidates with
   at least one supported boundary. Open **Details** to confirm which assay and
   boundary carry that support. It can be an inferred boundary without a called
   peak in that assay.
4. Compare peak heights, q-values and support in each replicate. If the batch
   declared a matched-condition comparison, inspect its exploratory fixed-margin
   support, common odds ratio and pair directions as a separate line of evidence.
   This read-depth result is not a population-level replicate test.
5. Inspect the annotation and WIG coverage at the candidate's coordinates.
   Look for the expected endpoint peak, neighboring peaks, gene-body coverage
   where RIBO data are available, and plausible in-frame boundaries.
6. Record the candidate's `Identifier`, the samples supporting it and any
   interpretation notes in your own shortlist. Use its sequence or GFF entry
   for subsequent analyses.

Change **Sort by** to reorder candidates and **Rows per page** to show more
rows. Report filters help explore the saved results; they do not change the
underlying tables or recompute q-values. For a reusable filtered table, open
its XLSX file and save a filtered copy under a new name.

### A worked example: peak calls, q-values and inferred boundaries

Consider this **illustrative**, positive-strand ORF at `chromosome:100-162:+`
from a combined TIS/TTS run with `fdr=0.05`:

| Measurement | TIS-WT-1 | TTS-WT-1 |
| --- | --- | --- |
| Called peak height | 12 | missing |
| Boundary q-value | 0.004 | 0.02 |
| Local statistical support | true | true |

The TIS peak proposed the start at position 100. Following its reading frame
to the next allowed stop gives the end at 162. The TTS sample did not contribute
a called peak for this ORF, so the merged table's `Evidence` contains
`TIS-WT-1` only. It is a **TIS-called ORF with an inferred stop**.

The TTS q-value describes raw endpoints at that inferred stop's codon window.
It can be supported even though no individual position passed the peak-height
threshold: the statistical test and peak caller use different rules. This row
therefore passes a filter for at least one supported boundary, but it does not
pass a filter requiring both TIS and TTS calls. If the TTS peak height were
present too, `Evidence` would list both assays. Neither pattern by itself proves
translation or establishes biological function.

A q-value of 0.004 is an adjusted p-value for a local enrichment test. It is not
a 99.6% probability that this ORF is real. See [the statistical guide](statistics.md)
for the testing family, thresholds and assumptions.

## Open and filter the spreadsheet

For batch results, open `<experiment>_final.xlsx`, sheet **CDS**, in Excel,
LibreOffice or another spreadsheet application. For one sample, open
`table_per_sample/<basename>.xlsx`. The accompanying `.csv` files contain the
same tables but are **tab-separated**, despite their extension. Select **Tab**
as the delimiter if importing a `.csv` manually.

Useful filters include `Type`, `Genome`, `Strand`, `Evidence` in merged results,
per-sample peak-height columns, per-sample `peak_significant` flags and any
declared matched-comparison fields. To select ORFs called in both TIS replicates,
require positive values in both `TIS-WT-1_peak_height` and
`TIS-WT-2_peak_height`. Requiring both corresponding `peak_significant` fields
adds a user-defined consistency filter; it is not itself a new statistical
test across replicates.

### Location, sequence and annotation columns

| Column | Meaning |
| --- | --- |
| `Identifier` | Coordinate identity, such as `chromosome:100-162:+`; parsing is from the right, so contig names may themselves contain `:` or `-` |
| `Genome` | Reference contig identifier |
| `Start`, `Stop`, `Strand` | 1-based, inclusive genomic coordinates and translation strand |
| `Reference_length` | Physical length of the candidate's FASTA sequence |
| `Is_circular` | Whether the input GFF3 explicitly declared the whole reference circular |
| `Type` | Relationship to the supplied CDS annotation; categories are below |
| `Locus_tag` | Associated annotation display label, or the candidate's coordinate identity when no gene is assigned; repeated labels at distinct loci remain distinct internally |
| `Codon_count` | ORF length in codons, **including the stop codon**; multiply by 3 for nucleotide length |
| `Start_codon`, `Stop_codon` | Boundary codons in translation orientation |
| `Nucleotide_Seq` | ORF DNA in 5′ to 3′ translation direction, including its stop codon |
| `Amino_Acid_Seq` | Translation with the chosen code, including a terminal `*` for the stop |
| `15nt_window` | Up to 15 bases immediately upstream of the translation start, in transcript direction; shorter at reference edges |
| `5'-distance`, `3'-distance` | Distances from the proposed translation start to the associated annotation's 5′ and 3′ ends; these are not measured UTR lengths |

`Start` is always the lower coordinate of the linearized interval. On the
reverse strand, translation begins at `Stop`, the higher coordinate. For
example, the ORF
`chromosome:100-162:-` begins with a start codon spanning 160–162 and ends with
a stop codon spanning 100–102. It is still 63 nucleotides long.

### Circular-reference coordinates

Circular candidates use a bounded virtual interval so an origin crossing is
unambiguous: `1 <= Start <= Reference_length`, `Stop >= Start`, and a crossing
candidate has `Stop > Reference_length`. Reduce a virtual position `p` to its
physical reference position with `((p - 1) % Reference_length) + 1`. For
example, `Start=4900`, `Stop=5100`, `Reference_length=5000` covers physical
4900–5000 followed by 1–100. The identifier retains the same virtual interval.
This is an output convention only; source GFF3 feature rows must stay within
physical `1..Reference_length` coordinates.

Sequence extraction, codon searches, offsets, annotation overlaps and
expression overlaps wrap across the origin. A complete ORF is bounded to at
most one revolution. The `15nt_window` also wraps, but never repeats a physical
base when a reference is shorter than 15 bases. Circular multipart annotations
sharing a GFF3 `ID` are treated as one logical CDS. Linear references keep the
ordinary edge-clipped behavior.

| `Type` | How to interpret it |
| --- | --- |
| `Annotated` | Coordinates exactly match a CDS on the same contig and strand |
| `Annotated_Complex` | Coordinates exactly match a multipart CDS or one marked with a ribosomal-slippage exception; this label preserves annotation complexity but does not model a programmed frameshift |
| `Near_Annotated` | Same annotated stop, with a translation start less than 10 bases from the annotated start |
| `Internal_Inframe` | Inside an annotated CDS, overlapping a real CDS segment and in its reading frame |
| `N-terminal_extension` | Extends upstream of an annotated start while sharing its stop |
| `Internal_OutofFrame` | Inside an annotated CDS, in another reading frame |
| `Unannotated` | No relationship in the categories above was assigned |

`Unannotated` does not mean that a candidate lies entirely outside all annotated
genes, or that it is a newly validated gene. Inspect overlapping annotations,
strand and coverage. A candidate lying only in the genomic gap between parts of
a multipart CDS is not classified as internal to that CDS. Categories depend on
the annotation you supplied.

### Peak and statistical columns

In these column names, `<sample>` is a library name such as `TIS-WT-1`.

| Column | Meaning |
| --- | --- |
| `<sample>_peak_height` | Height of a called codon-window peak, using the run's normalization and `max` or `sum` operator |
| `Evidence` | Merged tables: TIS/TTS samples contributing a called peak for this ORF |
| `<sample>_relative_density` | Peak height divided by signal over the associated annotated gene, where available |
| `<sample>_peak_count` | Raw endpoints in the boundary's peak window for local statistics |
| `<sample>_background_count` | Raw endpoints in the two adjacent background flanks |
| `<sample>_peak_fold_enrichment` | Width-normalized peak endpoint density divided by local-background density |
| `<sample>_peak_p_value` | One-sided local enrichment p-value |
| `<sample>_peak_q_value` | P-value after the selected multiple-testing correction |
| `<sample>_peak_significant` | Local support: q-value at or below the configured threshold and enrichment greater than 1 |
| `<sample>_peak_correction` | `by` or `bh` correction |
| `<sample>_peak_tests_in_family` | Number of eligible codon windows tested for this sample and assay |
| `<sample>_peak_fit_status` | `ok`; `p_value_floored` when a positive tail is below float64 output range; `numerical_limit` when a count, bounded fallback or probability backend could not produce a safely supported tail and p=1 was retained; or `correction_limit` when a valid raw p-value was retained but q=1 was used because correction failed |

The peak caller keeps positions strictly above `min_peak_height`; `max` takes
the highest passing position and `sum` adds passing positions. The statistical
`peak_count` adds **all raw endpoints** in the window. Thus `peak_height` and
`peak_count` can differ even with raw normalization. Local fold enrichment is a
descriptive effect size, not a population effect: it is infinite when a positive
peak has zero background endpoints, zero when only the background has endpoints,
and missing when the whole tested region has no endpoints or no background bases.

### Exploratory fixed-margin matched columns

When a batch configuration declares matched comparisons, merged tables include
fields prefixed by `<comparison>_<assay>`, for example
`treated_vs_control_TIS_q_value`. These fields describe a conditional endpoint-
count comparison, not a biological-dispersion or population-effect model. The
most useful fields are:

| Column suffix | Meaning |
| --- | --- |
| `_matched_pair_count` | Number of matched replicate IDs |
| `_informative_pair_count` | Pairs that contribute a non-fixed conditional table |
| `_concordant_pair_count` | Informative pairs favoring the numerator direction |
| `_discordant_pair_count` | Informative pairs favoring the denominator direction |
| `_tied_pair_count` | Informative pairs with equal observed cross products |
| `_common_odds_ratio`, `_log2_common_odds_ratio` | Descriptive matched common effect |
| `_p_value`, `_q_value` | One-sided fixed-margin result and its separately adjusted value |
| `_significant` | Exploratory support: q threshold, common OR above 1 and a strict majority of informative pairs favoring the numerator |
| `_fit_status` | Whether the test was computed, lacked two informative pairs, used a numerical probability floor, exceeded the declared hypergeometric bound or lacked safe probability intermediates, exceeded the exact-convolution safety limit, or retained a valid raw p-value with q=1 after correction failed |

Raw numerator and denominator peak/background totals are also retained. Open
`matched_statistics/<comparison>_<assay>.tsv` for all tested codon windows,
including candidates not attached to a called ORF. At least two informative
pairs are required for p/q/support values; missing values are not negative
results. Eligible candidates that hit a numeric-range or exact-work limit remain
in correction with p=1 and cannot provide support. Each comparison and assay is
corrected separately, so selecting support from any of several families has no
joint FDR guarantee.

The aggregate matched table's raw counts are totals. Its companion
`matched_statistics/<comparison>_<assay>_strata.tsv` has one row per candidate
and pair with `replicate_id`, all four peak/background counts, `informative` and
`pair_direction`, plus the assay and condition labels. Join the two tables by
`candidate_id` to inspect heterogeneity
or reconstruct the fixed-margin probability. Pair IDs must represent the real,
predeclared experimental blocks; arbitrary ordinal matching or technical file
splits are not independent replication. The common odds ratio can hide different
magnitudes or directions across pairs.

The run JSON records the model and alternative, replicate IDs, minimum pair and
support rules, numerical limits, FDR settings and paths to both matched tables.
Use these values when comparing or reproducing analyses rather than inferring
the rules from a `significant` cell alone.

This model compares local enrichment between conditions within one assay. Its
p-value is driven by endpoint depth under fixed margins; it does not estimate
biological dispersion, test a population-level condition effect, test total
abundance, use RIBO as a control, or combine TIS and TTS significance. Input
libraries need comparably calibrated read-length and offset choices that place
endpoints at the same biological boundary. See
[matched-condition inference](statistics.md#matched-condition-inference).

## Understand missing values

A missing measurement may appear as an empty spreadsheet cell, `NA` or an
unavailable indicator in the report. Keep these separate from measured zeros.
Integer audit counts outside JavaScript's exact numeric range are retained as
decimal text in the HTML report, so their displayed digits match the canonical
tab-separated output instead of being silently rounded by the browser.

| Value | Interpretation |
| --- | --- |
| Missing peak height | This sample did not contribute a called peak for this ORF, or the assay was not supplied |
| Missing q-value or support flag | Statistics were off, the sample was absent, or no eligible candidate window was available for the boundary |
| Missing local fold enrichment | The boundary was absent, the tested region had no endpoints, or no background bases were available; it is not fold enrichment zero |
| `peak_significant = false` | The boundary was tested and did not satisfy the local-support rule |
| Local `fit_status = p_value_floored` | The exact tail is positive but below float64 output range; the exported p/q values remain positive |
| Local `fit_status = numerical_limit` | The eligible window remains in correction with p=1 and cannot provide support because its count, bounded tail fallback or numerical backend could not produce a safely supported tail |
| Local `fit_status = correction_limit` | Correction failed; the valid raw p-value remains auditable, but q=1 is used and the window cannot provide support |
| Missing matched p/q/support | The comparison was absent from this boundary or fewer than two pairs had informative margins |
| Matched `fit_status = numeric_range_limit` or `computational_limit` | The eligible candidate remains in its correction family with p=1 and cannot provide support |
| Matched `fit_status = correction_limit` | Correction failed; the valid raw p-value remains auditable, but q=1 is used and the candidate cannot provide support |
| Matched `significant = false` | The candidate was tested but did not meet the exploratory q/effect/strict-majority support rule |
| q-value of 1 | No usable adjusted evidence under the test; also used for windows with no endpoints, no available background, a declared conservative numerical limit, or correction failure |
| Missing RPKM | Expression counting was not supplied for that sample/result |
| RPKM of 0 | No overlapping accepted reads, or the supplied sample has a zero normalization denominator |
| Missing TE | No matched RNA/assay pair was available, or the underlying matched RNA rate was zero because it had no overlapping read or no accepted-read denominator |

A failed support rule does not demonstrate absence of translation. Low read
coverage, an unsuitable local background, an incorrect offset or an inferred
boundary can affect the evidence. Do not replace all missing cells with zero
when comparing samples.

## Interpret expression and ratios

RPKM counts accepted alignments overlapping the ORF on its strand, accounting
for ORF length and accepted-read total. Each alignment is counted once even
when its CIGAR produces multiple overlapping blocks. `normalization_scope`
sets whether the denominator covers that contig or the whole library.
`rpkm_read_usage`/`--all_reads_rpkm` controls whether all lengths or selected
lengths contribute. Compare RPKM values only with these choices in mind.
Paired mates are independent alignment records; there is no fragment-counting
mode.

TE is the assay RPKM divided by the matched RNA RPKM for the same condition and
replicate: RIBO/RNA, TIS/RNATIS, or TTS/RNATTS. It is a descriptive ratio; it has
no associated replicate-aware significance test. TE uses the unrounded
read-count/normalization-denominator rates for the selected contig or library
scope; RPKM columns are displayed to two decimal places. A very low positive
RPKM can therefore display as zero while still contributing a finite, auditable
TE instead of being treated as absent.

Merged tables with matched RIBO data can also contain
`<TIS-or-TTS-sample>_<RIBO-sample>_log2FC`. This is the log2 ratio of peak
heights. When RIBO height is zero or absent, the calculation substitutes 0.9
times the smallest positive peak height in the merged table. A large value
can therefore reflect that fallback. Check the original heights before
interpreting the ratio; it is not a treatment-versus-control p-value.

## Inspect candidates in a genome browser

1. Load the same reference FASTA used for the analysis into your browser.
2. Add the original annotation GFF and the candidate GFF:
   `<experiment>_final.gff` for a batch or `gff_per_sample/<basename>.gff` for a
   single sample.
3. Add relevant WIG files from `coverage_files/<basename>/`. Names identify
   the sample, contig and `forward` or `reverse` strand. Reserved filename
   characters in a contig identifier, including `/`, are percent-encoded in the
   filename; the WIG `chrom` declaration retains the exact reference name.
4. Navigate to the candidate's `Genome`, `Start` and `Stop` and inspect the
   matching strand. Use the same vertical scale when comparing compatible
   tracks, and check the run's normalization.

The WIG tracks show the selected mapping, offsets and normalization; they are
not unmodified BAM coverage. Their coordinates are always physical `1..L`.
GFF coordinates are 1-based and inclusive. A circular origin-spanning candidate
is emitted as two physical CDS rows with one shared `ID`; phase 0 is on the row
containing the biological translation start and the continuation row carries
the appropriate phase. Inspect both sides of the origin in browsers that do not
draw multipart circular features continuously. Topology declarations remain in
an empty GFF even when no candidate was called. Browser-generated indexes may
be needed for your FASTA or BAM. If a candidate appears at the wrong location,
first check contig names, topology declarations and the selected assembly.

## Use nucleotide and protein sequences

The exports include:

| File | Contents and use |
| --- | --- |
| `sequence_per_sample/<basename>.fna` | Candidate ORF DNA for one direct-run or per-replicate assay set; suitable for sequence searches or further nucleotide analysis |
| `sequence_per_sample/<basename>.faa` | Corresponding protein sequences for protein searches or annotation |
| `<experiment>_final.fna` | Batch ORF nucleotide sequences, merged by coordinate identity |
| `<experiment>_final.faa` | Corresponding merged protein sequences |

Sequences are in translation orientation on both strands. DNA includes the
terminal stop codon. The protein FASTA omits the terminal `*` that is shown in
`Amino_Acid_Seq`; valid alternative initiation codons are translated as the
initial methionine. A 63-nucleotide complete ORF therefore has 21 codons
including the stop and a 20-amino-acid protein sequence. Translation uses the
run's `genetic_code`.

Use the coordinate identifier in the FASTA header to find its table row.
An unannotated candidate has no fabricated description after that identifier.
If a result row lacks either required sequence, FASTA export stops with an
error instead of writing the strings `nan` or `<NA>` as biological sequence.
The report's **Details** view also shows selectable sequences for an individual
candidate. Keep the original identifier if another tool adds annotation so
you can join the results back to the ORF table.

## Keep or share an analysis

Archive the whole output folder to preserve relative report links and the
per-sample evidence. A batch folder contains:

```text
<experiment>/<mapping>/<normalization>/
├── report.html
├── files.html                       # links to coverage and statistical tables
├── complete.json
├── <experiment>_final.xlsx / .csv / .gff / .fna / .faa
├── config.resolved.tsv
├── <basename>.run.json
├── table_per_sample/
├── gff_per_sample/
├── sequence_per_sample/
├── coverage_files/<basename>/
├── peak_statistics/                 # present when local statistics are enabled
└── matched_statistics/              # aggregate and per-pair strata TSVs for declared comparisons
```

Save a separate shortlist with your filtering criteria and notes. Report
filters do not alter the original exports, and merging samples does not
recalculate statistical families. Preserve the run records, reference and
calibration files so collaborators can understand how the results were produced.
