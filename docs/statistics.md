# Use statistical support

Enable local statistics when you want additional evidence about whether a TIS
or TTS boundary has more read endpoints than its immediate surroundings. Use
this evidence to prioritize candidates for inspection alongside peak calls,
annotation, coverage and biological replicates.

The single-sample test describes **local endpoint enrichment**. It does not
establish translation or treat RIBO measurements as controls. Batch analyses
can additionally calculate exploratory fixed-margin support for a predeclared,
same-assay comparison across experimental pairs; that separate model is
described below.

## Enable the analysis

For a direct run, add `--statistics local` and use the calibrated endpoint
mapping method, `fiveprime` or `threeprime`:

```bash
orfbounder \
  --alignment_file_tis alignments/TIS-WT-1.bam \
  --annotation_file annotation.gff --genome_file genome.fa \
  --offset_json offsets.json --mapping_method fiveprime \
  --statistics local --background-width 50 --fdr 0.05 --fdr-method by \
  --output_path results/local-WT-1 --output_basename WT-1
```

For a batch, set these columns in your TSV:

| Column | Value in the command above | What to choose |
| --- | --- | --- |
| `statistics` | `local` | `none` disables statistical testing |
| `background_width` | `50` | Number of reference bases on **each side** of the peak window |
| `fdr` | `0.05` | Adjusted p-value threshold used for the local-support flag |
| `fdr_method` | `by` | Correction for multiple tests; `bh` is also available |

Set mapping, offsets, read lengths and background width before reviewing the
results. `centered` and `global` mapping are unsupported for statistics because
they distribute a read across multiple positions. The local test needs one
integer endpoint contribution per accepted read.

## Find and use the values

Open `WT-1.report.html` for the direct example, or `report.html` in the batch
experiment folder. Use **Local statistics** to find candidates with a supported
boundary, then open **Details** to inspect the sample, local fold enrichment,
q-value and whether that sample also contributed a called peak.

In the spreadsheet, the main fields are:

| Field, prefixed by the sample name | What it tells you |
| --- | --- |
| `peak_count` | Raw endpoints in this ORF boundary's codon window |
| `background_count` | Raw endpoints in the adjacent flanks on the same strand |
| `peak_fold_enrichment` | Peak-window endpoint density divided by adjacent-background endpoint density |
| `peak_p_value` | Evidence against a uniform local endpoint rate, before correction |
| `peak_q_value` | P-value adjusted across that sample/assay's tested codon windows |
| `peak_significant` | True when the q-value is at most your `fdr` threshold and enrichment is greater than 1 |
| `peak_correction` | Which correction was used: `by` or `bh` |
| `peak_tests_in_family` | Number of tested codon windows in the correction family |
| `peak_fit_status` | `ok`, `p_value_floored`, `numerical_limit`, or `correction_limit`; see the numerical conventions below |

For example, `TIS-WT-1_peak_q_value` concerns the start boundary in sample
`TIS-WT-1`. A TTS field concerns the stop boundary. A small q-value is not the
probability that the proposed ORF is false, and a false support flag is not
proof that the ORF is untranslated.

Statistics **do not automatically remove ORFs**. A called peak can lack local
statistical support. Conversely, a boundary inferred from the other assay can
receive statistical support even when no peak was called there. Check both
`peak_height`/`Evidence` and the statistical fields. The
[worked ORF example](results.md#a-worked-example-peak-calls-q-values-and-inferred-boundaries)
shows how these situations appear in a result row.

The full candidate tables are in:

```text
peak_statistics/<output_basename>_<sample>.tsv
```

They include every eligible codon window, including those that did not produce
an ORF. Open them as tab-separated text to inspect the original counts, window
coordinates, widths and `fold_enrichment`. In an ORF table that effect size is
named `<sample>_peak_fold_enrichment`; the candidate table's `p_value`,
`q_value`, `significant`, `correction` and `tests_in_family` likewise correspond
to sample-prefixed fields. Its `fit_status` becomes
`<sample>_peak_fit_status` in an ORF table.

## What the comparison measures

Each eligible start codon in TIS data, or stop codon in TTS data, has a five-base
peak window. The comparison uses immediately adjacent reference bases on the
same contig and strand, excluding the peak itself. With the default setting,
there are 50 background bases on each side, for a total of 100.

On a linear reference, peak and background windows are clipped at reference
edges. On a reference explicitly declared circular in GFF3, they wrap across
the origin. A physical base can contribute at most once even when both flanks
overlap on a small circle, so peak plus background width never exceeds the
reference length. Actual widths are always used. In a candidate TSV,
`background_width` is the **total available background width**, normally 100
when the configured `background_width` is 50 on a sufficiently long reference.
`Codon_start` is the codon's forward-reference start coordinate (physical
`1..L`). Window coordinates are 1-based inclusive virtual intervals;
`Window_stop` can exceed
`Reference_length` only for a declared origin crossing. Use `Is_circular` and
`Reference_length` when converting those values to physical coordinates.

The test uses raw endpoints after the same alignment filters, read-length
selection and offsets as peak calling. It includes endpoints at positions below
the peak-height threshold. Changing raw/mil normalization or its contig/library
scope does not change these statistical counts or results for the same input
endpoints.

For a numerical example, suppose a five-base window contains 12 endpoints and
its 100 background bases contain 20. The observed densities are 2.4 and 0.2
endpoints per base, so `fold_enrichment` is **12**. Under the uniform local model,
only 5/105 of the combined 32 endpoints would be expected in the peak window
on average. The p-value measures how unusual a count of at least 12 would be
under that model. Statistical support still depends on the adjusted p-value
across the full testing family; fold enrichment alone is insufficient.

For an extreme tail below the normal floating-point range, the software uses a
bounded recurrence with an exact integer starting combination and an upward
rounding safeguard. This reduces cancellation seen in the ordinary probability
backend; it is numerical engineering, not a formal exact-arithmetic guarantee.
If the fallback cannot be evaluated within its declared limits, the eligible
test stays in the correction family with p=1 and `numerical_limit` status.

If background counts are zero and peak counts are positive, fold enrichment is
infinite because no pseudocount is added. When the entire region has no reads,
or no background bases are available, fold enrichment is missing and the
p-value is 1. A peak count of zero with background reads gives fold enrichment
zero. Infinity therefore does not by itself imply support, and a missing value
must not be replaced by zero.

## Multiple testing for local evidence

One **sample and assay** defines a testing family. It contains all eligible
codon windows across every reference contig and both strands, including windows
with zero reads. The family is defined before height filtering and ORF-boundary
searches. The selected genetic code and start/stop codon sets determine which
codons are eligible.

The default Benjamini–Yekutieli (`by`) correction accommodates dependence among
tests, provided their individual p-values are valid. Benjamini–Hochberg (`bh`)
is less conservative and requires independent tests or suitable positive
dependence for its FDR guarantee. Both procedures are described in
[SciPy's false discovery control documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html).

ORFBounder calls the resulting adjusted p-value a **q-value**. The local-support
flag uses `q_value <= fdr` together with enrichment above 1. The flag is evidence
under this particular model; it is not a demonstrated error rate for biological
ORF annotations.

SciPy's BH/BY arithmetic can undershoot the exact adjustment of its binary64
input p-values by one or two representable steps. Both local and matched paths
therefore validate the returned family and apply the same small, family-scaled
upward rounding guard. A result exactly on such a numerical threshold may
conservatively lose its support flag. This is an engineering tolerance verified
against exact-fraction regressions, not a formal interval-arithmetic proof.

Batch merging preserves each sample's original family and q-values. It neither
pools these read counts nor recomputes these p-values across replicates. A report
filter selecting ORFs with local support in any sample is an exploration aid;
it does not control a new false discovery rate for that combined list. Selecting
favorable results across multiple offsets, background widths, samples or
parameter runs also introduces comparisons outside the reported family.

## Matched-condition inference

A batch run can calculate **exploratory fixed-margin matched support** for one
condition against another within the same TIS or TTS assay. The comparison is
directional and must be declared before inspecting results. It retains each
predeclared experimental pair or block as a separate stratum instead of pooling
its margins with other pairs.

Copy [the comparison template](../templates/template_matched_comparisons.tsv),
or create a tab-separated comparison file with exactly these columns:

```text
comparison	assay	numerator_condition	denominator_condition
treated_vs_control	TIS	treated	control
```

Then set these batch configuration fields:

| Column | What to choose |
| --- | --- |
| `statistics` | `local`; matched inference uses its raw candidate counts |
| `matched_comparison_file_path` | Path to the comparison TSV |
| `matched_fdr` | Q-value threshold, default `0.05` |
| `matched_fdr_method` | `by` by default, or `bh` |

For the example, filenames must include both `TIS-treated-<replicate>` and
`TIS-control-<replicate>`. The two conditions must have exactly the same set of
at least two replicate IDs. An ID must denote a real pair or block defined by
the experimental design: numerator and denominator samples carrying that ID
must have been intentionally paired, and different IDs must be independent
experimental units. Do not give unrelated samples matching ordinal labels just
to satisfy the filename rule. Technical lane or file subdivisions of one
library do not become independent pairs. More independent pairs improve the
consistency assessment, but do not make this a biological-dispersion model.

The assay column accepts `TIS` or `TTS`. Both sides of a comparison always use
that assay: a RIBO library is not silently used as a TIS or TTS control. Add
separate rows to test separate assays. ORFBounder does not combine their
p-values into a cross-assay score.

For each candidate and replicate, the model forms this 2 × 2 table from raw
endpoint counts:

|  | Peak window | Local background |
| --- | ---: | ---: |
| Numerator condition | `a` | `b` |
| Denominator condition | `c` | `d` |

It conditions on the row and column margins, making `a` hypergeometric under
equal peak-to-background odds. ORFBounder convolves the stratum distributions
and reports the one-sided fixed-margin probability of a numerator peak-count sum
at least as large as observed. Thus every pair retains its own local depth.
Library normalization is neither estimated nor applied to these raw tables.

The fixed-margin law is evaluated numerically: a mode-centered hypergeometric
probability recurrence reduces cancellation for large margins, and a small
upward rounding guard avoids the downward bias seen in exact-rational checks.
This is an engineering safeguard, not a formal exact-arithmetic guarantee;
unsafe numerical or work-limit cases retain a conservative p-value of 1.

The sampling units in this conditional probability are read endpoints, not the
biological pairs themselves. More endpoint depth can therefore make the p-value
smaller even when the number of pairs is unchanged. The calculation assumes the
fixed-margin count model within pairs and independence between pairs; it neither
estimates extra-binomial biological variation nor tests a population-level
condition effect.

The output table is
`matched_statistics/<comparison>_<assay>.tsv`. The merged ORF table receives
the same result fields with `<comparison>_<assay>_` prefixes:

| Field | Meaning |
| --- | --- |
| `matched_pair_count` | Number of replicate pairs in the declared comparison |
| `informative_pair_count` | Pairs whose margins allow peak/background odds to vary |
| `concordant_pair_count` | Informative pairs with greater observed odds in the numerator |
| `discordant_pair_count` | Informative pairs with greater observed odds in the denominator |
| `tied_pair_count` | Informative pairs whose observed cross products are equal |
| `numerator_peak_count`, `numerator_background_count` | Raw counts summed for display across numerator replicates |
| `denominator_peak_count`, `denominator_background_count` | Corresponding denominator totals |
| `common_odds_ratio` | Descriptive Mantel–Haenszel common odds-ratio estimate |
| `log2_common_odds_ratio` | Base-2 transform of that estimate |
| `p_value`, `q_value`, `significant` | Fixed-margin directional result, adjusted result and exploratory support flag |
| `correction`, `tests_in_family` | Correction and number of eligible candidates with at least two informative pairs in this comparison |
| `fit_status` | `ok`; `no_informative_pairs` or `insufficient_informative_pairs` when fewer than two pairs can contribute; `p_value_floored` for a probability below floating-point output range; `numeric_range_limit` when a stratum exceeds the declared hypergeometric bound or its probability backend cannot produce safe finite terms; `computational_limit` when an exact convolution would exceed the safety limit; or `correction_limit` when correction fails. Numeric/work limits retain p=1; a correction limit retains the valid raw p-value but uses q=1 and cannot support a candidate. |

Despite its historical column name, `significant` is an exploratory support
flag. It requires all three of the following: the q-value threshold, a common
odds ratio greater than 1, and a strict majority of informative pairs favoring
the numerator (`2 * concordant_pair_count > informative_pair_count`). Ties do
not favor either direction. At least two informative pairs are required before
a p-value is calculated. Candidates with fewer retain missing p/q/support
values rather than zeros and are excluded from the correction family.

The common-effect direction around one is evaluated as the exact sign of the
sum of the strata's rational cross-product differences. An exactly neutral
effect therefore cannot round above one and pass the support gate; the displayed
common odds ratio is kept consistent with that exact direction.

Each comparison and assay is corrected as its own family. Every candidate with
at least two informative pairs stays in that family; a numeric-range or
convolution-work failure contributes p=1 and cannot provide support. This avoids
letting data-dependent computational geometry shrink the correction. Selecting
candidates supported by any of several comparisons or assays is a new selection
rule and has no joint FDR guarantee. If correction itself fails, valid raw
p-values remain available for audit, but every q-value becomes 1 and the family
cannot provide support. The
strict-majority requirement is a descriptive consistency safeguard, not a
separate replicate-level p-value or a biological error-rate guarantee.

This model tests a change in **peak relative to its adjacent background**. It is
not a negative-binomial test of total expression, does not estimate biological
dispersion or a population effect, and does not establish translation. The
Mantel–Haenszel odds ratio summarizes a common effect; different pair magnitudes
or directions can make that summary unsuitable, and ORFBounder does not fit an
effect-heterogeneity model.

Pairing also does not eliminate condition-specific mapping or protocol biases.
Use the same reference, codon definitions, background geometry, endpoint
mapping convention and alignment policy across conditions. Read-length choices
and offsets should be comparably calibrated for each library so they place
endpoints at the same biological boundary; their numeric values need not be
identical when library-specific calibration supports a difference.

The aggregate matched table contains counts summed across pairs for display.
Its companion
`matched_statistics/<comparison>_<assay>_strata.tsv` contains one row per
candidate and pair: `replicate_id`, the four numerator/denominator peak and
background counts, `informative`, and `pair_direction`; assay and condition
labels are repeated so the file remains self-describing. Join it to the aggregate
table by `candidate_id` to audit directions or reconstruct the conditional
test. Keep both tables, the per-sample candidate tables, comparison TSV,
resolved configuration and run JSON together.

Each run JSON's `matched_comparisons` record identifies the `model`, one-sided
`alternative`, numeric backend, exact `replicate_ids`, `minimum_informative_pairs`, textual
`support_rule`, hypergeometric grand-total bound, convolution work limit,
minimum representable output p-value,
FDR settings, and paths to both candidate tables. These fields make numerical
floors and candidates conservatively retained with p=1 at a numerical or
computational limit distinguishable from ordinary fitted results.

Each informative 2 × 2 stratum has a deliberately conservative engineering
limit of `2^31 - 1` on its grand total `a + b + c + d`. This bound is checked
before calculating the floating-point common effect or calling SciPy, so
accepted behavior does not depend on a platform's integer width. An otherwise
eligible candidate with a larger stratum remains in the correction family with
`p_value=1`, missing common-effect fields and
`fit_status=numeric_range_limit`. Its raw per-condition totals and companion
stratum counts remain exact Python integers in tab-separated audit output, even
when their sums exceed signed 64-bit range. Spreadsheet cells above binary64's
exact-integer range (`2^53`) are stored as decimal text rather than silently
rounded. This is a software-safety bound,
not a claim about an experimental depth at which the fixed-margin model changes.

## Check the evidence before drawing conclusions

Inspect local coverage and neighboring peaks, particularly near gene or
transcript boundaries. Nearby peaks count as background, and the windows do
not follow gene boundaries. Endpoint periodicity, transcript structure,
sequence/mapping bias, low mappability, read dependence, PCR duplication and
overdispersion can make the uniform local model unsuitable. The configured
MAPQ, duplicate and `NH` multimapper policy is applied consistently to peak
calls, statistical counts, expression counts and recovered normalization
totals. Excluding duplicate-flagged records does not guarantee that every
technical duplicate was identified.

BY correction addresses dependence between valid tests; it cannot repair an
unsuitable background model. Interpret local support with your annotation,
controls, independent evidence and consistency across biological replicates.
A stronger biological claim requires evidence beyond this enrichment result.

## Statistical definition

Let `X` be the raw peak endpoint count, `B` the raw background count, and `Lp`
and `Lb` the actual numbers of peak and background bases. Conditional on
`N = X + B`, a homogeneous local endpoint rate gives the null distribution:

```text
X ~ Binomial(N, Lp / (Lp + Lb))
p_value = P[Binomial(N, Lp / (Lp + Lb)) >= observed X]
fold_enrichment = (X / Lp) / (B / Lb)
```

This is a one-sided conditional binomial test, equivalent to the `greater`
alternative described in [SciPy's binomial test documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binomtest.html).
The no-read and no-background conventions above define cases in which the
comparison cannot provide local enrichment evidence.

The enrichment direction used by the support flag is the exact integer
comparison `X * Lb > B * Lp`, not a comparison of rounded display values. An
exactly neutral ratio is displayed as 1.0; values extremely close to one are
kept on the same side as that exact comparison.

Mathematical binomial tails are positive, but extremely small tails can
underflow in floating-point software. ORFBounder re-evaluates such tails in log
space. If the result is below the smallest positive float64 value, it exports
that minimum positive value rather than zero and sets `fit_status` to
`p_value_floored`. Adjusted q-values are likewise never written as zero. This
floor is conservative at the available output precision and does not make a
p-value more exact.

For ordinary direct tails, the exact width fraction can itself fall between
binary64 values and the numerical backend can undershoot by a few output ULPs.
ORFBounder evaluates the null probability one representable step upward and
applies a small trial-scaled upward engineering guard to the returned tail. A
bounded exact-rational regression covers production peak widths, representative
background widths and all possible observed counts. The same upward width
conversion is used by the log-space fallback. These safeguards are deliberately
conservative around a configured threshold; they are not a formal
interval-arithmetic guarantee.

The binomial trial count has a deliberately conservative engineering support
limit of `2^31 - 1`. This is far above realistic local endpoint depth; it is not
a theorem about SciPy accuracy. A larger total, a fallback that cannot finish
within its declared recurrence or exact-combination work bounds, or a numerical
backend that cannot produce safe finite tail terms retains the eligible
candidate in the correction family with `p_value=1`, `significant=false` and
`fit_status=numerical_limit`. A multiple-testing backend failure instead keeps
the valid raw p-values for audit, sets every q-value to 1, prevents support and
uses `fit_status=correction_limit`. This deliberately discards adjusted evidence
instead of risking an anti-conservative value or shrinking the multiple-testing
family. The run JSON records the local model, one-sided
alternative, positive probability floor, supported total-count bound and both
fallback work bounds. Raw counts above signed 64-bit range remain exact Python
integers in attached tab-separated output; ordinary counts retain nullable
`Int64` columns. The standalone merger parses local and matched integer evidence
as decimal integers rather than through floating point and rejects fractional,
scientific-notation or non-finite count fields.
