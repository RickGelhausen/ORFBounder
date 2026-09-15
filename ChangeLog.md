### Unreleased maintenance and container update

 + Rejected ambiguous reference inputs before analysis: FASTA and alignment
   reference identifiers must be nonempty and unique, primary mapped records
   must provide CIGAR, and GFF3 feature rows cannot repeat or omit attribute
   names. Multi-parent annotations remain supported through one comma-separated
   `Parent` value.
 + Normalized nullable local and matched support flags consistently when
   re-importing TSV or XLSX results, and rejected values outside the documented
   Boolean/blank schema instead of retaining unsupported scalar values.
 + Made standalone result-to-GFF import preserve literal and numeric-looking
   identity labels, reject duplicate headers and detect contradictions between
   canonical coordinate identifiers and their coordinate columns.
 + Prevented floating-point threshold crossings at exactly neutral effect
   boundaries: local enrichment uses an exact count/width cross product and
   matched common-effect direction uses an exact reduced-rational sign. Displayed
   ratios are kept consistent with those support decisions.
 + Rounded direct local-binomial null probabilities and tails conservatively;
   representative and wider exact-rational sweeps cover the observed SciPy
   few-ULP undershoots that could otherwise cross a one-test-family threshold.
 + Required every calling sample to appear in supplied `min` normalization
   counts, preventing omissions or misspellings from silently changing the
   normalization target while still allowing intentional reference samples.
 + Inferred query length from CIGAR when a valid SAM/BAM record omits `SEQ`, so
   read-length filters, offsets, expression counts and recovered denominators
   use the same length.
 + Kept candidates lying wholly in gaps between multipart CDS segments
   unannotated rather than classifying them as internal to the outer bounds.
 + Rejected local-statistics attachment onto an existing sample/family schema,
   preventing silent replacement of p/q values and testing-family metadata.
 + Applied one shared, small family-scaled upward rounding guard to local and
   matched BH/BY adjusted values, covering observed one-to-two-ULP backend
   undershoots while keeping threshold decisions conservative.
 + Omitted fabricated `<NA>` descriptions from unannotated FASTA headers and
   rejected missing sequence cells instead of exporting `nan`/`<NA>` as sequence.
 + Preserved leading zeros in numeric-looking reference IDs and annotation
   labels when re-importing TSV/XLSX results for merging and GFF export.
 + Rejected malformed or nonpositive `NH` alignment tags rather than silently
   treating them as unique alignments or failing with an opaque type error.
 + Removed a newly created empty output directory after a failed staged run
   only while its inode still belongs to that run; existing and changed
   directories remain untouched.
 + Reduced subnormal local binomial-tail cancellation using a bounded exact
   starting combination and an upward rounding guard. Both the ordinary and
   rare-fallback paths bound a rounded peak-window width ratio upward.
 + Reduced cancellation in matched-condition hypergeometric tails with a
   mode-centered probability recurrence and an upward floating-point rounding
   guard, while retaining conservative numerical/work limits and an exact
   neutral-boundary direction check for the common effect.
 + Required canonical read-length keys in offset JSON, rejected duplicate
   read-length/offset JSON keys and bounded expansion of accidental enormous
   read-length ranges before allocation.
 + Preserved literal annotation labels such as `NA`, `NULL` and `N/A` when
   re-importing TSV/XLSX result tables for merging, while retaining blank
   cells as missing measurements.
 + Preserved raw statistical audit integers across nullable TSV and Excel
   exports and through merges; offline reports encode values outside
   JavaScript's exact integer range as decimal text rather than rounding them.
 + Rejected Boolean or malformed local FDR thresholds and duplicate local
   boundary rows, preventing accidental support or overwritten evidence.
   Direct-run numeric settings now fail before output when they are not valid
   integers or finite real thresholds.
 + Kept transcript-mediated GFF3 `Parent` lineage and percent-encoded list
   members intact, validated every multipart CDS row, clipped linear terminal
   codon tracks, and made WIG filenames safe for sequence IDs with separators.
 + Aligned the standalone GFF helper's circular and legacy-linear topology
   handling with the main exporter, while preserving exclusive publication.
 + Tightened output recovery around symlinked ancestors, `..` aliases and
   untrusted stage controls; linked or nonregular lock/journal entries are not
   followed, and unrelated lookalike stages are not modified.
 + Closed FASTA and test reader handles deterministically and verified the
   repository suite with warnings promoted to errors.
 + Added exploratory fixed-margin, directional TIS/TTS condition comparisons
   for predeclared independent pairs. Support requires at least two informative
   pairs and strict-majority agreement; BY/BH families remain separate, and no
   biological-dispersion or population-effect claim is made. Per-pair strata
   are exported for audit, and exact calculations have explicit work and
   numeric-range safeguards; eligible limit cases remain in correction with a
   conservative p-value of 1. Hypergeometric grand totals are bounded before
   SciPy evaluation so oversized Python integers behave consistently across
   platforms.
 + Made local and matched statistical backend failures conservative through
   multiple-testing correction. Invalid support/probability arrays and unsafe
   adjusted values cannot escape; correction failures retain valid raw
   p-values but use q=1, an explicit `correction_limit` status and no support.
 + Added explicit MAPQ, duplicate and `NH`-multimapper policies with exclusion diagnostics; paired mates remain separate alignment records.
 + Added whole-library `min` normalization from per-sample totals.
 + Added opt-in circular-reference calling, wrapped statistics/counting and
   physical multipart GFF3 export with explicit topology metadata.
 + Added safe multipart annotation identities, repeated display labels and coordinate identifiers containing delimiters.
 + Indexed annotation candidates by chromosome, strand and interval instead of
   scanning every gene for every ORF, while preserving classification priority
   and circular virtual-coordinate matching.
 + Added recoverable publication journals, a cross-filesystem exclusive-copy fallback,
   final content re-verification of inputs and publication-boundary rejection
   of symlinked output roots.
 + Hardened publication and recovery with no-follow directory operations,
   immediate ownership journaling and retained journals after incomplete
   cleanup. Standalone merges and multi-experiment mapped-count recovery now
   publish their companion files as all-or-nothing sets, and count recovery
   rechecks content fingerprints before publication.
 + Added offline candidate reports, nucleotide/protein FASTA exports and completion records.
 + Reused immutable, content-keyed FASTA/GFF reference parsing across batch
   parameter grids and repeated configuration rows while keeping mutable
   density state isolated.
 + Reduced peak memory during offline-report creation by streaming its JSON
   payload; report search no longer retains a duplicate lower-cased table.
 + Serialized report candidates one row at a time instead of materializing a
   second prepared table, substantially reducing generation memory for large
   reports while preserving the same embedded JSON.
 + Exposed per-sample local fold enrichment alongside q-values so effect size
   and read-depth evidence can be inspected separately.
 + Prevented extreme local binomial tails from being exported as zero; added
   explicit probability-floor and conservative numerical-limit statuses and
   recorded their bounds in run metadata.
 + Made repeated-table merging order-independent: missing measurements coalesce,
   while conflicting measurements or shared ORF metadata are rejected.
 + Kept annotation strings literal in Excel, preserved tabs/newlines in quoted
   direct and merged tabular fields, and suppressed ANSI colors in captured
   logs or with `NO_COLOR`.
 + Required supplied read-length filters to cover every calling sample (and
   length-filtered expression sample) by exact key or a top-level default in
   direct, batch and count-recovery workflows, preventing sample-name typos
   from silently enabling all read lengths.
 + Required offset JSON coverage for every direct-run calling sample before
   any assay writes coverage, matching batch preflight behavior.
 + Calculated translation-efficiency ratios from unrounded read-count rates at
   the selected normalization scope while retaining two-decimal RPKM display,
   so low positive expression is not rounded into a missing ratio.
 + Vectorized merged TIS/TTS-to-RIBO fold-change contrasts while preserving
   missing, zero, negative and infinite-value behavior and output ordering.
 + Packaged the documented helper modules in wheels and made the toy-example
   invocation work from installed distributions.
 + Added a locked Ruff development dependency and CI linting; container CI now
   exercises both command-line entry points.
 + Removed unused combined-search and workflow-specific mapping-name helpers.
 + Added selectable genetic codes and whole-library RPM/RPKM scaling.
 + Protected output publication and recorded input fingerprints, including alignment aliases.
 + Reworked user guides around running analyses, inspecting evidence and using exported results.
 + Added installable commands, Python 3.11+ packaging and a standalone Docker image.
 + Added configuration defaults, validation-only mode, reference checks and run metadata.
 + Fixed alignment offsets/CIGAR handling, SAM input, missing coverage, TTS boundary searches and annotation coordinates.
 + Fixed expression/merged measurements and added opt-in local enrichment with BY/BH correction.
 + Updated documentation, templates, helper scripts and synthetic example; added regression and workflow tests.

### version 2.0.0 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 27.11.2025
 + Reworked/Modernized the code
 + Reworked output structure
 + Added unit tests + fixed small bugs
 + Added toy example
 + Added License
 + updated README

### version 1.2.0 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 21.02.2021

 + created bam input version of ORFBounder
 + updated input file parsing
 + updated parameter parsing
 + updated README

### version 1.1.0 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 11.11.2021

 + full restructering of ORFBounder code
 + removal of outside tools
 + added new script to perform multiple ORFBounder runs

### version 1.0.0 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 01.12.2020

 + initial commit
 + renamed to ORFBounder
