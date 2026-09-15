"""Write a self-contained, offline view of the exported ORF candidate table."""

from collections import Counter
from html import escape
import json
import math
from numbers import Integral, Real
from pathlib import Path
import re
from urllib.parse import urlsplit

import pandas as pd

from lib.matched_statistics import ATTACH_FIELDS


_SAMPLE_FIELDS = (
    "peak_height", "relative_density", "rpkm", "peak_count",
    "background_count", "peak_fold_enrichment", "peak_p_value",
    "peak_q_value", "peak_significant", "peak_correction",
    "peak_tests_in_family", "peak_fit_status",
)

_MATCHED_GENERIC_FIELDS = {"p_value", "q_value", "significant", "correction", "tests_in_family"}
_MATCHED_PREFIX_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
_JAVASCRIPT_SAFE_INTEGER = (1 << 53) - 1


def _json_value(value):
    """Retain missing values as null and serialize pandas/numpy scalar values."""
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if value is None or value is pd.NA or value is pd.NaT:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bool):
        return value
    if isinstance(value, Integral):
        value = int(value)
        # JSON itself supports arbitrary decimal integers, but the report is
        # consumed by JavaScript, whose Number type would silently round these
        # values. Decimal text preserves exact audit counts and provenance IDs.
        return value if abs(value) <= _JAVASCRIPT_SAFE_INTEGER else str(value)
    if isinstance(value, Real):
        if math.isnan(value):
            return None
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return float(value)
    return value if isinstance(value, str) else str(value)


def _positive(value):
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value > 0
    if isinstance(value, str) and re.fullmatch(r"[+-]?\d+", value):
        return int(value) > 0
    return False


def _true(value):
    # Also accept values read back from TSVs, without treating the text "False"
    # or a missing value as a true statistical support flag.
    return value is True or value == 1 or (isinstance(value, str) and value.lower() == "true")


def _local_link(target):
    target = str(target)
    parsed = urlsplit(target)
    if (not target or target != target.lstrip() or parsed.scheme or parsed.netloc
            or parsed.path.startswith("/") or target.startswith(("/", "\\"))
            or "\\" in target or any(ord(char) < 32 or ord(char) == 127 for char in target)):
        raise ValueError("Report links must be relative paths to local result files.")
    return escape(target, quote=True)


class _ScriptSafeJSONWriter:
    """Escape JSON chunks while writing them inside an HTML script element."""

    def __init__(self, output):
        self.output = output

    def write(self, chunk):
        # JSON script elements are still tokenized as HTML. Escaping tag
        # delimiters prevents values from terminating the element; escaping
        # ampersands also avoids character-reference surprises in consumers.
        chunk = chunk.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
        chunk = chunk.replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
        return self.output.write(chunk)


class _BufferedScriptSafeJSONWriter:
    """Buffer encoder fragments without retaining the complete JSON document."""

    def __init__(self, output, buffer_size=1024 * 1024):
        self.writer = _ScriptSafeJSONWriter(output)
        self.buffer_size = buffer_size
        self.chunks = []
        self.buffered_characters = 0

    def write(self, chunk):
        self.chunks.append(chunk)
        self.buffered_characters += len(chunk)
        if self.buffered_characters >= self.buffer_size:
            self.flush()

    def flush(self):
        if self.chunks:
            self.writer.write("".join(self.chunks))
            self.chunks = []
            self.buffered_characters = 0


def _json_encoder():
    return json.JSONEncoder(
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    )


def _write_json_value(value, encoder, writer):
    for chunk in encoder.iterencode(value):
        writer.write(chunk)


def _write_bounded_json_value(value, encoder, writer):
    """Use the C encoder for a value whose size is already explicitly bounded."""
    writer.write(encoder.encode(value))


def _write_script_safe_json(value, output):
    """Stream compact, standards-compliant JSON without retaining a second copy."""
    writer = _BufferedScriptSafeJSONWriter(output)
    _write_json_value(value, _json_encoder(), writer)
    writer.flush()


def _infer_matched_comparisons(columns):
    """Recognize attached comparison columns without confusing local samples.

    Generic suffixes such as ``q_value`` are not sufficient: a single-sample
    local field also ends in that text. At least one matched-specific field must
    establish the prefix before the remaining attached fields are associated
    with it.
    """
    signature_fields = set(ATTACH_FIELDS) - _MATCHED_GENERIC_FIELDS
    prefixes = []
    for column in columns:
        # Local sample fields can share a generic matched suffix (notably
        # ``fit_status``); they must never manufacture a comparison prefix.
        if any(column.endswith("_" + field) for field in _SAMPLE_FIELDS):
            continue
        for field in sorted(signature_fields, key=len, reverse=True):
            suffix = "_" + field
            if not column.endswith(suffix):
                continue
            prefix = column[:-len(suffix)]
            if _MATCHED_PREFIX_PATTERN.fullmatch(prefix) and prefix not in prefixes:
                prefixes.append(prefix)
            break

    result = []
    for prefix in prefixes:
        match = re.fullmatch(r"(.+)_(TIS|TTS)", prefix)
        comparison = match.group(1) if match else prefix
        assay = match.group(2) if match else "Unknown"
        result.append({
            "prefix": prefix,
            "comparison": comparison,
            "assay": assay,
            "numerator_condition": None,
            "denominator_condition": None,
        })
    return result


def _matched_comparisons(columns, metadata):
    """Deduplicate declared comparison metadata and add column-only fallbacks."""
    def present(value):
        return value is not None and value != ""

    comparisons = {}
    for run in metadata:
        declared = run.get("matched_comparisons") or []
        if not isinstance(declared, list):
            continue
        for item in declared:
            if not isinstance(item, dict):
                continue
            name = item.get("comparison")
            assay = item.get("assay")
            prefix = item.get("prefix")
            if prefix is None and present(name) and present(assay):
                prefix = f"{name}_{str(assay).upper()}"
            if not isinstance(prefix, str) or _MATCHED_PREFIX_PATTERN.fullmatch(prefix) is None:
                continue
            inferred = re.fullmatch(r"(.+)_(TIS|TTS)", prefix)
            record = {
                "prefix": prefix,
                "comparison": str(name) if present(name) else inferred.group(1) if inferred else prefix,
                "assay": str(assay).upper() if present(assay) else inferred.group(2) if inferred else "Unknown",
                "numerator_condition": item.get("numerator_condition") or None,
                "denominator_condition": item.get("denominator_condition") or None,
            }
            for field in (
                    "fdr", "fdr_method", "model", "numeric_backend", "alternative", "replicate_ids",
                    "minimum_informative_pairs", "support_rule",
                    "max_exact_convolution_work", "maximum_supported_hypergeometric_total",
                    "minimum_positive_p_value",
                    "candidate_table", "strata_table"):
                if field in item:
                    record[field] = item[field]
            if prefix not in comparisons:
                comparisons[prefix] = record
            else:
                # Sample manifests repeat experiment-level metadata. Fill
                # absent values, but never let manifest order silently choose
                # between contradictory labels or model settings.
                for field, value in record.items():
                    existing = comparisons[prefix].get(field)
                    if present(existing) and present(value) and existing != value:
                        raise ValueError(
                            f"Conflicting matched-comparison metadata for {prefix!r}, "
                            f"field {field!r}."
                        )
                    if not present(existing) and present(value):
                        comparisons[prefix][field] = value

    for inferred in _infer_matched_comparisons(columns):
        comparisons.setdefault(inferred["prefix"], inferred)
    return list(comparisons.values())


def _prepare_schema(result_df, run_metadata):
    """Prepare the small report metadata shared by every candidate row."""
    columns = [str(column) for column in result_df.columns]
    if len(set(columns)) != len(columns):
        raise ValueError("Report table must have unique column names.")
    metadata = _json_value(run_metadata or [])
    matched_comparisons = _matched_comparisons(columns, metadata)
    matched_columns = {
        f"{comparison['prefix']}_{field}"
        for comparison in matched_comparisons
        for field in ATTACH_FIELDS
    }
    assays = {}
    for run in metadata:
        for assay, filename in (run.get("inputs") or {}).items():
            if filename and assay in {"TIS", "TTS", "RIBO"}:
                assays[Path(filename).stem.split("_")[0]] = assay
        # Canonical input paths can differ from the symlink names used to label
        # the exported measurements. Explicit run labels take precedence.
        for assay, sample in (run.get("sample_names") or {}).items():
            if sample and assay in {"TIS", "TTS", "RIBO"}:
                assays[str(sample)] = assay
    sample_names = []
    for column in columns:
        if column in matched_columns:
            continue
        for field in _SAMPLE_FIELDS:
            if column.endswith("_" + field):
                sample = column[:-(len(field) + 1)]
                if sample not in sample_names:
                    sample_names.append(sample)
                break
    samples = []
    for sample in sample_names:
        if (sample in {"TIS", "TTS", "RIBO"} and sample not in assays
                and not any(result_df[column].notna().any() for column in result_df.columns
                            if str(column).startswith(sample + "_"))):
            continue
        prefix = sample.split("-", 1)[0].upper()
        assay = assays.get(sample, prefix if prefix in {"TIS", "TTS", "RIBO", "RNA", "RNATIS", "RNATTS"} else "Other")
        samples.append({"name": sample, "assay": assay})
    return {
        "columns": columns,
        "samples": samples,
        "matched_comparisons": matched_comparisons,
        "runs": metadata,
    }


def _classify_row(get_value, samples, matched_comparisons):
    """Return the compact evidence indexes used by the browser interface."""
    called, significant, tested = [], [], []
    for sample in samples:
        name = sample["name"]
        if sample["assay"] not in {"TIS", "TTS"}:
            continue
        if _positive(get_value(name + "_peak_height")):
            called.append(name)
        if any(get_value(name + "_" + field) is not None for field in (
                "peak_significant", "peak_p_value", "peak_q_value")):
            tested.append(name)
        if _true(get_value(name + "_peak_significant")):
            significant.append(name)

    matched_tested, matched_significant = [], []
    for comparison in matched_comparisons:
        prefix = comparison["prefix"]
        status = get_value(prefix + "_fit_status")
        if (any(get_value(prefix + "_" + field) is not None
                for field in ("p_value", "q_value", "significant"))
                or status == "ok"):
            matched_tested.append(prefix)
        if _true(get_value(prefix + "_significant")):
            matched_significant.append(prefix)
    return called, significant, tested, matched_tested, matched_significant


def _prepared_row(raw_values, index, schema):
    values = [_json_value(value) for value in raw_values]

    def get_value(column):
        position = index.get(column)
        return values[position] if position is not None else None

    evidence = _classify_row(
        get_value, schema["samples"], schema["matched_comparisons"],
    )
    return dict(zip(
        ("values", "called", "significant", "tested", "matched_tested", "matched_significant"),
        (values, *evidence),
    ))


def _iter_prepared_rows(result_df, schema):
    index = {column: position for position, column in enumerate(schema["columns"])}
    for values in result_df.itertuples(index=False, name=None):
        yield _prepared_row(values, index, schema)


def _prepare_data(result_df, run_metadata):
    """Build the complete report value for compatibility and focused tests."""
    schema = _prepare_schema(result_df, run_metadata)
    return {
        "columns": schema["columns"],
        "rows": list(_iter_prepared_rows(result_df, schema)),
        "samples": schema["samples"],
        "matched_comparisons": schema["matched_comparisons"],
        "runs": schema["runs"],
    }


def _summarize_report(result_df, schema):
    """Calculate report cards without copying all table values into Python lists."""
    index = {column: position for position, column in enumerate(schema["columns"])}
    category_index = index.get("Type")
    assays = {sample["name"]: sample["assay"] for sample in schema["samples"]}
    categories = Counter()
    tis_count = tts_count = local_count = matched_count = 0

    for values in result_df.itertuples(index=False, name=None):
        def get_value(column):
            position = index.get(column)
            return _json_value(values[position]) if position is not None else None

        evidence = _classify_row(
            get_value, schema["samples"], schema["matched_comparisons"],
        )
        called, significant, _, _, matched_significant = evidence
        category = _json_value(values[category_index]) if category_index is not None else None
        categories[str(category or "Unspecified")] += 1
        tis_count += any(assays[sample] == "TIS" for sample in called)
        tts_count += any(assays[sample] == "TTS" for sample in called)
        local_count += bool(significant)
        matched_count += bool(matched_significant)
    return categories, tis_count, tts_count, local_count, matched_count


def _write_report_data(result_df, schema, output):
    """Serialize rows incrementally so memory is bounded by one candidate row."""
    writer = _BufferedScriptSafeJSONWriter(output)
    encoder = _json_encoder()
    writer.write('{"columns":')
    _write_bounded_json_value(schema["columns"], encoder, writer)
    writer.write(',"rows":[')
    for row_number, row in enumerate(_iter_prepared_rows(result_df, schema)):
        if row_number:
            writer.write(",")
        _write_bounded_json_value(row, encoder, writer)
    writer.write('],"samples":')
    _write_bounded_json_value(schema["samples"], encoder, writer)
    writer.write(',"matched_comparisons":')
    _write_bounded_json_value(schema["matched_comparisons"], encoder, writer)
    writer.write(',"runs":')
    _write_bounded_json_value(schema["runs"], encoder, writer)
    writer.write("}")
    writer.flush()


def write_report(
    result_df: pd.DataFrame,
    output_path: Path,
    title: str,
    *,
    run_metadata: list[dict] | None = None,
    links: dict[str, str] | None = None,
) -> Path:
    """Write an offline HTML report without changing the table or its statistics.

    ``links`` maps user-facing labels to paths relative to this HTML file.
    Called support uses positive TIS/TTS peak heights; optional local and
    fixed-margin matched support use their saved flags, including at inferred
    ORF boundaries. No p-values or support scores are computed by the report.
    """
    schema = _prepare_schema(result_df, run_metadata)
    categories, tis_count, tts_count, local_count, matched_count = _summarize_report(
        result_df, schema,
    )
    statistics_enabled = any(run.get("statistics") == "local" for run in schema["runs"])
    statistics_present = any(column.endswith("_peak_q_value") or column.endswith("_peak_significant")
                             for column in schema["columns"]) or statistics_enabled
    card_values = [
        ("ORF candidates", len(result_df), "Rows in this result table"),
        ("With a called TIS peak", tis_count, "At least one TIS sample"),
        ("With a called TTS peak", tts_count, "At least one TTS sample"),
        ("With local statistical support", local_count, "At least one TIS/TTS boundary" if statistics_present
         else "Local statistics unavailable"),
    ]
    if schema["matched_comparisons"]:
        card_values.append((
            "With exploratory fixed-margin matched support",
            matched_count,
            "Configured q/effect/direction rule passed with at least two informative pairs",
        ))
    cards = "".join(
        f'<div class="stat"><span>{escape(label)}</span><strong>{value:,}</strong><small>{escape(note)}</small></div>'
        for label, value, note in card_values
    )
    category_html = "".join(
        f'<span class="category">{escape(name.replace("_", " "))} <b>{count:,}</b></span>'
        for name, count in sorted(categories.items())
    ) or '<span class="muted">No candidate categories to display.</span>'
    link_html = "".join(
        f'<a class="file-link" href="{_local_link(target)}">{escape(str(label))}<span aria-hidden="true"> ↗</span></a>'
        for label, target in (links or {}).items()
    )
    files = (f'<section class="files" aria-labelledby="files-heading"><h2 id="files-heading">Open result files</h2>'
             f'<p>Keep this report alongside the result folders so these links work.</p><div class="file-links">'
             f'{link_html}</div></section>') if link_html else ""
    notice = ("Local support means the original sample's boundary test met its configured criteria. "
              "Each sample retains its own multiple-testing correction. It can support an inferred boundary "
              "even without a called peak. This is exploratory local enrichment, not a combined test or "
              "a probability that an ORF is translated.") if statistics_present else (
              "Local statistics were not included in this table. Called peaks are based on the run's "
              "peak-height settings; they do not imply statistical significance.")
    if statistics_enabled and result_df.empty:
        notice = ("Local statistics were enabled, but this result contains no ORF candidates. "
                  "Use the statistical candidate TSV files to inspect the codon-window tests. " + notice)
    if schema["matched_comparisons"]:
        notice += (
            " Exploratory fixed-margin matched support compares local peak-versus-background "
            "endpoint odds within predeclared independent experimental blocks. Its conditional "
            "p/q-values are read-depth evidence. Support requires a q-value at or below the configured "
            "FDR, a common odds ratio above 1, at least two informative pairs and a strict majority "
            "favoring the numerator. It does not estimate biological "
            "dispersion or a population-level condition effect, and it is not the single-sample "
            "local test, differential abundance, or evidence that the ORF is translated."
        )
    replacements = {"TITLE": escape(str(title)), "STYLES": _CSS, "CARDS": cards,
                    "CATEGORIES": category_html, "NOTICE": escape(notice), "FILES": files,
                    "SCRIPT": _JS}
    # Replace only markers in the template, never ones supplied in user values.
    before_data, after_data = _HTML.split("__DATA__", 1)
    marker_pattern = r"__(TITLE|STYLES|CARDS|CATEGORIES|NOTICE|FILES|SCRIPT)__"
    before_data = re.sub(marker_pattern, lambda match: replacements[match.group(1)], before_data)
    after_data = re.sub(marker_pattern, lambda match: replacements[match.group(1)], after_data)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as output:
        output.write(before_data)
        _write_report_data(result_df, schema, output)
        output.write(after_data)
    return output_path


_HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ · ORFBounder</title><style>__STYLES__</style></head>
<body><a class="skip-link" href="#candidates">Skip to candidates</a><main>
<header><div class="brand">ORFBounder <span>RESULTS</span></div><h1>__TITLE__</h1>
<p class="intro">Explore candidate open reading frames, compare sample evidence, and take sequences into your next analysis.</p>
<p class="offline">Standalone report · Works offline · Filters only change this view</p></header>
<section aria-label="Candidate overview"><div class="stats">__CARDS__</div>
<div class="categories" aria-label="Candidate categories">__CATEGORIES__</div></section>
<aside class="interpretation"><strong>Read each evidence layer separately</strong>
<p>A called TIS or TTS peak records support from that assay. The opposite boundary may have been inferred from the sequence.
Support counts can overlap when the same ORF has evidence from several samples.</p><p>__NOTICE__</p></aside>
<section id="candidates" aria-labelledby="candidates-heading"><div class="section-title"><h2 id="candidates-heading">Explore candidates</h2>
<p>Coordinates are 1-based and inclusive. On the minus strand, translation runs from Stop toward Start. A circular origin crossing uses a virtual Stop greater than Reference_length.</p></div>
<noscript><p class="empty">Enable JavaScript in your browser to explore candidates. The overview above and linked result files remain available.</p></noscript>
<div id="controls" class="controls" hidden>
<label class="search">Search candidates<input id="search" type="search" placeholder="Locus, identifier, sequence or coordinate" autocomplete="off"></label>
<label>Category<select id="category"><option value="">Any category</option></select></label>
<label>Called support<select id="called"><option value="">Any called support</option><option value="TIS">TIS called</option>
<option value="TTS">TTS called</option><option value="both">TIS and TTS called</option></select></label>
<label>Called sample<select id="sample"><option value="">Any sample</option></select></label>
<label>Local statistics<select id="local"><option value="">Any local evidence</option><option value="supported">At least one supported boundary</option>
<option value="unsupported">Tested without local support</option><option value="unavailable">Local support unavailable</option></select></label>
<label>Sort by<select id="sort"><option value="coordinates">Genome and coordinates</option><option value="support">Most called samples</option>
<option value="length-desc">Longest ORF first</option><option value="length-asc">Shortest ORF first</option></select></label>
<label>Rows per page<select id="page-size"><option>25</option><option>50</option><option>100</option></select></label>
<button id="reset" class="secondary" type="button">Reset filters</button></div>
<p id="table-status" class="table-status" role="status" aria-live="polite"></p>
<div class="table-scroll"><table class="candidate-table"><caption class="sr-only">ORF candidates matching the selected filters; open Details for sample evidence and sequences.</caption>
<thead><tr><th scope="col">Candidate</th><th scope="col">Position</th><th scope="col">Category</th><th scope="col">Codons</th>
<th scope="col">Called samples</th><th scope="col">Local support</th><th scope="col">Exploratory matched support</th><th scope="col">Inspect</th></tr></thead><tbody id="candidate-rows"></tbody></table></div>
<p id="empty-state" class="empty" hidden></p>
<nav id="pagination" class="pagination" aria-label="Candidate pages" hidden><button id="previous" class="secondary" type="button">Previous</button>
<span id="page-label"></span><button id="next" class="secondary" type="button">Next</button></nav></section>
__FILES__
<section class="guide" aria-labelledby="reading-heading"><h2 id="reading-heading">Use these results</h2>
<div class="guide-grid"><div><h3>Inspect a candidate</h3><p>Open Details to compare called peaks, single-sample local effect sizes and tests, exploratory fixed-margin matched support and expression measurements.
An unavailable measurement is shown as “—”; it is not a zero or a negative result.</p></div>
<div><h3>Check the genomic context</h3><p>Load the GFF and coverage tracks with the same reference genome in a genome browser.
Review the strand, peak position and surrounding genes before selecting candidates for follow-up.</p></div>
<div><h3>Reuse the sequences</h3><p>Copy sequences from candidate details or use the FASTA exports. Nucleotide and protein sequences follow the predicted ORF's strand.
The codon count includes the stop codon; a terminal * in a protein denotes a stop.</p></div></div>
<details><summary>What the measurements mean</summary><dl class="glossary">
<dt>Category</dt><dd>How the predicted coordinates relate to the supplied annotation. It does not grade biological confidence.</dd>
<dt>Circular coordinates</dt><dd>Is_circular records an explicit whole-reference GFF3 declaration. An origin-spanning table interval has Stop greater than Reference_length; its GFF uses physical multipart rows and its WIG positions remain within 1..Reference_length.</dd>
<dt>Peak height</dt><dd>The run's maximum or summed codon-window signal after the selected normalization. Check Run settings before comparing values across runs.</dd>
<dt>Local fold enrichment</dt><dd>A descriptive ratio of endpoint density in the peak window to density in its adjacent background. Infinity means a positive peak count with zero background endpoints and no pseudocount; unavailable is not zero.</dd>
<dt>Local support and q-value</dt><dd>Support uses the saved sample-level flag. A q-value is the adjusted p-value from that sample's eligible codon tests; it is not a combined score across samples. A correction limit retains the valid raw p-value but uses q=1 and cannot provide support.</dd>
<dt>Exploratory fixed-margin matched support</dt><dd>A read-depth endpoint comparison within real, predeclared independent experimental blocks. Support requires a q-value at or below the configured FDR, a common odds ratio above 1, at least two informative pairs and strict-majority directional agreement. Separate comparisons and assays are corrected separately; the common odds ratio can conceal heterogeneous pair effects. Join the aggregate and companion strata tables by candidate_id to inspect the per-pair counts and directions. This is not a biological-dispersion or population-effect test.</dd>
<dt>Relative density, RPKM and TE</dt><dd>Descriptive peak-to-gene density, expression scaled by ORF length and the selected contig or library read total, and matched assay/RNA expression ratios. RPKM is displayed to two decimal places; TE uses the unrounded rates, so a very low positive displayed RPKM can still contribute to TE. These measurements do not constitute a statistical comparison.</dd>
</dl></details></section>
<section id="run-section" class="runs" aria-labelledby="run-heading" hidden><h2 id="run-heading">Run settings</h2>
<p>Use these records to check how reads, peaks and sequences were processed.</p><div id="run-records"></div></section>
<footer>Generated by ORFBounder. Candidate evidence supports review and follow-up; annotation categories and support counts do not establish translation.</footer>
</main><script id="report-data" type="application/json">__DATA__</script><script>__SCRIPT__</script></body></html>
"""


_CSS = """
:root{color-scheme:light;--ink:#153441;--muted:#516671;--teal:#086b65;--line:#dbe5e6;--paper:#fff;--wash:#f4f8f8}
*{box-sizing:border-box}body{margin:0;background:var(--wash);color:var(--ink);font:15px/1.6 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
main{max-width:1440px;margin:auto;padding:44px 38px 28px}header{margin-bottom:30px}.brand{font-weight:750;font-size:20px;letter-spacing:-.5px}.brand span{font-size:10px;letter-spacing:2px;margin-left:12px;color:var(--teal);vertical-align:middle}
h1{font-size:clamp(29px,4vw,44px);line-height:1.15;font-weight:650;letter-spacing:-1.4px;margin:23px 0 14px;overflow-wrap:anywhere}h2{font-size:22px;line-height:1.3;letter-spacing:-.4px;margin:0 0 8px}h3{font-size:16px;margin:0 0 6px}p{margin:8px 0}.intro{max-width:820px;font-size:17px;color:var(--muted)}.offline{font-size:12px;color:var(--muted);margin-top:14px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:14px}.stat{background:var(--paper);border:1px solid var(--line);border-radius:10px;padding:20px 22px;display:flex;flex-direction:column}.stat span{font-size:13px;font-weight:650}.stat strong{font-size:35px;font-weight:650;letter-spacing:-1px;line-height:1.6}.stat small{font-size:12px;color:var(--muted)}.categories{display:flex;gap:8px;flex-wrap:wrap;margin:15px 0 22px}.category{padding:4px 11px;border-radius:5px;background:#e6eeee;font-size:12px}.category b{margin-left:8px}
.interpretation{border-left:3px solid var(--teal);padding:15px 20px;background:#e9f3f1;margin:0 0 36px;border-radius:0 8px 8px 0;font-size:13px}.interpretation strong{font-size:14px}.interpretation p{max-width:1150px}
.section-title p,.files>p,.runs>p{color:var(--muted);font-size:13px}.controls{display:flex;gap:14px 12px;flex-wrap:wrap;align-items:end;margin-top:19px;padding:19px;background:var(--paper);border:1px solid var(--line);border-radius:9px}.controls label{flex:1 1 155px;display:flex;flex-direction:column;gap:5px;font-size:12px;font-weight:650}.controls .search{flex-basis:310px}
.controls label:last-of-type{flex:0 1 130px}input,select,button{font:inherit}input,select{width:100%;min-height:40px;padding:8px 10px;border:1px solid #aebfc4;border-radius:5px;background:#fff;color:var(--ink);font-size:13px}button{cursor:pointer;border:1px solid var(--teal);background:var(--teal);color:#fff;border-radius:5px;padding:8px 12px;font-size:13px;font-weight:600;min-height:38px}.secondary{background:var(--paper);border-color:#aebfc4;color:var(--ink)}button:hover{filter:brightness(.94)}button:disabled{opacity:.45;cursor:default}a{color:var(--teal);text-underline-offset:3px}:focus-visible{outline:3px solid #d27819;outline-offset:3px}
.table-status{font-size:13px;color:var(--muted);margin:19px 0 9px}.table-scroll{overflow-x:auto;border:1px solid var(--line);border-radius:8px;background:var(--paper)}table{border-collapse:collapse;width:100%;text-align:left}th{font-size:11px;letter-spacing:.4px;color:var(--muted);font-weight:700;background:#eaf0f0}th,td{padding:13px 15px;border-bottom:1px solid var(--line);vertical-align:top}td{font-size:13px}.candidate-table{min-width:910px}.candidate-table>tbody>tr:last-child>td{border-bottom:0}.candidate-name{font-weight:650}.subtext{display:block;color:var(--muted);font-size:11px;overflow-wrap:anywhere}.candidate-table .position{white-space:nowrap}.pill{display:inline-block;border-radius:4px;padding:2px 6px;background:#edf2f3;font-size:11px;margin:1px 4px 2px 0;white-space:nowrap}.pill.supported{background:#def0e8;color:#155341}.muted{color:var(--muted)}
.pagination{display:flex;justify-content:flex-end;align-items:center;gap:14px;margin:14px 0 32px;font-size:13px}.empty{padding:25px;border:1px dashed #aebfc4;border-radius:8px;background:var(--paper);font-size:14px}.detail-row>td{padding:0 20px 24px;background:#f8fbfa}.candidate-detail{border-top:2px solid var(--teal);padding-top:20px}.candidate-detail>p{font-size:13px}.candidate-detail h3{margin-top:18px}.sample-table{min-width:760px;background:#fff}.matched-table{min-width:1320px;background:#fff}.sample-table th,.sample-table td,.matched-table th,.matched-table td{padding:9px 11px}.sample-table th,.matched-table th{font-size:11px}
details{margin-top:16px}summary{cursor:pointer;font-weight:650;font-size:13px;color:var(--teal);padding:4px 0}.metrics{width:100%;max-width:1100px;table-layout:fixed}.metrics th{width:40%;background:transparent;overflow-wrap:anywhere;text-transform:none;letter-spacing:0;font-size:12px}.metrics td{overflow-wrap:anywhere;font-variant-numeric:tabular-nums}.sequence-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:15px}.sequence-box:last-child{grid-column:1/-1}.sequence-box label{display:block;font-size:12px;font-weight:650;margin:12px 0 6px}.sequence-box textarea{display:block;resize:vertical;width:100%;min-height:72px;padding:10px;border:1px solid #b7c9c8;border-radius:5px;font:12px/1.7 ui-monospace,SFMono-Regular,Consolas,monospace;color:var(--ink);background:#fff;word-break:break-all}
.files,.guide,.runs{margin-top:36px;border-top:1px solid var(--line);padding-top:26px}.file-links{display:flex;gap:10px;flex-wrap:wrap;margin-top:14px}.file-link{padding:8px 12px;background:#fff;border:1px solid var(--line);border-radius:5px;text-decoration:none;font-size:13px}.guide-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:28px;margin-top:19px}.guide p{color:var(--muted);font-size:13px}.glossary{font-size:13px;max-width:1100px}.glossary dt{font-weight:650;margin-top:12px}.glossary dd{margin:3px 0;color:var(--muted)}.run-card{padding:16px 20px;border:1px solid var(--line);background:#fff;border-radius:8px;margin:12px 0}.run-card h3{overflow-wrap:anywhere}.run-fields{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px 18px;font-size:12px}.run-fields dt{color:var(--muted)}.run-fields dd{margin:0;overflow-wrap:anywhere}.run-card pre{max-height:340px;overflow:auto;padding:12px;background:var(--wash);font-size:11px;white-space:pre-wrap;overflow-wrap:anywhere}footer{font-size:11px;color:var(--muted);margin-top:40px;padding-top:18px;border-top:1px solid var(--line)}
.sr-only,.skip-link:not(:focus){position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}.skip-link:focus{display:block;padding:12px;background:#fff}[hidden]{display:none!important}
@media(max-width:900px){main{padding:28px 20px}.stats{grid-template-columns:repeat(2,1fr)}.guide-grid{grid-template-columns:1fr;gap:16px}.run-fields{grid-template-columns:repeat(2,minmax(0,1fr))}}
@media(max-width:520px){main{padding:24px 12px}.stat{padding:15px}.stat span{font-size:11px}.stat strong{font-size:29px}.sequence-grid{grid-template-columns:1fr}.controls{padding:12px}.controls label{flex-basis:100%}}
@media print{body{background:#fff}main{padding:0}.controls,.pagination,button,.offline{display:none!important}.stats{grid-template-columns:repeat(auto-fit,minmax(150px,1fr))}.table-scroll{overflow:visible}.candidate-table{min-width:0}th,td{padding:7px}.interpretation{break-inside:avoid}}
"""


_JS = r"""
"use strict";
(() => {
  const data = JSON.parse(document.getElementById("report-data").textContent);
  const columns = new Map(data.columns.map((name, index) => [name, index]));
  const assays = new Map(data.samples.map(sample => [sample.name, sample.assay]));
  const matchedComparisons = new Map(data.matched_comparisons.map(comparison => [comparison.prefix, comparison]));
  const get = (row, key) => columns.has(key) ? row.values[columns.get(key)] : null;
  const text = value => value === null || value === undefined || value === "" ? "—" : String(value);
  const pretty = value => text(value).replaceAll("_", " ");
  const comparisonLabel = prefix => {
    const comparison = matchedComparisons.get(prefix);
    return comparison ? text(comparison.comparison) + " (" + text(comparison.assay) + ")" : prefix;
  };
  const number = value => typeof value === "number" ? (Number.isInteger(value) ? value.toLocaleString() :
    Math.abs(value) > 0 && Math.abs(value) < 0.001 ? value.toExponential(3) : String(Number(value.toPrecision(5)))) :
    typeof value === "string" && /^[+-]?\d+$/.test(value) ? BigInt(value).toLocaleString() : text(value);
  const el = (tag, value, className) => {
    const node = document.createElement(tag);
    if (value !== undefined) node.textContent = text(value);
    if (className) node.className = className;
    return node;
  };
  const ids = ["search", "category", "called", "sample", "local", "sort", "page-size"];
  const controls = Object.fromEntries(ids.map(id => [id, document.getElementById(id)]));
  let page = 0;
  let filtered = [];
  // Search values in place instead of retaining a second, lower-cased copy of
  // the full result table. Long nucleotide/protein columns otherwise roughly
  // double the browser memory needed merely to open a report.
  const includesQuery = (row, query) => row.values.some(value =>
    value !== null && String(value).toLocaleLowerCase().includes(query));
  const addOption = (select, value, label) => {
    const option = el("option", label);
    option.value = value;
    select.append(option);
  };
  [...new Set(data.rows.map(row => text(get(row, "Type"))))].sort().forEach(category =>
    addOption(controls.category, category, pretty(category)));
  data.samples.filter(sample => ["TIS", "TTS"].includes(sample.assay)).forEach(sample =>
    addOption(controls.sample, sample.name, sample.name));
  const hasAssay = (row, assay) => row.called.some(name => assays.get(name) === assay);
  const numeric = (row, column) => {
    const value = get(row, column);
    return typeof value === "number" ? value : Number.POSITIVE_INFINITY;
  };
  const crossesOrigin = row => get(row, "Is_circular") === true &&
    typeof get(row, "Stop") === "number" && typeof get(row, "Reference_length") === "number" &&
    get(row, "Stop") > get(row, "Reference_length");
  const coordinates = (a, b) => text(get(a, "Genome")).localeCompare(text(get(b, "Genome")), undefined, {numeric:true}) ||
    numeric(a, "Start") - numeric(b, "Start") || numeric(a, "Stop") - numeric(b, "Stop") ||
    text(get(a, "Strand")).localeCompare(text(get(b, "Strand")));
  const matches = row => {
    const query = controls.search.value.trim().toLocaleLowerCase();
    if (query && !includesQuery(row, query)) return false;
    if (controls.category.value && text(get(row, "Type")) !== controls.category.value) return false;
    if (controls.sample.value && !row.called.includes(controls.sample.value)) return false;
    const called = controls.called.value;
    if (called === "both" && !(hasAssay(row, "TIS") && hasAssay(row, "TTS"))) return false;
    if (["TIS", "TTS"].includes(called) && !hasAssay(row, called)) return false;
    const local = controls.local.value;
    if (local === "supported" && !row.significant.length) return false;
    if (local === "unsupported" && (!row.tested.length || row.significant.length)) return false;
    if (local === "unavailable" && row.tested.length) return false;
    return true;
  };
  const addCell = (tr, value, className) => {
    const td = el("td", value, className);
    tr.append(td);
    return td;
  };
  const supportCell = (tr, names, empty, supported) => {
    const td = addCell(tr);
    if (!names.length) td.append(el("span", empty, "muted"));
    names.forEach(name => td.append(el("span", name, "pill" + (supported ? " supported" : ""))));
  };
  const detail = (row, detailId) => {
    const panel = el("div", undefined, "candidate-detail");
    panel.append(el("h3", text(get(row, "Identifier"))));
    panel.append(el("p", "Genome " + text(get(row, "Genome")) + " · " + number(get(row, "Start")) + "–" +
      number(get(row, "Stop")) + " · strand " + text(get(row, "Strand")) + " · " + number(get(row, "Codon_count")) + " codons" +
      (crossesOrigin(row) ? " · crosses circular origin (reference length " + number(get(row, "Reference_length")) + ")" :
       get(row, "Is_circular") === true ? " · circular reference" : "")));
    panel.append(el("h3", "Sample evidence"));
    panel.append(el("p", "No call recorded means this sample has no positive called peak in this table. An inferred boundary can still have local statistical support. Local fold enrichment is a descriptive endpoint-density ratio: Infinity means a positive peak with zero background endpoints and no pseudocount, zero means no peak endpoints with positive background, and unavailable means the ratio could not be calculated. A floored p-value remains positive but is below the reportable float range; numerical limit means p=1 was retained conservatively; correction limit means the valid raw p-value was retained but q=1 was used. Neither limit can provide support.", "muted"));
    const scroll = el("div", undefined, "table-scroll");
    const table = el("table", undefined, "sample-table");
    table.append(el("caption", "Called peak and local statistical evidence for this candidate", "sr-only"));
    const head = el("thead");
    const headings = el("tr");
    ["Sample", "Assay", "Called peak", "Peak height", "Local fold enrichment", "Local boundary support", "q-value", "Fit status", "RPKM"].forEach(label => {
      const th = el("th", label); th.scope = "col"; headings.append(th);
    });
    head.append(headings); table.append(head);
    const body = el("tbody");
    data.samples.forEach(sample => {
      const tr = el("tr");
      const name = sample.name;
      const isBoundaryAssay = ["TIS", "TTS"].includes(sample.assay);
      const flag = get(row, name + "_peak_significant");
      addCell(tr, name); addCell(tr, sample.assay);
      addCell(tr, isBoundaryAssay ? row.called.includes(name) ? "Called" : "No call recorded" : "Context measurement");
      addCell(tr, number(get(row, name + "_peak_height")));
      addCell(tr, isBoundaryAssay ? number(get(row, name + "_peak_fold_enrichment")) : "—");
      addCell(tr, isBoundaryAssay ? row.significant.includes(name) ? "Supported" : flag !== null ? "Not supported" : "Unavailable" : "—");
      addCell(tr, number(get(row, name + "_peak_q_value")));
      addCell(tr, isBoundaryAssay ? pretty(get(row, name + "_peak_fit_status")) : "—");
      addCell(tr, number(get(row, name + "_rpkm")));
      body.append(tr);
    });
    table.append(body); scroll.append(table); panel.append(scroll);
    if (data.matched_comparisons.length) {
      panel.append(el("h3", "Exploratory fixed-margin matched evidence"));
      panel.append(el("p", "This conditional endpoint-count test retains each block's margins. Its saved support flag requires a q-value at or below the configured FDR, a common odds ratio above 1, at least two informative pairs and a strict majority of informative pairs favoring the numerator. A correction limit retains the valid raw p-value but uses q=1 and cannot provide support. Read depth drives its p/q-values; it does not estimate biological dispersion or a population-level condition effect. Inspect pair heterogeneity in the companion _strata.tsv table, joined by candidate_id, and treat this separately from single-sample local support, abundance, and translation evidence.", "muted"));
      const matchedScroll = el("div", undefined, "table-scroll");
      const matchedTable = el("table", undefined, "matched-table");
      matchedTable.append(el("caption", "Exploratory fixed-margin matched local-enrichment evidence for this candidate", "sr-only"));
      const matchedHead = el("thead");
      const matchedHeadings = el("tr");
      ["Comparison", "Assay", "Conditions", "Matched pairs", "Informative pairs", "Numerator-favoring pairs",
        "Denominator-favoring pairs", "Tied pairs", "Common OR", "log2 OR", "q-value", "Exploratory support", "Fit status"].forEach(label => {
        const th = el("th", label); th.scope = "col"; matchedHeadings.append(th);
      });
      matchedHead.append(matchedHeadings); matchedTable.append(matchedHead);
      const matchedBody = el("tbody");
      data.matched_comparisons.forEach(comparison => {
        const tr = el("tr");
        const prefix = comparison.prefix;
        const conditions = comparison.numerator_condition !== null && comparison.denominator_condition !== null ?
          text(comparison.numerator_condition) + " > " + text(comparison.denominator_condition) : "Not recorded";
        const tested = row.matched_tested.includes(prefix);
        const supported = row.matched_significant.includes(prefix);
        addCell(tr, comparison.comparison); addCell(tr, comparison.assay); addCell(tr, conditions);
        addCell(tr, number(get(row, prefix + "_matched_pair_count")));
        addCell(tr, number(get(row, prefix + "_informative_pair_count")));
        addCell(tr, number(get(row, prefix + "_concordant_pair_count")));
        addCell(tr, number(get(row, prefix + "_discordant_pair_count")));
        addCell(tr, number(get(row, prefix + "_tied_pair_count")));
        addCell(tr, number(get(row, prefix + "_common_odds_ratio")));
        addCell(tr, number(get(row, prefix + "_log2_common_odds_ratio")));
        addCell(tr, number(get(row, prefix + "_q_value")));
        addCell(tr, supported ? "Supported" : tested ? "Not supported" : "Unavailable");
        addCell(tr, pretty(get(row, prefix + "_fit_status")));
        matchedBody.append(tr);
      });
      matchedTable.append(matchedBody); matchedScroll.append(matchedTable); panel.append(matchedScroll);
    }
    const sequenceFields = [["Nucleotide_Seq", "Nucleotide sequence"], ["Amino_Acid_Seq", "Protein sequence"], ["15nt_window", "15 nt sequence window"]];
    const sequences = el("div", undefined, "sequence-grid");
    sequenceFields.forEach(([column, label], index) => {
      if (!columns.has(column)) return;
      const box = el("div", undefined, "sequence-box");
      const name = el("label", label + " · select text to copy");
      const area = el("textarea");
      area.id = detailId + "-sequence-" + index; name.htmlFor = area.id;
      area.readOnly = true; area.spellcheck = false; area.rows = column === "15nt_window" ? 2 : 3;
      area.value = text(get(row, column)); box.append(name, area); sequences.append(box);
    });
    panel.append(sequences);
    const all = el("details");
    all.append(el("summary", "All table fields and exact values"));
    const metrics = el("table", undefined, "metrics");
    metrics.append(el("caption", "Original exported column names and unrounded values", "sr-only"));
    const metricBody = el("tbody");
    data.columns.forEach((column, index) => {
      const tr = el("tr"); const th = el("th", column); th.scope = "row";
      tr.append(th); addCell(tr, text(row.values[index])); metricBody.append(tr);
    });
    metrics.append(metricBody); all.append(metrics); panel.append(all);
    return panel;
  };
  const render = () => {
    const size = Number(controls["page-size"].value);
    const pages = Math.max(1, Math.ceil(filtered.length / size));
    page = Math.max(0, Math.min(page, pages - 1));
    const body = document.getElementById("candidate-rows");
    body.replaceChildren();
    filtered.slice(page * size, (page + 1) * size).forEach((row, index) => {
      const tr = el("tr");
      const identity = addCell(tr);
      identity.append(el("span", get(row, "Locus_tag") || get(row, "Identifier"), "candidate-name"),
        el("span", get(row, "Identifier"), "subtext"));
      const position = addCell(tr, undefined, "position");
      position.append(el("span", number(get(row, "Start")) + "–" + number(get(row, "Stop")) + " (" + text(get(row, "Strand")) + ")"),
        el("span", text(get(row, "Genome")) + (crossesOrigin(row) ? " · crosses origin" : ""), "subtext"));
      addCell(tr, pretty(get(row, "Type"))); addCell(tr, number(get(row, "Codon_count")));
      supportCell(tr, row.called, "No call recorded", false);
      supportCell(tr, row.significant, row.tested.length ? "No supported boundary" : "Unavailable", true);
      supportCell(tr, row.matched_significant.map(comparisonLabel),
        row.matched_tested.length ? "No exploratory matched support" : data.matched_comparisons.length ? "Unavailable" : "Not configured", true);
      const actions = addCell(tr);
      const button = el("button", "Details"); button.type = "button";
      const detailId = "candidate-detail-" + (page * size + index);
      button.setAttribute("aria-expanded", "false"); button.setAttribute("aria-controls", detailId);
      button.setAttribute("aria-label", "Details for " + text(get(row, "Identifier")));
      actions.append(button); body.append(tr);
      const detailRow = el("tr", undefined, "detail-row"); detailRow.id = detailId; detailRow.hidden = true;
      const cell = addCell(detailRow); cell.colSpan = 8; body.append(detailRow);
      button.addEventListener("click", () => {
        if (!cell.childNodes.length) cell.append(detail(row, detailId));
        detailRow.hidden = !detailRow.hidden; button.setAttribute("aria-expanded", String(!detailRow.hidden));
        button.textContent = detailRow.hidden ? "Details" : "Close";
      });
    });
    const start = filtered.length ? page * size + 1 : 0;
    const stop = Math.min((page + 1) * size, filtered.length);
    document.getElementById("table-status").textContent = "Showing " + start.toLocaleString() + "–" + stop.toLocaleString() +
      " of " + filtered.length.toLocaleString() + " matching candidates (" + data.rows.length.toLocaleString() + " total).";
    const empty = document.getElementById("empty-state");
    empty.hidden = filtered.length > 0;
    empty.textContent = data.rows.length ? "No candidates match these filters. Try another search or choose Reset filters." :
      "No ORF candidates were reported. Check the read counts and settings below, then inspect the coverage tracks and input reference. A zero-candidate result does not by itself establish an absence of translation.";
    document.getElementById("pagination").hidden = filtered.length === 0;
    document.getElementById("page-label").textContent = "Page " + (page + 1) + " of " + pages;
    document.getElementById("previous").disabled = page === 0;
    document.getElementById("next").disabled = page >= pages - 1;
  };
  const update = () => {
    filtered = data.rows.filter(matches);
    const sort = controls.sort.value;
    filtered.sort((a, b) => {
      if (sort === "support") return b.called.length - a.called.length || coordinates(a, b);
      if (sort.startsWith("length")) {
        const av = get(a, "Codon_count"), bv = get(b, "Codon_count");
        if (typeof av !== "number") return typeof bv === "number" ? 1 : coordinates(a, b);
        if (typeof bv !== "number") return -1;
        return (sort === "length-asc" ? av - bv : bv - av) || coordinates(a, b);
      }
      return coordinates(a, b);
    });
    page = 0; render();
  };
  ids.forEach(id => controls[id].addEventListener(id === "search" ? "input" : "change", update));
  document.getElementById("reset").addEventListener("click", () => {
    ids.forEach(id => { controls[id].value = id === "sort" ? "coordinates" : id === "page-size" ? "25" : ""; });
    update(); controls.search.focus();
  });
  document.getElementById("previous").addEventListener("click", () => { page -= 1; render(); });
  document.getElementById("next").addEventListener("click", () => { page += 1; render(); });
  const runFields = [
    ["mapping", "Read mapping"], ["normalization", "Peak normalization"], ["normalization_scope", "Normalization scope"],
    ["genetic_code", "Genetic code"], ["min_peak_height", "Minimum peak height"], ["peak_height_operator", "Peak operator"],
    ["start_codons", "Start codons"], ["stop_codons", "Stop codons"], ["tts_start_selection", "TTS start selection"],
    ["statistics", "Local statistics"], ["fdr", "FDR threshold"], ["fdr_method", "Correction"],
    ["background_width", "Background flank (nt)"], ["all_reads_rpkm", "Use all read lengths for RPKM"],
    ["completed_at", "Recorded at"], ["orfbounder_version", "ORFBounder version"]
  ];
  data.runs.forEach((run, index) => {
    const card = el("article", undefined, "run-card");
    card.append(el("h3", run.output_basename || "Run " + (index + 1)));
    const fields = el("dl", undefined, "run-fields");
    runFields.forEach(([key, label]) => {
      if (run[key] === undefined || run[key] === null) return;
      const group = el("div"); group.append(el("dt", label), el("dd", Array.isArray(run[key]) ? run[key].join(", ") : text(run[key])));
      fields.append(group);
    });
    card.append(fields);
    if (run.alignment_diagnostics) {
      const reads = el("details"); reads.append(el("summary", "Read counts by sample and contig"));
      Object.entries(run.alignment_diagnostics).forEach(([name, diagnostics]) => {
        reads.append(el("h3", name));
        const counts = diagnostics.accepted_reads_by_contig;
        if (counts) {
          const list = el("dl", undefined, "run-fields");
          Object.entries(counts).forEach(([contig, count]) => {
            const group = el("div"); group.append(el("dt", contig), el("dd", number(count) + " accepted reads")); list.append(group);
          });
          reads.append(list);
        }
        const all = el("details"); all.append(el("summary", "All recorded read counts"), el("pre", JSON.stringify(diagnostics, null, 2))); reads.append(all);
      });
      card.append(reads);
    }
    const all = el("details"); all.append(el("summary", "Full run record, input paths and checksums"), el("pre", JSON.stringify(run, null, 2)));
    card.append(all); document.getElementById("run-records").append(card);
  });
  document.getElementById("run-section").hidden = data.runs.length === 0;
  document.getElementById("controls").hidden = false;
  update();
})();
"""
