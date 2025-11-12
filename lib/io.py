#!/usr/bin/env python

from pathlib import Path
from typing import Optional


import csv
import json
import collections
import pandas as pd

from Bio import SeqIO


import lib.messaging as msg


def generate_genome_dict(genome_file: Path) -> dict[str, str]:
    """
    read a genome fasta file into a dictionary
    """
    genome_file_parsed = SeqIO.parse(genome_file, "fasta")
    genome_dict = {}
    for entry in genome_file_parsed:
        genome_dict[str(entry.id)] = str(entry.seq)

    return genome_dict


def parse_read_lengths(read_length_json: Path) -> Optional[dict[str, list[str]]]:
    """
    Parse the read length input into a continuous list form.
    """

    if not read_length_json.is_file():
        msg.warning(
            "Warning: Empty read-lengths parameter given, using all available read lengths."
        )
        return None

    with open(read_length_json, "r", encoding="utf-8") as json_file:
        tmp_dict = json.load(json_file)

    read_length_dict = {}
    for file, read_length_string in tmp_dict.items():
        if isinstance(read_length_string, int):
            read_length_dict[file] = [str(read_length_string)]

        elif isinstance(read_length_string, str):
            parts = read_length_string.split(",")
            read_lengths = set()
            for part in parts:
                if "-" in part:
                    interval = part.split("-")
                    i1, i2 = int(interval[0]), int(interval[1])
                    if i1 > i2:
                        i1, i2 = i2, i1

                    for i in range(i1, i2 + 1):
                        read_lengths.add(i)
                else:
                    read_lengths.add(int(part))

            read_length_dict[file] = [str(x) for x in sorted(read_lengths)]
        else:
            raise ValueError("Error: Read-length JSON file is not in correct format!")

    return read_length_dict


def parse_total_reads(
    mapped_counts_file_path: Path,
    normalization_method: str
) -> Optional[dict[str, int]]:
    """
    Takes a tab seperated file of total read counts and determines the minimum for each chromosome
    Format:  sample chromosome total_reads
    """

    min_read_count_dict = {}
    if normalization_method == "min":
        error_msg = (
            "Error: min normalization method chosen but no mapped_counts_file_path given!\n"
            "Either use a different normalization method or provide a file containing total read counts for each sample and each chromosome.\n"
            "Consider using our helper script to create the required files."
        )


        if not mapped_counts_file_path.is_file():
            raise FileNotFoundError(error_msg)

        with open(mapped_counts_file_path, "r", encoding="utf-8") as f:
            lines = list(filter(None, [line.strip() for line in f.readlines()]))

            for line in lines:
                _, chrom, cur_count = line.split("\t")

                if chrom in min_read_count_dict:
                    if min_read_count_dict[chrom] > int(cur_count):
                        min_read_count_dict[chrom] = int(cur_count)
                else:
                    min_read_count_dict[chrom] = int(cur_count)

    else:
        return None

    return min_read_count_dict


def parse_alignment_input(
    alignment_file_tis: Path | None,
    alignment_file_tts: Path | None
) -> str:
    """
    Check whether the input alignment files are valid and determine the execution method for ORFBounder
    """

    if alignment_file_tis and alignment_file_tts:
        if not alignment_file_tis.is_file():
            raise FileNotFoundError(
                f"Error: Non-empty alignment file path given for TIS does not exist: {alignment_file_tis}"
            )

        if not alignment_file_tts.is_file():
            raise FileNotFoundError(
                f"Error: Non-empty alignment file path given for TTS does not exist: {alignment_file_tts}"
            )

        return "combined_methods"

    if alignment_file_tis:
        if not alignment_file_tis.is_file():
            raise FileNotFoundError(
                f"Error: Non-empty alignment file path given for TIS does not exist: {alignment_file_tis}"
            )
        return "TIS"

    if alignment_file_tts:
        if not alignment_file_tts.is_file():
            raise FileNotFoundError(
                f"Error: Non-empty alignment file path given for TTS does not exist: {alignment_file_tts}"
            )
        return "TTS"

    raise FileNotFoundError("Error: Please ensure to either provide a TIS file, a TTS file or both!")


def check_alignment_path_input(
    alignment_file_path: Path | None,
    alignment_file_tis: Path | None,
    alignment_file_tts: Path | None
) -> Optional[set[Path]]:
    """
    Check alignment input path.
    Ensure that there is:
     - one sam/bam file corresponding to each input method (TIS, TTS)
     - (optional) one RNA bam file corresponding to each input method (TIS, TTS)
    """

    if not alignment_file_path or not alignment_file_path.is_dir():
        return None

    valid_bam = set()

    alignment_file_list = [
        file for file in alignment_file_path.iterdir()
        if file.suffix in {".bam", ".sam"}
    ]

    condition, replicate = "", ""
    if alignment_file_tis:
        tis_path = Path(alignment_file_tis)
        tis_prefix = tis_path.stem
        condition, replicate = tis_prefix.split("-")[1:]
        rnatis_prefix = f"RNATIS-{condition}-{replicate}"
        for file in alignment_file_list:
            if tis_prefix in file.name or rnatis_prefix in file.name:
                valid_bam.add(file)

    if alignment_file_tts:
        tts_path = Path(alignment_file_tts)
        tts_prefix = tts_path.stem
        condition, replicate = tts_prefix.split("-")[1:]
        rnatts_prefix = f"RNATTS-{condition}-{replicate}"
        for file in alignment_file_list:
            if tts_prefix in file.name or rnatts_prefix in file.name:
                valid_bam.add(file)

    if condition != "" and replicate != "":
        for file in alignment_file_list:
            if (
                f"RIBO-{condition}-{replicate}" in file.name
                or f"RNA-{condition}-{replicate}" in file.name
            ):
                valid_bam.add(file)

    if len(valid_bam) == 0:
        return None

    return valid_bam


def parse_offset_json(offset_json: Path) -> dict:
    """
    Read offset JSON file into a dictionary
    """

    if not offset_json.is_file():
        raise FileNotFoundError(
            f"Error: Offset JSON file does not exist! {offset_json}"
        )

    if offset_json.stat().st_size == 0:
        raise ValueError(
            f"Error: Offset JSON file is empty or invalid! {offset_json}"
        )

    with open(offset_json, "r", encoding="utf-8") as json_file:
        return json.load(json_file)


def write_gff_file(
    dataframe_out: pd.DataFrame,
    output_path: Path,
    output_filename: Path
) -> None:
    """
    write a dataframe to a gff file
    """

    file_path = output_path / output_filename
    file_path.parent.mkdir(parents=True, exist_ok=True)

    msg.message(f"Writing: {file_path}")
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("##gff-version 3\n")
    with open(file_path, "a", encoding="utf-8") as f:
        dataframe_out.to_csv(
            f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE
        )


def write_codon_interval_gff(
    output_path: Path,
    output_basename: str,
    codon_dict: dict
) -> None:
    """
    Create a gff3 file with all codon intervals.
    """

    if len(codon_dict) == 0:
        raise ValueError("Error: Codon dictionary is empty, cannot create codon interval gff file!")

    nTuple_gff = collections.namedtuple(
        "Pandas",
        [
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
    )

    rows = []
    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue
        chrom, mid, strand = key.split(":")
        start, stop = mid.split("-")

        attribute = f"ID={key};Peak_height={val[1]};Name={val[0]};Start_codon={val[0]}"

        rows.append(
            nTuple_gff(
                chrom,
                "ORFBounder",
                "codon_interval",
                int(start) + 1,
                int(stop) + 1,
                ".",
                strand,
                ".",
                attribute,
            )
        )

    df = pd.DataFrame.from_records(
        rows,
        columns=[
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
    )

    write_gff_file(df, output_path, output_basename)


def excel_writer(out_file_name: Path, data_frames: dict[str, pd.DataFrame]) -> None:
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only = [
        "Nucleotide_Seq",
        "Amino_Acid_Seq",
        "Start_codon",
        "Stop_codon",
        "Strand",
        "Codon_count",
    ]
    writer = pd.ExcelWriter(out_file_name, engine="xlsxwriter")
    for sheetname, df in data_frames.items():
        df.to_excel(writer, sheet_name=sheetname, index=False)
        worksheet = writer.sheets[sheetname]
        worksheet.freeze_panes(1, 0)
        for idx, col in enumerate(df):
            series = df[col]
            if col in header_only:
                max_len = len(str(series.name)) + 2
            else:
                max_len = (
                    max((series.astype(str).str.len().max(), len(str(series.name)))) + 1
                )

            worksheet.set_column(idx, idx, max_len)
    writer.close()


def write_results_to_gff(
    result_df: pd.DataFrame,
    output_path: Path,
    output_basename: str,
    split_gff: bool
) -> None:
    """
    write a gff file comtaining the ORFs from the csv,

    if split_gff == True then write one additional gff file for each gene_type
    """

    nTuple_gff = collections.namedtuple(
        "Pandas",
        [
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
    )

    gff_all = []
    gff_annotated = []
    gff_unannotated = []
    gff_near_annotated = []
    gff_internal_inframe = []
    gff_n_terminal = []
    gff_internal_out = []

    for row in result_df.itertuples(index=False, name=None):
        (
            gene_type,
            identifier,
            chrom,
            start,
            stop,
            strand,
            locus_tag,
            codon_count,
            rpm_tis,
            rpm_tts,
            rpm_ribo,
            start_codon,
            stop_codon,
        ) = row[0:13]

        attribute = f"ID={identifier};Name={locus_tag};Peak_height_TIS={rpm_tis};Peak_height_TTS={rpm_tts};Peak_height_RIBO={rpm_ribo};Start_codon={start_codon};Stop_codon={stop_codon};Codon_count={codon_count};Type={gene_type};"
        cur_tuple = nTuple_gff(
            chrom,
            "ORFBounder",
            "CDS",
            int(start),
            int(stop),
            ".",
            strand,
            ".",
            attribute,
        )

        gff_all.append(cur_tuple)
        if split_gff:
            if gene_type == "Annotated":
                gff_annotated.append(cur_tuple)
            elif gene_type == "Unannotated":
                gff_unannotated.append(cur_tuple)
            elif gene_type == "Near_Annotated":
                gff_near_annotated.append(cur_tuple)
            elif gene_type == "Internal_Inframe":
                gff_internal_inframe.append(cur_tuple)
            elif gene_type == "N-terminal_extension":
                gff_n_terminal.append(cur_tuple)
            elif gene_type == "Internal_OutofFrame":
                gff_internal_out.append(cur_tuple)

    msg.message("Generating gff files...")

    df_all = pd.DataFrame.from_records(
        gff_all,
        columns=[
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
    )
    write_gff_file(
        df_all, output_path, Path("result_gffs") / f"{output_basename}.gff"
    )

    if split_gff:
        df_annotated = pd.DataFrame.from_records(
            gff_annotated,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_annotated,
            output_path,
            Path("result_gffs") / f"{output_basename}_annotated.gff",
        )

        df_unannotated = pd.DataFrame.from_records(
            gff_unannotated,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_unannotated,
            output_path,
            Path("result_gffs") / f"{output_basename}_unannotated.gff",
        )

        df_near_annotated = pd.DataFrame.from_records(
            gff_near_annotated,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_near_annotated,
            output_path,
            Path("result_gffs") / f"{output_basename}_near_annotated.gff",
        )

        df_internal_inframe = pd.DataFrame.from_records(
            gff_internal_inframe,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_internal_inframe,
            output_path,
            Path("result_gffs") / f"{output_basename}_internal_inframe.gff",
        )

        df_n_terminal = pd.DataFrame.from_records(
            gff_n_terminal,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_n_terminal,
            output_path,
            Path("result_gffs") / f"{output_basename}_n_terminal.gff",
        )

        df_internal_out = pd.DataFrame.from_records(
            gff_internal_out,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_internal_out,
            output_path,
            Path("result_gffs") / f"{output_basename}_internal_out.gff",
        )
    msg.success("Done")


def write_results_to_table(
    df_results: pd.DataFrame,
    output_path: Path,
    output_basename: str
) -> None:
    """
    write a csv and xlsx file containing all information,
    """

    msg.message("Generating output_tables...")

    out_csv = output_path / f"{output_basename}.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    df_results.to_csv(out_csv, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    out_xlsx = output_path / f"{output_basename}.xlsx"
    df_dict = {"CDS": df_results}
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)

    excel_writer(out_xlsx, df_dict)
    msg.success("Done")
