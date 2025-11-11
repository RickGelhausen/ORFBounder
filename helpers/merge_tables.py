#!/usr/bin/env python3
"""
Script to merge multiple ORFBounder result tables into a single table.
"""
from pathlib import Path

import argparse
import pandas as pd
import numpy as np

SINGLE_FILE_VARIABLE_COLUMNS = [
    8,
    9,
    10,
    11,
    12,
    13,
    21,
    22,
    23,
    24,
]  # When dropping evidence column


def concat_cols(row, cols):
    """
    Concatenate column names with non-zero and non-NaN values
    """
    col_values = [
        "_".join(str(col).split("_")[:2])
        for col in cols
        if row[col] != 0 and not np.isnan(row[col])
    ]
    return ",".join(col_values)


def main():
    """
    Script to merge multiple ORFBounder result tables into a single table.
    """
    parser = argparse.ArgumentParser(description="Merge tables")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input tables")
    parser.add_argument("-o", "--output", required=True, help="Output table")
    args = parser.parse_args()

    table_paths = [Path(path) for path in args.input]
    name_info = [
        (path.stem.split("_")[0], path.stem.split("_")[1]) for path in table_paths
    ]

    dataframe_dict = {}
    for i, (name_prefix, name_suffix) in enumerate(name_info):
        dataframe_dict[f"{name_prefix}-{name_suffix}"] = pd.read_excel(
            table_paths[i], sheet_name=0
        )

    for key, df in dataframe_dict.items():
        df.drop(df.columns[14], axis=1, inplace=True)
        df_cols = [
            col
            for col in df.columns
            if col.endswith("_log2FC")
            or col.endswith("_peak_height")
            or col.endswith("_relative_density")
        ]
        dataframe_dict[key] = df.rename(
            columns={c: f"{key}_{c}" for c in df.columns if c in df_cols}
        )

    data_frames = [dataframe_dict[key] for key in list(dataframe_dict.keys())]

    new_df = pd.merge(
        data_frames[0],
        data_frames[1],
        on=[
            "Type",
            "Identifier",
            "Genome",
            "Start",
            "Stop",
            "Strand",
            "Locus_tag",
            "Codon_count",
            "Start_codon",
            "Stop_codon",
            "15nt_window",
            "Nucleotide_Seq",
            "Amino_Acid_Seq",
            "5'-distance",
            "3'-distance",
        ],
        how="outer",
    )
    drop_cols = [col for col in new_df.columns if col.endswith("_y")]
    new_df.drop(drop_cols, axis=1, inplace=True)
    rename_cols = [col for col in new_df.columns if col.endswith("_x")]
    new_df.rename(columns={c: c[:-2] for c in rename_cols}, inplace=True)

    for i in range(2, len(data_frames)):
        new_df = pd.merge(
            new_df,
            data_frames[i],
            on=[
                "Type",
                "Identifier",
                "Genome",
                "Start",
                "Stop",
                "Strand",
                "Locus_tag",
                "Codon_count",
                "Start_codon",
                "Stop_codon",
                "15nt_window",
                "Nucleotide_Seq",
                "Amino_Acid_Seq",
                "5'-distance",
                "3'-distance",
            ],
            how="outer",
        )
        drop_cols = [col for col in new_df.columns if col.endswith("_y")]
        new_df.drop(drop_cols, axis=1, inplace=True)
        rename_cols = [col for col in new_df.columns if col.endswith("_x")]
        new_df.rename(columns={c: c[:-2] for c in rename_cols}, inplace=True)

    peak_height_cols = [col for col in new_df.columns if col.endswith("_peak_height")]
    log2fc_cols = [col for col in new_df.columns if col.endswith("_log2FC")]
    relative_density_cols = [
        col for col in new_df.columns if col.endswith("_relative_density")
    ]
    rpkm_cols = [col for col in new_df.columns if col.endswith("_rpkm")]
    te_cols = [col for col in new_df.columns if col.endswith("_TE")]

    new_df["Evidence"] = new_df.apply(lambda row: concat_cols(row, log2fc_cols), axis=1)

    new_header = (
        [
            "Identifier",
            "Type",
            "Genome",
            "Start",
            "Stop",
            "Strand",
            "Locus_tag",
            "Codon_count",
        ]
        + [col for col in log2fc_cols]
        + [
            "Evidence",
            "Start_codon",
            "Stop_codon",
            "15nt_window",
            "Nucleotide_Seq",
            "Amino_Acid_Seq",
            "5'-distance",
            "3'-distance",
        ]
        + [col for col in rpkm_cols]
        + [col for col in te_cols]
        + [col for col in peak_height_cols]
        + [col for col in relative_density_cols]
    )

    # sort by Genome, Start, Stop, Strand
    new_df.sort_values(by=["Genome", "Start", "Stop", "Strand"], inplace=True)
    new_df = new_df[new_header]

    writer = pd.ExcelWriter(args.output, engine="xlsxwriter")
    new_df.to_excel(writer, index=False, freeze_panes=(1, 1))

    for column in new_df:
        if column not in ["Nucleotide_Seq", "Amino_Acid_Seq", "Evidence"]:
            column_width = (
                max(new_df[column].astype(str).map(len).max(), len(column)) + 2
            )
        else:
            column_width = len(column) + 2
        col_idx = new_df.columns.get_loc(column)
        writer.sheets["Sheet1"].set_column(col_idx, col_idx, column_width)
    writer.close()


if __name__ == "__main__":
    main()
