"""
Create a gff out of the final table and color it according to the TIS/TTS status
"""
from pathlib import Path

import argparse
import pandas as pd
import numpy as np


COLOR_DICT = { "TIS_only" : "0,0,255", "TTS_only" : "255,0,255", "both" : "180,180,180"}

def create_gff(df, output):
    """
    Create a gff out of the final table and color it according to the TIS/TTS status
    """
    out_rows = []
    for row in df.itertuples(index=False):
        identifier = getattr(row, "Identifier")
        genome = getattr(row, "Genome")
        start = getattr(row, "Start")
        stop = getattr(row, "Stop")
        strand = getattr(row, "Strand")
        evidence = getattr(row, "Evidence").lower()

        if "tis" in evidence and "tts" in evidence:
            color = COLOR_DICT["both"]
        elif "tis" in evidence:
            color = COLOR_DICT["TIS_only"]
        elif "tts" in evidence:
            color = COLOR_DICT["TTS_only"]
        else:
            raise ValueError("Error. Empty evidence row.")

        log2fc_cols = [col for col in df.columns if col.endswith("_log2FC")]
        log2fc_col_indices = [df.columns.get_loc(col) for col in log2fc_cols]
        log2fc_attr = ""
        for i, col in enumerate(log2fc_col_indices):
            if row[col] != 0 and not np.isnan(row[col]):
                condition = log2fc_cols[i].split('_')[1].split("-")[1] + log2fc_cols[i].split('_')[1].split("-")[2]
                col_name = f"{log2fc_cols[i].split('-')[0].lower()}_{log2fc_cols[i].split('-')[1].split('_')[0]}_{condition}_{log2fc_cols[i].split('_')[-1]}"
                log2fc_attr += f"{col_name}={float(row[col]):.4f};"

        attributes = f"ID={identifier};Name={identifier};color={color};" + log2fc_attr
        out_rows.append([genome, "ORFBounder", "CDS", start, stop, ".", strand, ".", attributes])

    out_df = pd.DataFrame(out_rows, columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])
    with open(output, 'w', encoding='utf-8') as f:
        f.write("##gff-version 3\n")
        out_df.to_csv(f, sep="\t", index=False, header=False)


def main():
    """
    Main function to parse arguments and create the gff file.
    """
    parser = argparse.ArgumentParser(description='Create a gff out of the final table and color it according to the TIS/TTS status')
    parser.add_argument('-i', '--input', type=Path, required=True, help='Input table')
    parser.add_argument('-o', '--output', type=Path, required=True, help='Output table')
    args = parser.parse_args()

    df = pd.read_excel(args.input, sheet_name=0)

    create_gff(df, args.output)

if __name__ == '__main__':
    main()