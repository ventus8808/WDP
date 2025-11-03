#!/usr/bin/env python3
"""
Stratification_extract.py

Reads all CSV files in `Result/brms_stratified/`, extracts rows where the RUCC column is
missing/empty, and writes a combined output file with blocks grouped by source file.
Each block begins with a single-line label derived from the source filename
(e.g. "Race White" or "Sex Male"), followed by the CSV header and rows that had no RUCC value.

Output path (explicitly requested):
    /Users/ventus/Repository/WDP/Result/brms_stratified/Combined_Res.csv

Usage (from project root):
    python Analysis/Stratification_extract.py
Or with explicit input/output:
    python Analysis/Stratification_extract.py --input Result/brms_stratified --output Result/brms_stratified/Combined_Res.csv

Behavior:
 - Skips non-CSV files and the Combined_Res.csv itself.
 - Treats RUCC as missing if the RUCC column is absent, NA/NaN, or an empty/whitespace string.
 - If a file has no rows with missing RUCC, the script writes the filename-derived label and a comment line
   indicating no rows were written for that block (so it's clear which files were processed).
 - Prints progress to stdout.

Note: This script reads CSVs with pandas and treats all columns as strings to preserve empty values.
"""

from pathlib import Path
import argparse
import sys

try:
    import pandas as pd
except Exception as e:
    print(
        "This script requires pandas. Install it with `pip install pandas` and retry."
    )
    raise

# Default paths (derived relative to this script's location)
HERE = Path(__file__).resolve().parent  # WDP/Code/Analysis
# Move up to the repository root (WDP) from WDP/Code/Analysis
DEFAULT_REPO_ROOT = HERE.parents[1]  # WDP
DEFAULT_INPUT_DIR = DEFAULT_REPO_ROOT / "Result" / "brms_stratified"
DEFAULT_OUTPUT_FILE = Path(
    "/Users/ventus/Repository/WDP/Result/brms_stratified/Combined_Res.csv"
)


def derive_label_from_filename(path: Path) -> str:
    """
    Derive a human-friendly label from the filename.

    Examples:
      C00_C97_brms_Race_White.csv -> Race White
      C00_C97_brms_Sex_Male.csv   -> Sex Male
    If `brms_` token is present, use the substring after it; otherwise fallback to the stem
    with underscores replaced by spaces.
    """
    stem = path.stem
    if "brms_" in stem:
        label = stem.split("brms_", 1)[1]
    else:
        # Fallback: remove leading tokens that look like codes (e.g., C00, C97)
        parts = stem.split("_")
        start = 0
        for i, tok in enumerate(parts):
            # consider token a code if it starts with 'C' followed by digits, or is purely numeric
            if tok.startswith("C") and tok[1:].isdigit():
                start = i + 1
                continue
            # otherwise assume we've reached the descriptive part
            start = i
            break
        label = "_".join(parts[start:]) if start < len(parts) else stem
    return label.replace("_", " ").strip()


def find_rucc_column(columns):
    """
    Find a RUCC-like column name among `columns` (case-insensitive).
    Returns the exact column name if found, otherwise None.
    """
    for col in columns:
        if col is None:
            continue
        if "rucc" == str(col).strip().lower():
            return col
    # fallback: any column that contains 'rucc' e.g. 'RUCC2013' or 'rucc_code'
    for col in columns:
        if col is None:
            continue
        if "rucc" in str(col).lower():
            return col
    return None


def process_file(path: Path):
    """
    Read CSV at `path` and return a tuple (label, filtered_df).
    - label: derived label string for the file
    - filtered_df: DataFrame of rows where RUCC is missing/empty; an empty DataFrame if none.
    If the file can't be read, returns (None, None)

    Additionally: after selecting rows with missing RUCC, drop any rows where the Model
    column contains the substring 'RUCC' (case-insensitive).
    """
    try:
        df = pd.read_csv(path, dtype=str)
    except Exception as e:
        print(f"Skipping {path.name}: failed to read CSV ({e})", file=sys.stderr)
        return None, None

    label = derive_label_from_filename(path)
    rucc_col = find_rucc_column(df.columns)
    if rucc_col is None:
        # No RUCC column found; return empty dataframe with same columns to keep headers consistent
        print(
            f"File {path.name}: RUCC column not found; will write a note for this file."
        )
        return label, pd.DataFrame(columns=df.columns)

    # Consider RUCC missing if NA/NaN or empty/whitespace
    mask_missing = df[rucc_col].isna() | (df[rucc_col].astype(str).str.strip() == "")
    filtered = df[mask_missing].copy()

    # Drop rows where the Model column contains 'RUCC' (case-insensitive), if a Model column exists.
    model_col = None
    for col in df.columns:
        if col is None:
            continue
        if "model" == str(col).strip().lower():
            model_col = col
            break
    # fallback: any column name containing 'model'
    if model_col is None:
        for col in df.columns:
            if col is None:
                continue
            if "model" in str(col).lower():
                model_col = col
                break

    if model_col is not None and not filtered.empty:
        # Align mask with filtered rows by using filtered's index on the original df
        model_values = df.loc[filtered.index, model_col].astype(str)
        mask_rucc_in_model = model_values.str.lower().str.contains("rucc", na=False)
        before_count = len(filtered)
        filtered = filtered[~mask_rucc_in_model]
        dropped = before_count - len(filtered)
        if dropped:
            print(
                f"{path.name}: dropped {dropped} rows where {model_col} contains 'RUCC'"
            )

    return label, filtered


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Combine all stratified BRMS CSVs and keep only rows whose Model first letter is 'E'."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default=str(DEFAULT_INPUT_DIR),
        help="Directory containing stratified CSV files (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=str(DEFAULT_OUTPUT_FILE),
        help="Combined output file path (default: %(default)s)",
    )
    parser.add_argument(
        "--model-first-letter",
        "-m",
        type=str,
        default="E",
        help="Keep only rows where the first letter of the Model column (after stripping) matches this value (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    input_dir = Path(args.input)
    output_file = Path(args.output)

    if not input_dir.exists() or not input_dir.is_dir():
        print(
            f"Input directory does not exist or is not a directory: {input_dir}",
            file=sys.stderr,
        )
        sys.exit(2)

    csv_files = sorted(
        [
            p
            for p in input_dir.iterdir()
            if p.suffix.lower() == ".csv" and p.name != output_file.name
        ]
    )
    if not csv_files:
        print(f"No CSV files found in {input_dir}", file=sys.stderr)
        return

    # Read all files and append a Source column
    frames = []
    for p in csv_files:
        try:
            df = pd.read_csv(p, dtype=str)
        except Exception as e:
            print(f"Skipping {p.name}: failed to read CSV ({e})", file=sys.stderr)
            continue
        df["Source"] = derive_label_from_filename(p)
        frames.append(df)
        print(f"Read {p.name}: {len(df)} rows")

    if not frames:
        print("No data frames read. Exiting.")
        return

    combined = pd.concat(frames, ignore_index=True, sort=False)
    print(f"Combined rows: {len(combined)}")

    # Identify Model column (case-insensitive)
    model_col = None
    for col in combined.columns:
        if col is None:
            continue
        if "model" == str(col).strip().lower():
            model_col = col
            break
    if model_col is None:
        for col in combined.columns:
            if col is None:
                continue
            if "model" in str(col).lower():
                model_col = col
                break

    # Filter rows: keep only those where the first non-whitespace character of the Model value is the target letter
    if model_col is None:
        print(
            "Warning: 'Model' column not found; writing combined file without filtering."
        )
    else:
        target = args.model_first_letter.upper()
        model_series = combined[model_col].astype(str).str.strip().str[:1].str.upper()
        mask_keep = model_series == target
        before = len(combined)
        combined = combined[mask_keep].copy()
        after = len(combined)
        print(
            f"Kept {after} rows where first letter of {model_col} is '{target}' (dropped {before - after})."
        )

    # Ensure parent dir of output exists
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Write blocks in the requested order and format:
    # Column headers at the top, then each block begins with the Source label on its own line,
    # followed by the block rows (no header) with columns separated by tabs. Blocks are separated by a single empty line.
    ordered_labels = [
        "Sex Male",
        "Sex Female",
        "Race White",
        "Race Black or African American",
        "Race Asian or Pacific Islander",
        "Race American Indian or Alaska Native",
    ]

    # Normalize Source column to string for matching
    combined["Source"] = combined["Source"].astype(str)

    with output_file.open("w", encoding="utf-8", newline="") as out_f:
        # Write column headers at the top, tab-separated. Exclude the 'Source' column.
        cols = [c for c in combined.columns if c != "Source"]
        out_f.write("\t".join(cols) + "\n\n")

        first_block = True
        for label in ordered_labels:
            # select rows where Source exactly matches the label (strip whitespace)
            block_df = combined[combined["Source"].str.strip() == label]
            if block_df.empty:
                # If no rows for this label, still write the label and skip rows,
                # or skip entirely? We'll skip empty blocks to keep output concise.
                continue

            if not first_block:
                out_f.write("\n")  # blank line between blocks
            first_block = False

            # Write the label line
            out_f.write(f"{label}\n\n")

            # Write rows, tab-separated.
            for _, row in block_df.iterrows():
                vals = ["" if (pd.isna(row.get(c))) else str(row.get(c)) for c in cols]
                out_f.write("\t".join(vals) + "\n")

    print(f"Combined results written to: {output_file}")


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
