#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import argparse

def csvs_to_parquet(folder_path: Path):
    """
    Convert all CSV files in a folder to Parquet files.

    Parameters
    ----------
    folder_path : Path
        Path to the folder containing CSV files.
    """

    if not folder_path.is_dir():
        raise FileNotFoundError(f"Folder does not exist: {folder_path}")

    csv_files = list(folder_path.glob("*.csv"))
    if not csv_files:
        print(f"No CSV files found in {folder_path}")
        return

    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        parquet_file = csv_file.with_suffix(".parquet")
        df.to_parquet(parquet_file, engine = "pyarrow", index=False)
        print(f"Converted {csv_file.name} → {parquet_file.name}")


def main():
    parser = argparse.ArgumentParser(
        description = "Convert all CSV files in a folder to Parquet."
    )
    parser.add_argument("folder", type = Path, 
                        help = "Path to folder containing CSV files.")
    args = parser.parse_args()

    csvs_to_parquet(args.folder)


if __name__ == "__main__":
    main()