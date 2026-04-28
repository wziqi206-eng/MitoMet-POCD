"""Extract and validate GSE95426 expression matrix.

This script validates sample metadata, extracts the expression matrix
from the GEO series matrix file, and saves processed CSV files.
"""

import gzip
from io import StringIO
from pathlib import Path

import pandas as pd


def get_series_line(prefix: str, lines: list[str]) -> list[str]:
    """Return lines that start with a given prefix."""
    return [line for line in lines if line.startswith(prefix)]


def parse_geo_metadata_line(line: str) -> list[str]:
    """Parse a GEO metadata line separated by tabs."""
    return [item.strip().strip('"') for item in line.split("\t")[1:]]


def main() -> None:
    """Run metadata validation and expression matrix extraction."""
    metadata_path = Path("data/metadata/GSE95426_metadata.csv")
    series_matrix_path = Path("data/raw/GSE95426/GSE95426_series_matrix.txt.gz")
    output_dir = Path("data/processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    if not series_matrix_path.exists():
        raise FileNotFoundError(f"Series matrix file not found: {series_matrix_path}")

    metadata = pd.read_csv(metadata_path)
    print("Metadata shape:", metadata.shape)

    group_counts = metadata["group"].value_counts()
    print("\nGroup counts:")
    print(group_counts)

    assert group_counts.get("control", 0) == 6, "Control group sample count is not 6."
    assert group_counts.get("POCD", 0) == 6, "POCD group sample count is not 6."

    with gzip.open(series_matrix_path, "rt", encoding="utf-8", errors="replace") as file:
        series_lines = file.read().splitlines()

    geo_accession_line = get_series_line("!Sample_geo_accession", series_lines)[0]
    sample_title_line = get_series_line("!Sample_title", series_lines)[0]

    geo_accessions = parse_geo_metadata_line(geo_accession_line)
    sample_titles = parse_geo_metadata_line(sample_title_line)

    metadata_sample_ids = metadata["sample_id"].tolist()
    metadata_sample_titles = metadata["sample_title"].tolist()

    assert metadata_sample_ids == geo_accessions, (
        "Sample IDs in metadata do not match GEO series matrix."
    )
    assert metadata_sample_titles == sample_titles, (
        "Sample titles in metadata do not match GEO series matrix."
    )

    print("\nMetadata validation passed.")
    print("Sample IDs and sample titles match GEO series matrix.")

    table_start = series_lines.index("!series_matrix_table_begin") + 1
    table_end = series_lines.index("!series_matrix_table_end")

    expression_text = "\n".join(series_lines[table_start:table_end])
    expression_matrix = pd.read_csv(StringIO(expression_text), sep="\t")

    print("\nExpression matrix shape:", expression_matrix.shape)

    expression_sample_columns = expression_matrix.columns.tolist()[1:]
    assert expression_sample_columns == metadata_sample_ids, (
        "Expression matrix sample columns do not match metadata sample IDs."
    )

    print("Expression matrix sample columns match metadata sample IDs.")

    expression_output_path = output_dir / "GSE95426_expression_matrix_raw.csv"
    metadata_output_path = output_dir / "GSE95426_metadata_verified.csv"

    expression_matrix.to_csv(expression_output_path, index=False)
    metadata.to_csv(metadata_output_path, index=False)

    print("\nSaved outputs:")
    print(expression_output_path)
    print(metadata_output_path)


if __name__ == "__main__":
    main()