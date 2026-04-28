"""Prepare and QC GSE95426 probe-level expression matrix.

This script performs basic quality control on the extracted GEO series
matrix, checks numeric expression values, handles duplicated probe IDs,
and saves a clean probe-level expression matrix for downstream analysis.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    """Run QC and prepare clean probe-level expression matrix."""
    expression_path = Path("data/processed/GSE95426_expression_matrix_raw.csv")
    metadata_path = Path("data/processed/GSE95426_metadata_verified.csv")

    processed_dir = Path("data/processed")
    table_dir = Path("results/tables")
    figure_dir = Path("results/figures")

    processed_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    if not expression_path.exists():
        raise FileNotFoundError(f"Expression matrix not found: {expression_path}")

    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    expression = pd.read_csv(expression_path)
    metadata = pd.read_csv(metadata_path)

    print("Raw expression shape:", expression.shape)
    print("Metadata shape:", metadata.shape)

    if "ID_REF" not in expression.columns:
        raise ValueError("Expected column 'ID_REF' not found in expression matrix.")

    sample_ids = metadata["sample_id"].tolist()
    expression_sample_cols = expression.columns.tolist()[1:]

    if expression_sample_cols != sample_ids:
        raise ValueError("Expression sample columns do not match metadata sample IDs.")

    print("Sample columns match metadata.")

    # Convert expression values to numeric.
    expression_numeric = expression.copy()
    for col in sample_ids:
        expression_numeric[col] = pd.to_numeric(expression_numeric[col], errors="coerce")

    missing_values = int(expression_numeric[sample_ids].isna().sum().sum())
    duplicated_probe_count = int(expression_numeric["ID_REF"].duplicated().sum())

    print("Missing expression values:", missing_values)
    print("Duplicated probe IDs:", duplicated_probe_count)

    # If duplicated probe IDs exist, average them.
    clean_expression = (
        expression_numeric
        .groupby("ID_REF", as_index=False)[sample_ids]
        .mean()
    )

    print("Clean probe-level expression shape:", clean_expression.shape)

    value_min = float(clean_expression[sample_ids].min().min())
    value_max = float(clean_expression[sample_ids].max().max())
    value_mean = float(clean_expression[sample_ids].stack().mean())
    value_median = float(clean_expression[sample_ids].stack().median())

    # Heuristic log2-scale check.
    if value_max <= 30:
        scale_interpretation = "likely_log2_or_normalized_expression"
    else:
        scale_interpretation = "possibly_raw_or_non_log_expression"

    qc_summary = pd.DataFrame(
        [
            {"metric": "raw_rows", "value": expression.shape[0]},
            {"metric": "raw_columns", "value": expression.shape[1]},
            {"metric": "clean_rows", "value": clean_expression.shape[0]},
            {"metric": "clean_columns", "value": clean_expression.shape[1]},
            {"metric": "sample_count", "value": len(sample_ids)},
            {"metric": "missing_values", "value": missing_values},
            {"metric": "duplicated_probe_ids", "value": duplicated_probe_count},
            {"metric": "expression_min", "value": value_min},
            {"metric": "expression_max", "value": value_max},
            {"metric": "expression_mean", "value": value_mean},
            {"metric": "expression_median", "value": value_median},
            {"metric": "scale_interpretation", "value": scale_interpretation},
        ]
    )

    clean_output_path = processed_dir / "GSE95426_expression_matrix_clean_probe_level.csv"
    qc_summary_path = table_dir / "GSE95426_qc_summary.csv"
    boxplot_path = figure_dir / "GSE95426_sample_boxplot_probe_level.png"

    clean_expression.to_csv(clean_output_path, index=False)
    qc_summary.to_csv(qc_summary_path, index=False)

    # Sample-level boxplot for expression distribution QC.
    plt.figure(figsize=(12, 6))
    clean_expression[sample_ids].boxplot(rot=45)
    plt.title("GSE95426 Probe-level Expression Distribution")
    plt.ylabel("Expression value")
    plt.xlabel("Sample")
    plt.tight_layout()
    plt.savefig(boxplot_path, dpi=300)
    plt.close()

    print("\nQC summary:")
    print(qc_summary)

    print("\nSaved outputs:")
    print(clean_output_path)
    print(qc_summary_path)
    print(boxplot_path)


if __name__ == "__main__":
    main()