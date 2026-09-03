"""Command-line entry point for MethylPhylo."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .core import analyze, read_matrix


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="methylphylo",
        description="Build an exploratory UPGMA tree from a DNA-methylation beta-value matrix.",
    )
    parser.add_argument(
        "matrix", type=Path, help="CSV or TSV matrix; rows are CpGs and columns are samples"
    )
    parser.add_argument(
        "--top-loci", type=int, default=500, help="most-variable loci to retain (default: 500)"
    )
    parser.add_argument(
        "--metric",
        choices=("euclidean", "manhattan", "correlation"),
        default="euclidean",
        help="sample distance metric (default: euclidean)",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("methylphylo-results"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    matrix = read_matrix(args.matrix)
    result = analyze(matrix, top_loci=args.top_loci, metric=args.metric)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "tree.nwk").write_text(f"{result.newick}\n", encoding="utf-8")
    report = {
        "input": str(args.matrix),
        "metric": result.metric,
        "samples": list(result.processed.matrix.samples),
        "input_loci": result.processed.input_loci,
        "selected_loci": len(result.processed.matrix.loci),
        "missing_values": result.processed.missing_values,
        "missing_rate": result.processed.missing_rate,
        "dropped_loci": result.processed.dropped_loci,
        "newick": result.newick,
        "distance_matrix": result.distances.round(6).tolist(),
    }
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {args.output_dir / 'tree.nwk'} and {args.output_dir / 'report.json'}")


if __name__ == "__main__":
    main()
