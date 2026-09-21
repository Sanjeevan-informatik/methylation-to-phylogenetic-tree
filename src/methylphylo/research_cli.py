"""Research CLI; the original methylphylo command remains available."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio import Phylo

from .datasets import load_dataset, prepare
from .inference import RECIPES, catalog, compare_trees
from .research import benchmark, run_analysis, simulate, write_json


def parser():
    p = argparse.ArgumentParser(
        prog="methylphylo-research", description="Multi-data phylogenetic research workbench"
    )
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("catalog", help="List supported model/data combinations")
    run = sub.add_parser("infer", help="Infer a tree with explicit data/model assumptions")
    run.add_argument("input", type=Path)
    run.add_argument(
        "--data-type",
        required=True,
        choices=["dna", "vcf", "beta", "counts", "binary", "binary-fasta"],
    )
    run.add_argument("--recipe", required=True, choices=RECIPES)
    run.add_argument("--missing", choices=["complete", "retain"], default="complete")
    run.add_argument("--sites", choices=["all", "variable"], default="all")
    run.add_argument("--genotype-policy", choices=["homozygous", "phased"], default="homozygous")
    run.add_argument("--min-depth", type=int, default=5)
    run.add_argument("--low", type=float)
    run.add_argument("--high", type=float)
    run.add_argument("--allow-uncorrected", action="store_true")
    run.add_argument("--bootstrap", type=int, default=0)
    run.add_argument("--block-size", type=int, default=1)
    run.add_argument("--seed", type=int, default=1)
    run.add_argument("--threads", type=int, default=1)
    run.add_argument("--timeout", type=int, default=3600)
    run.add_argument("--iqtree")
    run.add_argument("--truth", type=Path)
    run.add_argument("--output", type=Path, required=True)
    sim = sub.add_parser("simulate", help="Generate independent-site data on a known tree")
    sim.add_argument("--kind", choices=["dna", "binary"], default="dna")
    sim.add_argument("--sites", type=int, default=2000)
    sim.add_argument("--seed", type=int, default=1)
    sim.add_argument("--missing", type=float, default=0.0)
    sim.add_argument("--error", type=float, default=0.0)
    sim.add_argument("--gain", type=float, default=0.2)
    sim.add_argument("--loss", type=float, default=0.1)
    sim.add_argument("--tree", type=Path)
    sim.add_argument("--output", type=Path, required=True)
    bench = sub.add_parser(
        "benchmark", help="Run a JSON manifest; record failures and comparable trees"
    )
    bench.add_argument("manifest", type=Path)
    bench.add_argument("--output", type=Path, required=True)
    bench.add_argument("--iqtree")
    compare = sub.add_parser(
        "compare", help="Compare unrooted topologies with identical tip labels"
    )
    compare.add_argument("left", type=Path)
    compare.add_argument("right", type=Path)
    compare.add_argument("--output", type=Path)
    serve = sub.add_parser("serve", help="Open the local research browser interface")
    serve.add_argument("--port", type=int, default=8765)
    return p


def main():
    p = parser()
    a = p.parse_args()
    try:
        if a.command == "catalog":
            print(json.dumps(catalog(), indent=2))
        elif a.command == "infer":
            data = load_dataset(
                a.input, a.data_type, genotype_policy=a.genotype_policy, min_depth=a.min_depth
            )
            data = prepare(data, missing=a.missing, sites=a.sites, low=a.low, high=a.high)
            report = run_analysis(
                data,
                a.recipe,
                a.output,
                seed=a.seed,
                threads=a.threads,
                executable=a.iqtree,
                timeout=a.timeout,
                allow_uncorrected=a.allow_uncorrected,
                replicates=a.bootstrap,
                block_size=a.block_size,
                truth=a.truth,
            )
            print(f"Wrote {a.output / 'report.json'} ({report['data']['sites']} retained sites)")
        elif a.command == "simulate":
            simulate(
                a.output,
                kind=a.kind,
                sites=a.sites,
                seed=a.seed,
                missing=a.missing,
                error=a.error,
                gain=a.gain,
                loss=a.loss,
                tree_path=a.tree,
            )
            print(f"Wrote known-tree simulation to {a.output}")
        elif a.command == "benchmark":
            rows = benchmark(a.manifest, a.output, a.iqtree)
            failed = sum(row["status"] != "ok" for row in rows)
            print(
                f"Wrote {a.output / 'summary.csv'}; {len(rows) - failed} succeeded, {failed} failed"
            )
            if failed:
                raise SystemExit(1)
        elif a.command == "compare":
            result = compare_trees(Phylo.read(a.left, "newick"), Phylo.read(a.right, "newick"))
            if a.output:
                write_json(a.output, result)
            print(json.dumps(result, indent=2))
        elif a.command == "serve":
            from .server import serve

            serve(a.port)
    except (ValueError, OSError) as exc:
        p.error(str(exc))


if __name__ == "__main__":
    main()
