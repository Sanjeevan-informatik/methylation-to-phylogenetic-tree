"""Reproducible runs, simulations, and benchmark manifests."""

from __future__ import annotations

import csv
import html
import io
import json
import platform
import time
from datetime import datetime, timezone
from pathlib import Path

import Bio
import numpy as np
import scipy
from Bio import Phylo
from scipy.linalg import expm

from . import __version__
from .datasets import Dataset, load_dataset, prepare, write_fasta
from .inference import bootstrap, compare_trees, infer, newick, split_set


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def tree_svg(tree):
    """Topology-only preview; Newick preserves the actual branch lengths."""
    tips = tree.get_terminals()
    width, height = 1000, max(240, len(tips) * 28 + 90)
    ys = {tip: 55 + i * 28 for i, tip in enumerate(tips)}
    depths = tree.depths(unit_branch_lengths=True)
    max_depth = max(depths.values()) or 1
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">',
        "<title>Tree topology (branch lengths not to scale)</title>",
        '<rect width="100%" height="100%" fill="#fff"/>',
        '<g font-family="sans-serif" font-size="13" fill="#163042">',
        '<text x="24" y="24">Tree topology · branch lengths not to scale</text>',
    ]

    def draw(node):
        x = 25 + depths[node] / max_depth * 700
        if node.is_terminal():
            elements.append(f'<text x="{x + 8}" y="{ys[node] + 4}">{html.escape(node.name)}</text>')
        else:
            for child in node.clades:
                draw(child)
            ys[node] = sum(ys[c] for c in node.clades) / len(node.clades)
            lo, hi = min(ys[c] for c in node.clades), max(ys[c] for c in node.clades)
            elements.append(f'<path d="M{x} {lo} V{hi}" stroke="#247c88" fill="none"/>')
            for child in node.clades:
                cx = 25 + depths[child] / max_depth * 700
                elements.append(f'<path d="M{x} {ys[child]} H{cx}" stroke="#247c88"/>')
            if node.confidence is not None:
                elements.append(
                    f'<text x="{x + 3}" y="{ys[node] - 5}">{node.confidence:.0f}</text>'
                )

    draw(tree.root)
    elements.append("</g></svg>")
    return "".join(elements)


def run_analysis(
    data,
    recipe,
    output,
    *,
    seed=1,
    threads=1,
    executable=None,
    timeout=3600,
    allow_uncorrected=False,
    replicates=0,
    block_size=1,
    truth=None,
):
    output = Path(output)
    if replicates < 0 or block_size < 1:
        raise ValueError("Replicates must be nonnegative and block size positive")
    if output.exists() and any(output.iterdir()):
        raise ValueError(
            "Output directory is not empty; choose a fresh directory to preserve provenance"
        )
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    engine_options = {
        "threads": threads,
        "executable": executable,
        "timeout": timeout,
        "allow_uncorrected": allow_uncorrected,
    }
    tree, details = infer(data, recipe, output / "engine", seed=seed, **engine_options)
    inference_seconds = time.perf_counter() - started
    if replicates:
        tree, support, trees = bootstrap(
            data,
            recipe,
            tree,
            output / "bootstrap",
            replicates,
            seed=seed,
            block_size=block_size,
            **engine_options,
        )
        details["bootstrap"] = support
        (output / "bootstrap.nwk").write_text("\n".join(trees) + "\n")
        if support["failures"]:
            details["warnings"].append(
                "Some bootstrap replicates failed; support is conditional on successful runs."
            )
    if truth is not None:
        truth_tree = Phylo.read(str(truth), "newick")
        details["reference_comparison"] = compare_trees(tree, truth_tree)
        details["reference_tree"] = str(Path(truth).resolve())
        import hashlib

        details["reference_sha256"] = hashlib.sha256(Path(truth).read_bytes()).hexdigest()
    # Export the exact retained observations, not an imputed or undocumented approximation.
    if data.kind == "beta":
        with (output / "retained.tsv").open("w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["locus", *data.names])
            for name, row in zip(data.loci, data.values.T):
                writer.writerow([name, *["NA" if np.isnan(x) else float(x) for x in row]])
    else:
        write_fasta(data, output / "retained.fasta")
    (output / "retained-sites.tsv").write_text(
        "alignment_column\tsource_site\n"
        + "".join(f"{i + 1}\t{loc}\n" for i, loc in enumerate(data.loci))
    )
    nwk = newick(tree)
    (output / "tree.nwk").write_text(nwk + "\n")
    (output / "tree.svg").write_text(tree_svg(tree))
    report = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "software": {
            "methylphylo": __version__,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "biopython": Bio.__version__,
            "scipy": scipy.__version__,
        },
        "data": {
            "kind": data.kind,
            "samples": list(data.names),
            "sites": len(data.loci),
            "missing_fraction": float(data.missing.mean()),
            "sha256": data.fingerprint(),
            "provenance": data.provenance,
        },
        "settings": {
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout,
            "allow_uncorrected": allow_uncorrected,
            "bootstrap_replicates": replicates,
            "block_size": block_size,
        },
        "inference_seconds": inference_seconds,
        "total_seconds": time.perf_counter() - started,
        "rooting": "UPGMA ultrametric"
        if tree.rooted
        else "unrooted; Newick root is representational",
        "newick": nwk,
        **details,
    }
    write_json(output / "report.json", report)
    return report


def simulate(
    output,
    *,
    kind="dna",
    sites=2000,
    seed=1,
    missing=0.0,
    error=0.0,
    gain=0.2,
    loss=0.1,
    tree_path=None,
):
    if kind not in {"dna", "binary"} or sites < 1 or seed < 1:
        raise ValueError("Simulation requires dna/binary, positive sites and seed")
    if not 0 <= missing < 1 or not 0 <= error < 1:
        raise ValueError("Missingness and error must lie in [0, 1)")
    if not np.isfinite([gain, loss]).all() or gain <= 0 or loss <= 0:
        raise ValueError("Gain and loss rates must be finite and positive")
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Simulation output must be an empty directory")
    tree = (
        Phylo.read(str(tree_path), "newick")
        if tree_path
        else Phylo.read(io.StringIO("((A:0.05,B:0.05):0.15,(C:0.05,D:0.05):0.15);"), "newick")
    )
    names = tuple(t.name for t in tree.get_terminals())
    split_set(tree)
    for c in tree.find_clades():
        if c is not tree.root and (
            c.branch_length is None or not np.isfinite(c.branch_length) or c.branch_length < 0
        ):
            raise ValueError("Every non-root simulation branch needs a finite nonnegative length")
    if kind == "dna":
        q = np.full((4, 4), 1 / 3)
        np.fill_diagonal(q, -1)
        pi = np.full(4, 0.25)
    else:
        q = np.array([[-gain, gain], [loss, -loss]])
        pi = np.array([loss, gain]) / (gain + loss)
    rng = np.random.default_rng(seed)
    observations = {}

    def evolve(node, states):
        if node.is_terminal():
            observations[node.name] = states.copy()
        for child in node.clades:
            p = expm(q * child.branch_length)
            p = np.clip(p, 0, 1)
            p /= p.sum(axis=1, keepdims=True)
            draws = rng.random(sites)
            next_states = (draws[:, None] > np.cumsum(p[states], axis=1)).sum(axis=1)
            evolve(child, np.minimum(next_states, len(pi) - 1))

    evolve(tree.root, rng.choice(len(pi), sites, p=pi))
    values = np.array([observations[n] for n in names])
    corrupt = rng.random(values.shape) < error
    offsets = rng.integers(1, len(pi), values.shape)
    values[corrupt] = (values[corrupt] + offsets[corrupt]) % len(pi)
    values[rng.random(values.shape) < missing] = -1
    data = Dataset(names, tuple(str(i + 1) for i in range(sites)), values, kind)
    output.mkdir(parents=True, exist_ok=True)
    write_fasta(data, output / "alignment.fasta")
    (output / "truth.nwk").write_text(newick(tree) + "\n")
    write_json(
        output / "simulation.json",
        {
            "kind": kind,
            "model": "JC69" if kind == "dna" else "asymmetric binary CTMC",
            "sites": sites,
            "seed": seed,
            "missing": missing,
            "error": error,
            "gain": gain if kind == "binary" else None,
            "loss": loss if kind == "binary" else None,
            "branch_units": "expected substitutions/site"
            if kind == "dna"
            else "time in specified rate units",
            "root_distribution": pi.tolist(),
            "rate_matrix": q.tolist(),
            "truth": newick(tree),
            "limitations": "Independent sites, shared topology, stationary root; no recombination, mixtures or selection",
        },
    )
    return data


def benchmark(manifest_path, output, executable=None):
    manifest_path, output = Path(manifest_path), Path(output)
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest.get("runs"), list) or not manifest["runs"]:
        raise ValueError("Manifest needs a nonempty runs list")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Benchmark output must be empty")
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "manifest.json", manifest)
    rows, reports = [], []
    ids = set()
    for entry in manifest["runs"]:
        run_id = entry.get("id", "")
        from .datasets import SAFE_NAME

        if not SAFE_NAME.fullmatch(run_id) or run_id in {".", ".."} or run_id in ids:
            raise ValueError("Benchmark run IDs must be unique safe names")
        ids.add(run_id)
        row = {"id": run_id, "recipe": entry.get("recipe", ""), "status": "failed"}
        try:
            path = (manifest_path.parent / entry["input"]).resolve()
            data = load_dataset(
                path,
                entry["data_type"],
                genotype_policy=entry.get("genotype_policy", "homozygous"),
                min_depth=entry.get("min_depth", 5),
            )
            data = prepare(
                data,
                missing=entry.get("missing", "complete"),
                sites=entry.get("sites", "all"),
                low=entry.get("low"),
                high=entry.get("high"),
            )
            truth = manifest_path.parent / entry["truth"] if entry.get("truth") else None
            report = run_analysis(
                data,
                entry["recipe"],
                output / run_id,
                seed=entry.get("seed", 1),
                threads=entry.get("threads", 1),
                executable=executable,
                timeout=entry.get("timeout", 3600),
                allow_uncorrected=entry.get("allow_uncorrected", False),
                replicates=entry.get("replicates", 0),
                block_size=entry.get("block_size", 1),
                truth=truth,
            )
            row.update(
                status="ok",
                data_sha256=report["data"]["sha256"],
                sites=len(data.loci),
                inference_seconds=report["inference_seconds"],
                normalized_rf=report.get("reference_comparison", {}).get("normalized_rf", ""),
            )
            reports.append((run_id, report))
        except (ValueError, OSError, KeyError) as exc:
            row["error"] = str(exc)
        rows.append(row)
    fields = [
        "id",
        "recipe",
        "status",
        "data_sha256",
        "sites",
        "inference_seconds",
        "normalized_rf",
        "error",
    ]
    with (output / "summary.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    comparisons = []
    for i, (left_id, left) in enumerate(reports):
        for right_id, right in reports[:i]:
            if set(left["data"]["samples"]) != set(right["data"]["samples"]):
                continue
            comparisons.append(
                {
                    "left": left_id,
                    "right": right_id,
                    "same_observations": left["data"]["sha256"] == right["data"]["sha256"],
                    **compare_trees(
                        Phylo.read(io.StringIO(left["newick"]), "newick"),
                        Phylo.read(io.StringIO(right["newick"]), "newick"),
                    ),
                }
            )
    write_json(output / "comparisons.json", comparisons)
    write_json(output / "summary.json", rows)
    return rows
