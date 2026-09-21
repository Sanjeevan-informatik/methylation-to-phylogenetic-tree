"""Distance/parsimony inference and a validated IQ-TREE adapter."""

from __future__ import annotations

import copy
import math
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import (
    DistanceMatrix,
    DistanceTreeConstructor,
    NNITreeSearcher,
    ParsimonyScorer,
    ParsimonyTreeConstructor,
)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .core import upgma
from .datasets import Dataset, write_fasta


@dataclass(frozen=True)
class Recipe:
    label: str
    kind: str
    method: str
    model: str
    interpretation: str


RECIPES = {
    "beta-upgma": Recipe("Beta UPGMA", "beta", "upgma", "euclidean", "similarity"),
    "beta-nj": Recipe("Beta neighbour joining", "beta", "nj", "euclidean", "similarity"),
    "beta-nj-manhattan": Recipe("Beta NJ Manhattan", "beta", "nj", "manhattan", "similarity"),
    "dna-nj-p": Recipe("DNA NJ p-distance", "dna", "nj", "p", "sequence relationships"),
    "dna-nj-jc69": Recipe("DNA NJ JC69", "dna", "nj", "jc69", "sequence relationships"),
    "dna-nj-k2p": Recipe("DNA NJ K2P", "dna", "nj", "k2p", "sequence relationships"),
    "dna-parsimony": Recipe(
        "DNA parsimony (NNI heuristic)", "dna", "parsimony", "fitch", "sequence relationships"
    ),
    "dna-ml-jc69": Recipe("DNA ML JC69", "dna", "ml", "JC", "sequence relationships"),
    "dna-ml-hky": Recipe("DNA ML HKY+G4", "dna", "ml", "HKY+G4", "sequence relationships"),
    "dna-ml-gtr": Recipe("DNA ML GTR+G4", "dna", "ml", "GTR+G4", "sequence relationships"),
    "dna-ml-gtr-asc": Recipe(
        "SNP ML GTR+G4+ASC", "dna", "ml", "GTR+G4+ASC", "sequence relationships"
    ),
    "binary-nj": Recipe("Binary NJ Hamming", "binary", "nj", "p", "state relationships"),
    "binary-parsimony": Recipe(
        "Binary parsimony (NNI heuristic)", "binary", "parsimony", "fitch", "state relationships"
    ),
    "binary-ml-jc2": Recipe("Binary symmetric ML", "binary", "ml", "JC2", "state relationships"),
    "binary-ml-gtr2": Recipe("Binary asymmetric ML", "binary", "ml", "GTR2", "state relationships"),
    "binary-ml-jc2-asc": Recipe(
        "Binary symmetric ML + ASC", "binary", "ml", "JC2+ASC", "state relationships"
    ),
    "binary-ml-gtr2-asc": Recipe(
        "Binary asymmetric ML + ASC", "binary", "ml", "GTR2+ASC", "state relationships"
    ),
}


def catalog():
    return {key: asdict(value) for key, value in RECIPES.items()}


def distances(data: Dataset, model: str):
    n = len(data.names)
    matrix, overlap = np.zeros((n, n)), np.zeros((n, n), dtype=int)
    np.fill_diagonal(overlap, (~data.missing).sum(axis=1))
    for i in range(n):
        for j in range(i):
            valid = ~(data.missing[i] | data.missing[j])
            if not valid.any():
                raise ValueError(
                    f"No shared observed sites for {data.names[i]} and {data.names[j]}"
                )
            a, b = data.values[i, valid], data.values[j, valid]
            overlap[i, j] = overlap[j, i] = len(a)
            if model == "euclidean":
                d = float(np.sqrt(np.mean((a - b) ** 2)))
            elif model == "manhattan":
                d = float(np.mean(np.abs(a - b)))
            else:
                p = float(np.mean(a != b))
                if model == "p":
                    d = p
                elif model == "jc69":
                    if p >= 0.75:
                        raise ValueError(
                            "JC69 saturation: mismatch fraction >= 0.75; no finite distance"
                        )
                    d = -0.75 * math.log1p(-4 * p / 3)
                elif model == "k2p":
                    transitions = ((a == 0) & (b == 2)) | ((a == 2) & (b == 0))
                    transitions |= ((a == 1) & (b == 3)) | ((a == 3) & (b == 1))
                    ts = float(np.mean(transitions))
                    tv = p - ts
                    if 1 - 2 * ts - tv <= 0 or 1 - 2 * tv <= 0:
                        raise ValueError("K2P saturation: no finite distance")
                    d = -0.5 * math.log(1 - 2 * ts - tv) - 0.25 * math.log(1 - 2 * tv)
                else:
                    raise ValueError(f"Unsupported distance: {model}")
            matrix[i, j] = matrix[j, i] = d
    return matrix, overlap


def distance_tree(data, model, method="nj"):
    matrix, overlap = distances(data, model)
    dm = DistanceMatrix(list(data.names), [list(matrix[i, : i + 1]) for i in range(len(matrix))])
    constructor = DistanceTreeConstructor()
    if method == "nj":
        tree = constructor.nj(dm)
    else:
        # Use sample-weighted cluster averages, including unequal cluster sizes.
        root = upgma(data.names, matrix)

        def convert(node, parent_height=None):
            return Phylo.BaseTree.Clade(
                branch_length=0.0 if parent_height is None else parent_height - node.height,
                name=node.name,
                clades=[]
                if node.children is None
                else [convert(c, node.height) for c in node.children],
            )

        tree = Phylo.BaseTree.Tree(root=convert(root), rooted=True)
    tree.rooted = method == "upgma"
    return tree, {"distance_matrix": matrix.tolist(), "pairwise_observed_sites": overlap.tolist()}


def parsimony_tree(data):
    if data.missing.any():
        raise ValueError("Parsimony requires complete sites; choose --missing complete")
    alphabet = "ACGT" if data.kind == "dna" else "01"
    alignment = MultipleSeqAlignment(
        [
            SeqRecord(Seq("".join(alphabet[int(x)] for x in row)), id=name)
            for name, row in zip(data.names, data.values)
        ]
    )
    start, _ = distance_tree(data, "p")
    scorer = ParsimonyScorer()
    tree = ParsimonyTreeConstructor(NNITreeSearcher(scorer), start).build_tree(alignment)
    score = scorer.get_score(tree, alignment)
    # Parsimony does not estimate evolutionary branch lengths.
    for clade in tree.find_clades():
        clade.branch_length = None
    tree.rooted = False
    return tree, {
        "parsimony_score": int(score),
        "search": "NNI from NJ start; not a global optimum guarantee",
    }


def iqtree_binary(explicit=None):
    if explicit:
        found = shutil.which(explicit)
        if not found:
            raise ValueError(f"IQ-TREE executable not found: {explicit}")
        return found
    for name in ("iqtree2", "iqtree", "iqtree3"):
        if found := shutil.which(name):
            return found
    raise ValueError(
        "IQ-TREE is required for ML recipes. Install it or select a NJ/UPGMA/parsimony recipe."
    )


def validate_recipe(data, recipe_id, allow_uncorrected=False):
    if recipe_id not in RECIPES:
        raise ValueError(f"Unknown recipe: {recipe_id}")
    recipe = RECIPES[recipe_id]
    if data.kind != recipe.kind:
        raise ValueError(f"{recipe_id} needs {recipe.kind} data, received {data.kind}")
    variable_only = all(
        len(set(col[~miss])) > 1 for col, miss in zip(data.values.T, data.missing.T)
    )
    if "+ASC" in recipe.model:
        if data.missing.any():
            raise ValueError(
                "ASC currently requires complete sites to avoid ambiguous constant patterns"
            )
        if not variable_only:
            raise ValueError("ASC requires variable-only data; use --sites variable explicitly")
    elif recipe.method == "ml" and variable_only and not allow_uncorrected:
        raise ValueError(
            "Variable-only ML needs an ASC recipe or --allow-uncorrected for a labelled bias control"
        )
    return recipe, variable_only


def infer(
    data,
    recipe_id,
    output,
    *,
    seed=1,
    threads=1,
    executable=None,
    timeout=3600,
    allow_uncorrected=False,
):
    if seed < 1 or threads < 1 or timeout <= 0:
        raise ValueError("Seed, threads and timeout must be positive")
    recipe, variable_only = validate_recipe(data, recipe_id, allow_uncorrected)
    warnings = list(data.warnings)
    if data.kind == "beta":
        warnings.append("Similarity tree: no claim of cellular ancestry or evolutionary time.")
    if data.kind != "beta" and variable_only and "+ASC" not in recipe.model:
        warnings.append(
            "Variable-only uncorrected analysis; distances/branch lengths may be biased."
        )
    details = {}
    if recipe.method in {"nj", "upgma"}:
        tree, details = distance_tree(data, recipe.model, recipe.method)
    elif recipe.method == "parsimony":
        tree, details = parsimony_tree(data)
    else:
        binary = iqtree_binary(executable)
        output = Path(output).resolve()
        output.mkdir(parents=True, exist_ok=True)
        alignment = output / "alignment.fasta"
        write_fasta(data, alignment)
        prefix = output / "iqtree"
        if any(output.glob("iqtree.*")):
            raise ValueError("IQ-TREE output already exists; use a fresh output directory")
        version = subprocess.run(
            [binary, "--version"], capture_output=True, text=True, timeout=15, check=True
        ).stdout.strip()
        cmd = [
            binary,
            "-s",
            str(alignment),
            "-st",
            "DNA" if data.kind == "dna" else "BIN",
            "-m",
            recipe.model,
            "-nt",
            str(threads),
            "-seed",
            str(seed),
            "-pre",
            str(prefix),
            "-keep-ident",
        ]
        with (output / "engine.log").open("w") as log:
            try:
                subprocess.run(
                    cmd, stdout=log, stderr=subprocess.STDOUT, timeout=timeout, check=True
                )
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
                raise ValueError(
                    f"IQ-TREE failed or timed out; inspect {output / 'engine.log'}"
                ) from exc
        tree_path = output / "iqtree.treefile"
        if not tree_path.exists():
            raise ValueError("IQ-TREE did not produce a treefile")
        tree = Phylo.read(tree_path, "newick")
        tree.rooted = False
        details = {
            "command": cmd,
            "engine_version": version,
            "engine_report": str(output / "iqtree.iqtree"),
        }
    if {t.name for t in tree.get_terminals()} != set(data.names):
        raise ValueError("Inferred tree has different sample labels than the input")
    if any(c.branch_length is not None and c.branch_length < -1e-10 for c in tree.find_clades()):
        warnings.append(
            "Negative NJ branch lengths retained in Newick; indicates non-additive distances."
        )
    return tree, {"recipe": recipe_id, **asdict(recipe), **details, "warnings": warnings}


def newick(tree):
    """Keep absent parsimony branch lengths absent, and retain support and negative NJ edges."""

    def render(node):
        if node.is_terminal():
            label = "'" + node.name.replace("'", "''") + "'"
        else:
            label = "(" + ",".join(render(c) for c in node.clades) + ")"
            if node.confidence is not None:
                label += f"{node.confidence:.8g}"
        if node.branch_length is not None:
            label += f":{node.branch_length:.10g}"
        return label

    return render(tree.root) + ";"


def split_set(tree):
    names = [tip.name for tip in tree.get_terminals()]
    if None in names or len(set(names)) != len(names):
        raise ValueError("Tree tips must have unique names")
    all_names = frozenset(names)
    splits = set()
    for node in tree.get_nonterminals():
        side = frozenset(t.name for t in node.get_terminals())
        other = all_names - side
        if min(len(side), len(other)) > 1:
            splits.add(min(tuple(sorted(side)), tuple(sorted(other))))
    return all_names, splits


def compare_trees(left, right):
    names, a = split_set(left)
    names_b, b = split_set(right)
    if names != names_b:
        raise ValueError("Tree comparison requires exactly the same tip labels")
    rf = len(a ^ b)
    denom = len(a) + len(b)
    return {
        "rf": rf,
        "normalized_rf": rf / denom if denom else 0.0,
        "shared_splits": len(a & b),
        "left_splits": len(a),
        "right_splits": len(b),
        "definition": "Unrooted RF / sum of observed nontrivial splits; branch lengths ignored",
    }


def bootstrap(data, recipe_id, reference, output, replicates, *, seed=1, block_size=1, **kwargs):
    if replicates < 1 or block_size < 1:
        raise ValueError("Bootstrap replicates and block size must be positive")
    rng = np.random.default_rng(seed)
    sites = data.values.shape[1]
    # Blocks follow the input order and never wrap across the ends.
    blocks = [np.arange(i, min(i + block_size, sites)) for i in range(0, sites, block_size)]
    counts, failures, trees = {}, [], []
    for i in range(replicates):
        indices = np.concatenate([blocks[j] for j in rng.integers(0, len(blocks), len(blocks))])
        try:
            sample = data.subset(indices)
            tree, _ = infer(
                sample, recipe_id, Path(output) / f"replicate-{i + 1}", seed=seed + i + 1, **kwargs
            )
            trees.append(newick(tree))
            for split in split_set(tree)[1]:
                counts[split] = counts.get(split, 0) + 1
        except ValueError as exc:
            failures.append({"replicate": i + 1, "error": str(exc)})
    result = copy.deepcopy(reference)
    names = split_set(result)[0]
    for node in result.get_nonterminals():
        side = frozenset(t.name for t in node.get_terminals())
        other = names - side
        if min(len(side), len(other)) > 1 and trees:
            split = min(tuple(sorted(side)), tuple(sorted(other)))
            node.confidence = 100 * counts.get(split, 0) / len(trees)
    return (
        result,
        {
            "requested": replicates,
            "successful": len(trees),
            "failures": failures,
            "block_size": block_size,
            "seed": seed,
            "support_denominator": "successful replicates; review failures before interpreting support",
            "resampling": "sites" if block_size == 1 else "contiguous blocks in input order",
        },
        trees,
    )
