"""Strict, explicit data representations for comparative phylogenetics."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

SAFE_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")
DNA = set("ACGT")
MISSING = {"", "na", "nan", "null", ".", "?", "-", "n"}


@dataclass
class Dataset:
    names: tuple[str, ...]
    loci: tuple[str, ...]
    values: np.ndarray  # samples x sites; -1 for missing categorical states, NaN for beta
    kind: str  # dna, binary, beta
    provenance: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def __post_init__(self):
        if self.kind not in {"dna", "binary", "beta"}:
            raise ValueError("Unknown dataset kind")
        if len(self.names) < 4 or len(set(self.names)) != len(self.names):
            raise ValueError("Research inference requires at least four uniquely named samples")
        if any(not SAFE_NAME.fullmatch(name) for name in self.names):
            raise ValueError("Sample names must contain only letters, digits, _, . or -")
        if not self.loci or len(set(self.loci)) != len(self.loci):
            raise ValueError("Site identifiers must be nonempty and unique")
        if self.values.shape != (len(self.names), len(self.loci)):
            raise ValueError("Data shape does not match sample/site identifiers")
        if self.kind == "beta":
            valid = np.isnan(self.values) | (
                np.isfinite(self.values) & (self.values >= 0) & (self.values <= 1)
            )
        else:
            valid = np.isin(self.values, [-1, 0, 1, 2, 3] if self.kind == "dna" else [-1, 0, 1])
        if not np.all(valid):
            raise ValueError(f"Invalid values for {self.kind} data")
        if np.any(np.all(self.missing, axis=1)):
            raise ValueError("A sample has no observed sites")

    @property
    def missing(self):
        return np.isnan(self.values) if self.kind == "beta" else self.values < 0

    def subset(self, indices):
        # Bootstrap resamples may repeat sites; give them unique identifiers.
        indices = np.asarray(indices, dtype=int)
        return Dataset(
            self.names,
            tuple(f"{i}:{self.loci[j]}" for i, j in enumerate(indices)),
            self.values[:, indices].copy(),
            self.kind,
            dict(self.provenance),
            list(self.warnings),
        )

    def fingerprint(self):
        h = hashlib.sha256()
        h.update(self.kind.encode())
        h.update("\n".join(self.names).encode())
        h.update("\n".join(self.loci).encode())
        h.update(np.asarray(self.values, dtype="<f8").tobytes())
        return h.hexdigest()


def read_text(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        return handle.read()


def read_fasta(text, kind="dna"):
    names, seqs = [], []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            names.append(line[1:].strip())
            seqs.append("")
        elif not seqs:
            raise ValueError("FASTA must begin with a >sample header")
        else:
            seqs[-1] += line.upper()
    if not seqs or len(set(map(len, seqs))) != 1 or not seqs[0]:
        raise ValueError("FASTA must be a nonempty, equal-length alignment")
    alphabet = "ACGT" if kind == "dna" else "01"
    allowed = set(alphabet + "N?-")
    if any(set(s) - allowed for s in seqs):
        raise ValueError("Unsupported FASTA states; resolve IUPAC ambiguity explicitly or use N")
    values = np.array([[alphabet.find(c) for c in s] for s in seqs], dtype=int)
    return Dataset(tuple(names), tuple(str(i + 1) for i in range(len(seqs[0]))), values, kind)


def read_matrix_text(text, kind):
    first = text.splitlines()[0] if text.splitlines() else ""
    rows = list(csv.reader(io.StringIO(text), delimiter="\t" if "\t" in first else ","))
    if len(rows) < 2:
        raise ValueError("Matrix requires a header and site rows")
    names = tuple(x.strip() for x in rows[0][1:])
    loci, values = [], []
    for line, row in enumerate(rows[1:], 2):
        if not row or not any(x.strip() for x in row):
            continue
        if len(row) != len(names) + 1 or not row[0].strip():
            raise ValueError(f"Malformed matrix row {line}")
        loci.append(row[0].strip())
        values.append([np.nan if x.strip().lower() in MISSING else float(x) for x in row[1:]])
    array = np.array(values, dtype=float).T
    if kind == "binary":
        array = np.where(np.isnan(array), -1, array)
    return Dataset(names, tuple(loci), array, kind)


def read_counts(text, min_depth=5):
    """Long TSV: locus, sample, methylated, total. Never invent zero-depth observations."""
    if min_depth < 1:
        raise ValueError("min_depth must be positive")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not {"locus", "sample", "methylated", "total"} <= set(reader.fieldnames or []):
        raise ValueError("Counts TSV needs locus, sample, methylated, total columns")
    names, loci, records = [], [], {}
    filtered = 0
    for row in reader:
        sample, locus = row["sample"], row["locus"]
        if not locus:
            raise ValueError("Empty locus in counts")
        key = (sample, locus)
        if key in records:
            raise ValueError(f"Duplicate count for {sample}, {locus}")
        m, n = int(row["methylated"]), int(row["total"])
        if not 0 <= m <= n:
            raise ValueError("Counts must satisfy 0 <= methylated <= total")
        if sample not in names:
            names.append(sample)
        if locus not in loci:
            loci.append(locus)
        records[key] = m / n if n >= min_depth else np.nan
        filtered += n < min_depth
    values = np.array([[records.get((s, loc), np.nan) for loc in loci] for s in names])
    return Dataset(
        tuple(names),
        tuple(loci),
        values,
        "beta",
        {"min_depth": min_depth, "low_depth_observations": filtered},
        ["Counts converted to beta values; downstream distance models do not model depth."],
    )


def read_vcf(text, genotype_policy="homozygous"):
    """SNP-only conversion; never fill unlisted reference sites or invent phase."""
    if genotype_policy not in {"homozygous", "phased"}:
        raise ValueError("VCF policy must be homozygous or phased")
    names, loci, columns = None, [], []
    skipped, masked = 0, 0
    for line in text.splitlines():
        if line.startswith("#CHROM\t"):
            names = tuple(line.split("\t")[9:])
            continue
        if line.startswith("#") or not line.strip():
            continue
        if names is None:
            raise ValueError("VCF requires a #CHROM header")
        f = line.split("\t")
        if len(f) != 9 + len(names):
            raise ValueError("VCF sample column count mismatch")
        alleles = [f[3]] + f[4].split(",")
        if f[6] not in {"PASS", "."} or any(a not in DNA for a in alleles):
            skipped += 1
            continue
        if int(f[1]) < 1:
            raise ValueError("VCF positions must be positive")
        fields = f[8].split(":")
        if "GT" not in fields:
            raise ValueError("VCF record has no GT field")
        gt_index = fields.index("GT")
        states = []
        for cell in f[9:]:
            parts = cell.split(":")
            gt = parts[gt_index] if gt_index < len(parts) else "."
            tokens = re.split(r"[/|]", gt)
            if any(t == "." for t in tokens):
                states.extend([-1, -1] if genotype_policy == "phased" else [-1])
                continue
            indices = [int(t) for t in tokens]
            if any(i < 0 or i >= len(alleles) for i in indices):
                raise ValueError("VCF genotype allele index out of range")
            bases = ["ACGT".index(alleles[i]) for i in indices]
            if genotype_policy == "phased":
                if len(bases) != 2 or ("|" not in gt and bases[0] != bases[1]):
                    raise ValueError("Phased mode requires diploid phased GTs (or homozygous GTs)")
                states.extend(bases)
            elif len(set(bases)) == 1:
                states.append(bases[0])
            else:
                states.append(-1)
                masked += 1
        loci.append(f"{f[0]}:{f[1]}")
        columns.append(states)
    if names is None or not columns:
        raise ValueError("No usable SNP records in VCF")
    warnings = [
        "VCF is treated as SNP-only; unlisted sites are not assumed callable reference.",
        "Site filters (e.g. MAF) may introduce bias beyond simple +ASC correction.",
    ]
    if genotype_policy == "phased":
        names = tuple(f"{name}_h{h}" for name in names for h in (1, 2))
        warnings.append("Phase continuity is not validated; restrict to a consistent phase block.")
    else:
        warnings.append(f"Heterozygous genotypes masked as missing: {masked}.")
    return Dataset(
        names,
        tuple(loci),
        np.array(columns).T,
        "dna",
        {"snp_only": True, "genotype_policy": genotype_policy, "skipped_vcf_records": skipped},
        warnings,
    )


def load_dataset(path, data_type, **options):
    path = Path(path)
    text = read_text(path)
    if data_type in {"dna", "binary-fasta"}:
        data = read_fasta(text, "dna" if data_type == "dna" else "binary")
    elif data_type in {"beta", "binary"}:
        data = read_matrix_text(text, data_type)
    elif data_type == "counts":
        data = read_counts(text, options.get("min_depth", 5))
    elif data_type == "vcf":
        data = read_vcf(text, options.get("genotype_policy", "homozygous"))
    else:
        raise ValueError("Unknown input data type")
    data.provenance.update(
        {
            "input_type": data_type,
            "input_path": str(path.resolve()),
            "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    )
    return data


def prepare(data, missing="complete", sites="all", low=None, high=None):
    if missing not in {"complete", "retain"} or sites not in {"all", "variable"}:
        raise ValueError("Invalid missing-data or site-selection policy")
    before = data.values.shape[1]
    if low is not None or high is not None:
        if data.kind != "beta" or low is None or high is None or not 0 <= low < high <= 1:
            raise ValueError("Binarisation needs beta data and 0 <= low < high <= 1")
        arr = np.full(data.values.shape, -1, dtype=int)
        arr[data.values <= low] = 0
        arr[data.values >= high] = 1
        data = Dataset(
            data.names,
            data.loci,
            arr,
            "binary",
            dict(data.provenance),
            data.warnings + ["Thresholded beta values are not direct cell-level states."],
        )
        data.provenance.update({"binary_low": low, "binary_high": high})
    keep = ~np.all(data.missing, axis=0)
    if missing == "complete":
        keep &= ~np.any(data.missing, axis=0)
    if sites == "variable":
        keep &= np.array(
            [len(set(col[~miss])) > 1 for col, miss in zip(data.values.T, data.missing.T)]
        )
    if not keep.any():
        raise ValueError("No sites remain after filtering; change the explicit missing/site policy")
    result = Dataset(
        data.names,
        tuple(np.array(data.loci)[keep]),
        data.values[:, keep].copy(),
        data.kind,
        dict(data.provenance),
        list(data.warnings),
    )
    result.provenance.update(
        {
            "input_sites": before,
            "retained_sites": int(keep.sum()),
            "missing_policy": missing,
            "site_selection": sites,
        }
    )
    if missing == "retain":
        result.warnings.append(
            "Missing states retained; distances use pairwise deletion, ML marginalises."
        )
    if sites == "variable":
        result.warnings.append("Variable-site selection changes the sampling distribution.")
    return result


def write_fasta(data, path):
    if data.kind == "beta":
        raise ValueError("Continuous beta values cannot be written as a character alignment")
    alphabet = "ACGT" if data.kind == "dna" else "01"
    Path(path).write_text(
        "".join(
            f">{name}\n" + "".join(alphabet[int(x)] if x >= 0 else "?" for x in row) + "\n"
            for name, row in zip(data.names, data.values)
        ),
        encoding="utf-8",
    )
