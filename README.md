# MethylPhylo Studio

[![CI](https://github.com/Sanjeevan-informatik/methylation-to-phylogenetic-tree/actions/workflows/ci.yml/badge.svg)](https://github.com/Sanjeevan-informatik/methylation-to-phylogenetic-tree/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776ab)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-6ee7d8)](LICENSE)

MethylPhylo Studio v0.2 turns a DNA-methylation beta-value matrix into a transparent exploratory analysis: input validation, locus-level imputation, variance-based feature selection, pairwise sample distances, a rooted UPGMA tree, PCA, evolutionary-model probabilities, and reproducible exports.

**[Open live application](https://methylphylo-studio.sanjeevanvive.chatgpt.site)** — select **Model probability** to explore the JC69 and two-state models.

> This is exploratory research software, not a clinical diagnostic tool. A similarity tree based on methylation does not by itself establish cellular ancestry.

## Version 0.2 highlights

- **Jukes–Cantor (JC69)** nucleotide substitution probabilities.
- A simplified **methylated ↔ unmethylated** symmetric two-state model.
- Same-state, change-state, and individual substitution probabilities.
- Reusable Python functions with numerical and scientific validation.
- JC69 correction from an observed nucleotide mismatch proportion.
- **17 automated tests** covering the analysis pipeline and probability models.

For JC69 at an evolutionary distance of `d = 0.10` substitutions per site:

| Outcome | Probability |
| --- | ---: |
| Same nucleotide | **90.64%** |
| Any different nucleotide | **9.36%** |
| Each particular alternative nucleotide | **3.12%** |

## Why this project

This project connects my MSc research in SNP-based phylogenetic inference with current biomedical data engineering and methylation analysis. I built it as a focused return-to-bioinformatics project in 2026, with an emphasis on the skills expected in production-facing research software: a tested scientific core, reproducible interfaces, clear assumptions, privacy-aware interaction, automation, and concise documentation.

## What it demonstrates

- **Bioinformatics reasoning:** beta-value QC, missing-data handling, unsupervised feature selection, distance metrics, UPGMA, PCA, Newick, and continuous-time Markov models.
- **Software engineering:** typed Python models, a reusable package, deterministic tests, errors at data boundaries, CLI outputs, and separation of concerns.
- **User-centred delivery:** an accessible browser workspace with CSV/TSV input, interactive plots, quality reporting, and local-only processing.
- **Reproducibility and DevOps:** Docker, pinned runtime requirements, GitHub Actions, JSON provenance, and an example dataset.

## Try it locally

### Browser application

The browser app has no build step. Serve the repository root and open `/dist/`:

```bash
python -m http.server 8000
```

Then visit `http://localhost:8000/dist/`. Choose the synthetic cohort or upload your own matrix.

### Python command line

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .

methylphylo examples/demo_methylation.tsv \
  --top-loci 10 \
  --metric euclidean \
  --output-dir methylphylo-results
```

The CLI writes:

- `tree.nwk` — rooted UPGMA tree with branch lengths;
- `report.json` — parameters, quality metrics, sample names, the distance matrix, and Newick.

The Docker interface uses the same command:

```bash
docker build -t methylphylo .
docker run --rm -v "$PWD:/data" methylphylo \
  /data/examples/demo_methylation.tsv --top-loci 10 --output-dir /data/results
```

## Input format

CSV and TSV are accepted. The first column contains locus identifiers; each remaining column is one sample. Beta values must be in the closed interval `[0, 1]`. Empty cells and `NA`, `NaN`, `null`, or `.` are treated as missing.

```text
cpg_id      Control_1  Control_2  Primary_1
cg10200000  0.10       0.12       0.62
cg10200137  0.20       NA         0.58
```

## Method

1. Validate the matrix shape, unique sample names, numeric types, and beta-value range.
2. Drop loci that are entirely missing and report them.
3. Median-impute missing values independently within each locus.
4. Rank loci by population variance and retain the requested number.
5. Compute RMS Euclidean, mean Manhattan, or correlation distance between samples.
6. Build a rooted average-linkage UPGMA dendrogram and serialize it as Newick.
7. In the browser, calculate a two-component sample-space PCA for an orthogonal view of cohort structure.

UPGMA assumes an ultrametric structure. The output is therefore best interpreted as an exploratory similarity tree. For biological inference, compare metrics, inspect confounders and batch effects, and validate clusters with an appropriate domain model.

## Transition probabilities

The **Model probability** tab provides two forward probability calculators:

- **Jukes–Cantor (JC69)** for nucleotide substitutions. At evolutionary distance `d` (expected substitutions per site), `P(same) = 1/4 + 3/4 exp(-4d/3)`. Each named substitution has probability `1/4 - 1/4 exp(-4d/3)`; the probability of any different nucleotide is three times that value.
- **Symmetric two-state model** for a simplified methylated ↔ unmethylated process. If `d = rate × time`, then `P(same) = 1/2 + 1/2 exp(-2d)` and `P(switch) = 1/2 - 1/2 exp(-2d)`.

The JC69 correction for an observed nucleotide mismatch proportion `p` is `d = -3/4 ln(1 - 4p/3)`, defined for `0 ≤ p < 0.75`.

```python
from methylphylo import (
    binary_state_probabilities,
    jukes_cantor_distance,
    jukes_cantor_probabilities,
)

jc69 = jukes_cantor_probabilities(distance=0.1)
corrected_distance = jukes_cantor_distance(p_distance=0.1)
binary = binary_state_probabilities(distance=0.1)
```

These model distances have defined evolutionary units. The application intentionally does not treat a beta-value tree distance as a JC69 substitution distance.

## Quality checks

```bash
pip install -e ".[dev]"
ruff check src tests
pytest
node --check dist/app.js
```

The tests cover missing-value imputation, distance-matrix invariants, expected within-group similarity, Newick completeness, invalid beta values, deterministic reruns, JC69 inversion, probability normalization, model limits, and invalid model distances. CI executes on Python 3.10 and 3.12.

## Repository map

```text
dist/                    browser application
src/methylphylo/         reusable Python package and CLI
tests/                   scientific and reproducibility tests
examples/                synthetic input matrix
docs/architecture.md     data flow and design trade-offs
.github/workflows/       continuous integration
```

See [the architecture note](docs/architecture.md) for component boundaries, trade-offs, and a production-scale extension path.

## Responsible use

- The bundled dataset is synthetic and contains no patient data.
- Browser uploads are processed in memory and never leave the device.
- Exported reports contain sample labels; treat them according to the source dataset's governance rules.
- This software is for education, portfolio demonstration, and exploratory research—not clinical decision-making.

## Author

**Sanjeevan Vivekanantha** — MSc Bioinformatics; Python, phylogenetics, RNA-seq, single-cell data systems, and reproducible research software.

I am currently looking for bioinformatics, computational biology, scientific software, and biomedical data roles in Germany.
