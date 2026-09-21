# Implementation validation — 2026-09-21

This is software validation, not a publication benchmark or a claim of biological novelty.

- 60 pytest cases passed on Python 3.12, including the original 17 tests.
- All eight ML presets were exercised against the real IQ-TREE 2.4.0 Linux executable, rather than a mocked engine.
- The bundled benchmark completed all 18 runs: 17 presets plus a labelled uncorrected SNP control. It uses small synthetic inputs and does not establish superiority over other methods.
- Tests cover exact distance formulae, saturation, unequal-cluster UPGMA weighting, missing observations, VCF.gz, heterozygote/phase policies, methylation count validation, ASC compatibility, RF comparison, known-tree recovery, seeded bootstrap and output preservation.
- The local HTTP interface was tested for input/result round trips, ZIP contents and rejection of cross-origin requests.
- Chromium browser checks exercised DNA NJ, DNA HKY maximum likelihood, and beta UPGMA; no JavaScript page errors were observed. Desktop (1440 px) and mobile (390 px) views were inspected. A mobile overflow was corrected and verified absent. Result downloads were available.
- Ruff, JavaScript syntax validation, and git whitespace checks passed.

Runtime libraries used during validation: NumPy 2.3.5, SciPy 1.17.0 and Biopython 1.88. CI runs Python 3.10 and 3.12 and installs a checksum-verified IQ-TREE 2.4.0 binary. Library versions supported by each interpreter are resolved during CI; each run records its actual versions.

Reproduce from a checkout with IQ-TREE on PATH:

```bash
pip install -e '.[dev]'
ruff check src tests
pytest
node --check src/methylphylo/web/app.js
methylphylo-research benchmark examples/research/benchmark-all.json --output results/validation
```

Real datasets, memory/scaling benchmarks, model misspecification beyond the bundled simulators, and publication-scale uncertainty analyses remain study work. The UI is a local Python application. Updating this repository does not deploy that backend to the original static public website.
