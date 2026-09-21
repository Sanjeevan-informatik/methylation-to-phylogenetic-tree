# MethylPhylo research workbench 0.3

This release implements a reproducible comparison framework for DNA, SNP and methylation data. It supports the experimental work needed to extend *Efficient Phylogenetic Inference Using SNP-Based Approaches: A Comparison with Full Sequence Data*. Software implementation is not evidence of a new biological result or a guarantee of publication.

## Install and launch

```bash
python -m venv .venv
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -e '.[dev]'
methylphylo-research serve
```

Open http://127.0.0.1:8765. The local interface accepts text files and runs Python on your own machine. It is different from the original static beta-value explorer in `dist/`. It binds only to loopback. Uploaded input and run files are temporary: download the ZIP before closing the page. The ZIP includes the original upload, retained observations, site mapping, Newick, topology SVG, reports and engine logs. Browser runs are limited to 12 MB request bodies, 150 samples, 100,000 sites, 25 bootstrap replicates and 120 seconds per inference. Bootstrap is a series of inferences and may take longer overall. The CLI supports gzip and larger datasets, subject to available RAM; readers currently load files into memory.

ML recipes require a separately installed [IQ-TREE](https://iqtree.github.io/doc/Quickstart) executable on PATH (`iqtree2`, `iqtree` or `iqtree3`), or `--iqtree /absolute/path/to/iqtree2` on CLI runs. Integration was checked with IQ-TREE 2.4.0. The program records the detected version and exact command; no ML result is fabricated if the engine is unavailable. NJ and parsimony use Biopython. UPGMA uses the existing sample-weighted cluster-average implementation, including when cluster sizes differ. Simulation uses SciPy matrix exponentiation.

The original command remains compatible:

```bash
methylphylo examples/demo_methylation.tsv --output-dir legacy-results
```

## Inputs

All research datasets need at least four samples. Labels must be unique and contain only letters, digits, underscores, periods and hyphens. Site identifiers must also be unique. Genome assemblies, coordinates, tissue, batch, assay and sample matching must be checked by the investigator before loading.

| `--data-type` | Representation | Policy |
|---|---|---|
| `dna` | Aligned FASTA, A/C/G/T with N/?/- missing | Equal-length aligned homologous sites; IUPAC ambiguity rejected explicitly |
| `vcf` | VCF or VCF.gz with sample GT fields | SNP records only; PASS or unfiltered records; indels, symbolic alleles and failed records skipped and counted |
| `beta` | CSV/TSV: first column locus, remaining columns samples | Finite observed values in [0,1]; blank/NA/NaN/null/./?/-/N missing |
| `counts` | Long TSV: `locus`, `sample`, `methylated`, `total` | Nonnegative integer counts; `0 <= methylated <= total`; below `--min-depth` (default 5) becomes missing |
| `binary` | CSV/TSV matrix | Only 0/1 or missing |
| `binary-fasta` | Aligned FASTA | Only 0/1, with N/?/- missing |

Counts are converted to beta values. This is **not a coverage-aware evolutionary likelihood model**. Binary thresholds can be applied explicitly to beta or counts input with `--low 0.2 --high 0.8`; intermediate values become missing. Bulk beta values can represent cell/allele mixtures, so thresholding does not establish a cellular methylation state.

VCF default `--genotype-policy homozygous` retains haploid/homozygous calls and masks heterozygotes. `--genotype-policy phased` produces two haplotypes per diploid sample and rejects unphased heterozygotes. It does not infer phase or check continuity across phase sets. Restrict input to a consistent phase block before using phased mode. QUAL, depth, genotype quality, allele balance and MAF thresholds are not applied by this converter: perform and document appropriate upstream QC. No unlisted VCF position is filled from a reference. Reconstructing full callable sequence requires separate callability evidence.

## Presets

`methylphylo-research catalog` provides a machine-readable catalog. These are 17 combinations of data, algorithm and model, not 17 novel algorithms.

| Recipe | Input | Algorithm/model |
|---|---|---|
| beta-upgma | beta | UPGMA, RMS Euclidean distance |
| beta-nj | beta | NJ, RMS Euclidean distance |
| beta-nj-manhattan | beta | NJ, mean Manhattan distance |
| dna-nj-p | DNA | NJ, observed mismatch proportion |
| dna-nj-jc69 | DNA | NJ, JC69-corrected distance |
| dna-nj-k2p | DNA | NJ, Kimura two-parameter distance |
| dna-parsimony | DNA | Fitch parsimony, NNI heuristic from NJ start |
| dna-ml-jc69 | DNA | IQ-TREE ML, JC69 |
| dna-ml-hky | DNA | IQ-TREE ML, HKY+G4 |
| dna-ml-gtr | DNA | IQ-TREE ML, GTR+G4 |
| dna-ml-gtr-asc | variable-only DNA | IQ-TREE ML, GTR+G4+ASC |
| binary-nj | binary | NJ, Hamming mismatch proportion |
| binary-parsimony | binary | Fitch parsimony, NNI heuristic from NJ start |
| binary-ml-jc2 | binary | IQ-TREE ML, symmetric JC2 |
| binary-ml-gtr2 | binary | IQ-TREE ML, asymmetric GTR2 |
| binary-ml-jc2-asc | variable-only binary | IQ-TREE ML, JC2+ASC |
| binary-ml-gtr2-asc | variable-only binary | IQ-TREE ML, GTR2+ASC |

JC69 and K2P distances reject saturation instead of clipping it. Negative NJ branches are retained and flagged. Parsimony is heuristic, not a guaranteed global optimum, and has no estimated evolutionary branch lengths. NJ and reversible ML trees are unrooted; the Newick root is representational. UPGMA produces an ultrametric rooted similarity tree, not evidence of a biological root. The SVG uses equal edge spacing and is labelled as a topology preview, not a branch-length plot.

## Missingness and ascertainment

- `--missing complete` (default) drops any site with missing observations. No median imputation or high-variance feature selection is performed in the research pipeline.
- `--missing retain` keeps partially observed sites. Distances use pairwise deletion and report pairwise overlap counts. ML treats missing states as unknown; parsimony requires complete observations.
- All-missing sites are removed, and all-missing samples are rejected.
- `--sites variable` explicitly removes constant sites. `retained-sites.tsv` maps every retained column back to its source site.
- ASC recipes require complete variable-only observations in this release; ambiguous constant patterns are rejected.
- Variable-only ML without ASC is rejected unless `--allow-uncorrected` is supplied for an explicitly labelled bias control. This override is CLI-only.
- Simple ASC correction conditions on variable sites. It does not undo MAF filtering, array ascertainment, selected CpGs, or arbitrary upstream filters.

Raw likelihoods/AIC/BIC across different observation sets must not be used to rank methods. Engine likelihood reports are retained, but the summary deliberately does not rank different datasets by likelihood.

## A runnable SNP versus full-sequence experiment

```bash
# Bundled fixtures are synthetic, with known true trees.
methylphylo-research benchmark examples/research/benchmark-builtin.json --output results/builtin

# All 17 presets plus a labelled uncorrected SNP control; requires IQ-TREE.
methylphylo-research benchmark examples/research/benchmark-all.json --output results/all-models
```

The manifest compares full DNA, variable-only DNA with ASC, and variable-only DNA without ASC. All DNA versions derive from the same simulated alignment. Beta and binary examples are separate data representations; their different units and assumptions must not be conflated. Output paths must be fresh: the program refuses to overwrite prior results.

Each manifest run specifies `id`, `input`, `data_type`, `recipe` and optional `seed`, `threads`, `timeout`, `missing`, `sites`, `truth`, `replicates`, `block_size`, `min_depth`, `genotype_policy`, `low`, `high`, `allow_uncorrected`. Input and truth paths are relative to the original manifest. The copied manifest preserves the original specification; retain source data and paths when rerunning it elsewhere.

Outputs:

- `summary.csv` / `summary.json`: success/failure, retained-data hash, site count, inference seconds, optional RF to reference.
- `comparisons.json`: pairwise unrooted RF for trees with the same tip labels, and whether their observation hashes match.
- Per-run `report.json`, `tree.nwk`, `tree.svg`, retained observations and site mapping.
- ML runs retain IQ-TREE alignment, log and report; bootstrap runs retain replicate trees and engine files.

Failed runs remain in the summary and make the benchmark command exit nonzero. Comparisons against an empirical reference indicate agreement, not accuracy. Only an independently justified truth supports accuracy claims. RF is the symmetric difference of unrooted nontrivial splits; normalized RF divides by the sum of observed splits. Unresolved trees require careful interpretation. No branch-length accuracy or quartet metric is currently implemented.

## Simulations and uncertainty

```bash
methylphylo-research simulate --kind dna --sites 5000 --seed 19 --output results/sim-dna
methylphylo-research infer results/sim-dna/alignment.fasta --data-type dna \
  --recipe dna-nj-jc69 --bootstrap 100 --seed 3 \
  --truth results/sim-dna/truth.nwk --output results/dna-nj

methylphylo-research simulate --kind binary --gain 0.2 --loss 0.1 \
  --missing 0.1 --error 0.01 --seed 9 --output results/sim-binary
methylphylo-research infer results/sim-binary/alignment.fasta --data-type binary-fasta \
  --missing retain --recipe binary-ml-gtr2 --output results/binary-ml

methylphylo-research infer examples/research/counts.tsv --data-type counts \
  --low 0.2 --high 0.8 --recipe binary-nj --output results/counts-binary
```

Simulations use an explicit known tree (optional `--tree truth.nwk`), a stationary root and independent sites. DNA uses JC69 with branch lengths in expected substitutions/site. Binary uses a gain/loss CTMC with branches measured in the corresponding time units. Seeds, Q matrix and noise settings are saved. Simulations do not model recombination, LD, mixtures or selection. Do not claim general performance based only on these simple generating models.

Ordinary bootstrap resamples retained sites. `--block-size N` samples contiguous blocks in **input order**, preserving each sampled block. The final short block may produce variable replicate lengths. Blocks are not automatically separated by chromosome; pre-split chromosomes/windows. This resamples already filtered observations, so uncertainty is conditional on preprocessing. Failed replicates are listed and support uses successful replicates; failures can bias support and must be reviewed. A full paper should rerun relevant preprocessing within resamples and use justified genomic blocks when sites are dependent.

## Study design and remaining work

A defensible study should use multiple simulation settings/replicate seeds, real callable sequence data, external cohorts, and independent lineage evidence for methylation. Benchmark accuracy and speed with matched samples/sites and identical hardware. Inference wall time is measured; process peak memory is **not** measured. Record CPU, RAM, OS and engine build separately. Do not interpret unrelated patients as a cellular ancestry tree.

Matched DNA and methylation trees can be compared using `compare` or a manifest with matching sample IDs. Joint genetic–epigenetic likelihood, Bayesian inference, protein/codon models, coverage-aware methylation likelihood, ancestral-state reconstruction, controlled-access data download and public job hosting are not implemented. Do not treat an unsupported model as implemented because a related preset exists.

Before a new paper: choose and document real dataset accessions, prepare callability/QC scripts, evaluate alternative generating models and established baselines, freeze dependencies and parameters, archive outputs and code with a DOI, and write results from actual experiments. This release is the research infrastructure for that work.

## References

- [Original SNP/full-sequence preprint](https://doi.org/10.1101/2025.03.10.642383)
- [IQ-TREE models and ASC](https://iqtree.github.io/doc/Substitution-Models)
- [IQ-TREE command reference](https://iqtree.github.io/doc/Command-Reference)
- [Biopython tree construction](https://biopython.org/docs/latest/api/Bio.Phylo.TreeConstruction.html)
- [MethylTree baseline](https://github.com/ShouWenWang-Lab/MethylTree)
- [EPI-Clone study](https://doi.org/10.1038/s41586-025-09041-8)
