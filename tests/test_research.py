"""Scientific invariants, input failure modes, and complete research workflows."""

import gzip
import io
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from Bio import Phylo

from methylphylo.datasets import (
    Dataset,
    load_dataset,
    prepare,
    read_counts,
    read_fasta,
    read_matrix_text,
    read_vcf,
)
from methylphylo.inference import (
    RECIPES,
    bootstrap,
    compare_trees,
    distances,
    infer,
    newick,
    validate_recipe,
)
from methylphylo.research import benchmark, run_analysis, simulate

DNA = ">A\nAACCAA\n>B\nAACCAT\n>C\nAATTAA\n>D\nAATTAT\n"
VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD
1\t10\t.\tA\tG\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t1/1
1\t20\t.\tC\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1
"""


def test_distance_formula_and_overlap():
    data = read_fasta(">A\nAAAA\n>B\nGAAA\n>C\nCAAA\n>D\nAAAN\n")
    p, overlap = distances(data, "p")
    jc, _ = distances(data, "jc69")
    k2p, _ = distances(data, "k2p")
    assert p[0, 1] == 0.25
    assert jc[0, 1] == pytest.approx(-0.75 * np.log(2 / 3))
    assert k2p[0, 1] == pytest.approx(-0.5 * np.log(0.5))
    assert k2p[0, 2] == pytest.approx(-0.5 * np.log(0.75) - 0.25 * np.log(0.5))
    assert overlap[0, 3] == 3
    assert np.allclose(jc, jc.T)


@pytest.mark.parametrize("model", ["jc69", "k2p"])
def test_saturation_is_not_clipped(model):
    data = read_fasta(">A\nAAAA\n>B\nGGGG\n>C\nCCCC\n>D\nTTTT\n")
    with pytest.raises(ValueError, match="saturation"):
        distances(data, model)


@pytest.mark.parametrize(
    "text",
    [
        ">A\nA\n>B\nAA\n>C\nAA\n>D\nAA\n",
        ">A\nAA\n>A\nAA\n>C\nAA\n>D\nAA\n",
        ">A\nAR\n>B\nAA\n>C\nAA\n>D\nAA\n",
        ">A\nNN\n>B\nAA\n>C\nAA\n>D\nAA\n",
        ">A bad\nAA\n>B\nAA\n>C\nAA\n>D\nAA\n",
    ],
)
def test_invalid_alignment_rejected(text):
    with pytest.raises(ValueError):
        read_fasta(text)


def test_vcf_never_imputes_reference_and_gzip(tmp_path):
    data = read_vcf(VCF)
    assert data.loci == ("1:10", "1:20")
    assert data.values[1, 0] == -1
    path = tmp_path / "snp.vcf.gz"
    with gzip.open(path, "wt") as f:
        f.write(VCF)
    loaded = load_dataset(path, "vcf")
    np.testing.assert_array_equal(data.values, loaded.values)
    assert loaded.provenance["input_sha256"]
    assert prepare(data).values.shape[1] == 1


def test_phasing_is_never_invented():
    with pytest.raises(ValueError, match="phased"):
        read_vcf(VCF, "phased")
    data = read_vcf(VCF.replace("0/1", "0|1"), "phased")
    assert len(data.names) == 8
    assert data.values[2, 0] == 0 and data.values[3, 0] == 2


def test_vcf_bad_allele_and_duplicate_site():
    with pytest.raises(ValueError, match="out of range"):
        read_vcf(VCF.replace("0/1", "0/8"))
    with pytest.raises(ValueError, match="unique"):
        read_vcf(VCF.replace("1\t20", "1\t10"))


def test_counts_depth_and_duplicate_validation():
    text = "locus\tsample\tmethylated\ttotal\n" + "".join(
        f"{loc}\t{s}\t{m}\t{n}\n" for loc, m, n in [("cg1", 2, 10), ("cg2", 0, 0)] for s in "ABCD"
    )
    data = read_counts(text)
    assert np.all(data.values[:, 0] == 0.2)
    assert np.all(np.isnan(data.values[:, 1]))
    assert data.provenance["low_depth_observations"] == 4
    with pytest.raises(ValueError, match="Duplicate"):
        read_counts(text + "cg1\tA\t2\t10\n")
    with pytest.raises(ValueError, match="Counts must"):
        read_counts(text.replace("2\t10", "11\t10"))


def test_thresholds_and_missing_are_explicit():
    data = read_matrix_text("locus,A,B,C,D\ncg1,0.1,0.2,0.8,0.9\ncg2,0.1,0.5,0.9,0.9", "beta")
    binary = prepare(data, low=0.2, high=0.8)
    assert binary.kind == "binary" and len(binary.loci) == 1
    np.testing.assert_array_equal(binary.values[:, 0], [0, 0, 1, 1])
    retained = prepare(data, missing="retain", low=0.2, high=0.8)
    assert retained.values[1, 1] == -1
    with pytest.raises(ValueError):
        prepare(data, low=0.8, high=0.2)


def test_binary_values_cannot_be_continuous():
    with pytest.raises(ValueError, match="Invalid values"):
        read_matrix_text("locus,A,B,C,D\nx,0,0.5,1,1", "binary")


def test_recipe_and_asc_compatibility():
    dna = read_fasta(DNA)
    with pytest.raises(ValueError, match="needs beta"):
        validate_recipe(dna, "beta-upgma")
    with pytest.raises(ValueError, match="variable-only"):
        validate_recipe(dna, "dna-ml-gtr-asc")
    snps = prepare(dna, sites="variable")
    validate_recipe(snps, "dna-ml-gtr-asc")
    with pytest.raises(ValueError, match="ASC recipe"):
        validate_recipe(snps, "dna-ml-gtr")
    validate_recipe(snps, "dna-ml-gtr", allow_uncorrected=True)
    missing = prepare(
        read_fasta(DNA.replace("AACCAA", "AANCAA")), missing="retain", sites="variable"
    )
    with pytest.raises(ValueError, match="complete sites"):
        validate_recipe(missing, "dna-ml-gtr-asc")


def test_rf_root_invariance_and_tip_mismatch():
    a = Phylo.read(io.StringIO("((A,B),(C,D));"), "newick")
    b = Phylo.read(io.StringIO("(A,B,(C,D));"), "newick")
    c = Phylo.read(io.StringIO("((A,C),(B,D));"), "newick")
    assert compare_trees(a, b)["rf"] == 0
    assert compare_trees(a, c)["normalized_rf"] == 1
    with pytest.raises(ValueError, match="same tip"):
        compare_trees(a, Phylo.read(io.StringIO("((A,B),(C,E));"), "newick"))


def test_simulation_reproducible_and_true_tree_recovered(tmp_path):
    data = simulate(tmp_path / "sim1", sites=5000, seed=19)
    again = simulate(tmp_path / "sim2", sites=5000, seed=19)
    np.testing.assert_array_equal(data.values, again.values)
    report = run_analysis(
        data,
        "dna-nj-jc69",
        tmp_path / "run",
        replicates=10,
        seed=3,
        truth=tmp_path / "sim1/truth.nwk",
    )
    assert report["reference_comparison"]["rf"] == 0
    assert report["bootstrap"]["successful"] == 10
    assert (tmp_path / "run/tree.svg").exists()
    assert report["data"]["sha256"]
    with pytest.raises(ValueError, match="not empty"):
        run_analysis(data, "dna-nj-p", tmp_path / "run")


@pytest.mark.parametrize("recipe", [k for k, v in RECIPES.items() if v.method != "ml"])
def test_builtin_recipes_end_to_end(tmp_path, recipe):
    kind = RECIPES[recipe].kind
    if kind == "dna":
        data = read_fasta(DNA)
    else:
        array = np.array([[0, 0, 1, 0], [0, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 1]])
        data = Dataset(tuple("ABCD"), tuple("wxyz"), array, kind)
    tree, report = infer(data, recipe, tmp_path)
    assert {t.name for t in tree.get_terminals()} == set("ABCD")
    assert report["model"] == RECIPES[recipe].model
    if RECIPES[recipe].method == "parsimony":
        assert report["parsimony_score"] >= 0
        assert all(c.branch_length is None for c in tree.find_clades())
        assert ":" not in newick(tree)


def test_block_bootstrap_reproducible(tmp_path):
    data = read_fasta(DNA)
    tree, _ = infer(data, "dna-nj-p", tmp_path)
    a, ra, _ = bootstrap(data, "dna-nj-p", tree, tmp_path / "a", 12, seed=9, block_size=2)
    b, rb, _ = bootstrap(data, "dna-nj-p", tree, tmp_path / "b", 12, seed=9, block_size=2)
    assert newick(a) == newick(b)
    assert ra == rb and ra["successful"] == 12


def test_parsimony_missing_not_treated_as_extra_character(tmp_path):
    data = read_fasta(DNA.replace("AACCAA", "AACCAN"))
    with pytest.raises(ValueError, match="complete sites"):
        infer(data, "dna-parsimony", tmp_path)


def test_benchmark_preserves_failures_and_data_identity(tmp_path):
    (tmp_path / "input.fa").write_text(DNA)
    manifest = {
        "runs": [
            {"id": "nj", "input": "input.fa", "data_type": "dna", "recipe": "dna-nj-p"},
            {"id": "jc", "input": "input.fa", "data_type": "dna", "recipe": "dna-nj-jc69"},
            {"id": "wrong", "input": "input.fa", "data_type": "dna", "recipe": "beta-upgma"},
        ]
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    rows = benchmark(tmp_path / "manifest.json", tmp_path / "out")
    assert [r["status"] for r in rows] == ["ok", "ok", "failed"]
    assert rows[0]["data_sha256"] == rows[1]["data_sha256"]
    comparisons = json.loads((tmp_path / "out/comparisons.json").read_text())
    assert comparisons[0]["same_observations"]


def test_cli_error_is_actionable(tmp_path):
    path = tmp_path / "invalid.fa"
    path.write_text(DNA)
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "methylphylo.research_cli",
            "infer",
            str(path),
            "--data-type",
            "dna",
            "--recipe",
            "beta-upgma",
            "--output",
            str(tmp_path / "out"),
        ],
        capture_output=True,
        check=False,
        text=True,
    )
    assert run.returncode == 2
    assert "needs beta data" in run.stderr


IQTREE = shutil.which("iqtree2") or shutil.which("iqtree")


@pytest.mark.skipif(not IQTREE, reason="IQ-TREE executable not installed")
@pytest.mark.parametrize("recipe", [k for k, v in RECIPES.items() if v.method == "ml"])
def test_real_iqtree_integration(tmp_path, recipe):
    data = simulate(tmp_path / "sim", kind=RECIPES[recipe].kind, sites=1000, seed=8, gain=1, loss=1)
    if "asc" in recipe:
        data = prepare(data, sites="variable")
    report = run_analysis(data, recipe, tmp_path / "run", executable=IQTREE, timeout=90)
    assert "IQ-TREE" in report["engine_version"]
    assert Path(report["engine_report"]).exists()
    assert report["rooting"].startswith("unrooted")


def test_upgma_uses_sample_weights_for_unequal_clusters(tmp_path):
    data = Dataset(
        tuple("ABCDE"), ("site",), np.array([[0.0], [0.02], [0.1], [0.4], [0.9]]), "beta"
    )
    tree, _ = infer(data, "beta-upgma", tmp_path)
    # Final split: ABCD vs E. Mean intercluster distance = .9 - mean(0,.02,.1,.4).
    expected_height = (0.9 - np.mean([0, 0.02, 0.1, 0.4])) / 2
    assert all(tree.distance(tip) == pytest.approx(expected_height) for tip in tree.get_terminals())
