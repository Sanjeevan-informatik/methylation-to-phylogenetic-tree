"""Regenerate bundled synthetic fixtures; use a fresh temporary directory."""
import shutil
import tempfile
from pathlib import Path

from methylphylo.research import simulate

base = Path(__file__).parent
with tempfile.TemporaryDirectory() as directory:
    for kind in ('dna', 'binary'):
        folder = Path(directory) / kind
        simulate(folder, kind=kind, sites=1000, seed=8, gain=1, loss=1)
        shutil.copyfile(folder / 'alignment.fasta', base / f'{kind}.fasta')
        shutil.copyfile(folder / 'truth.nwk', base / f'{kind}-truth.nwk')
        shutil.copyfile(folder / 'simulation.json', base / f'{kind}-simulation.json')
