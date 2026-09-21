#!/usr/bin/env python3
"""Preserve A1 source/protocol bytes; adapt only IO and R library isolation."""
from pathlib import Path
import hashlib
import json
import shutil

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
OLD = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
SOURCE = OLD / 'source_snapshots/embedding_v8'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    original = CODE / 'original'
    protocol = CODE / 'protocol'
    original.mkdir(exist_ok=True)
    protocol.mkdir(exist_ok=True)
    manifest = json.loads((SOURCE / 'SOURCE_MANIFEST.json').read_text())
    names = ['run.py', 'adaptive_ica.py', 'export_geometry.R', 'native_geometry.R', 'score_candidates.R']
    originals = {}
    for name in names:
        assert sha(SOURCE / name) == manifest[name]
        shutil.copy2(SOURCE / name, original / name)
        originals[str((SOURCE / name).relative_to(ROOT))] = sha(SOURCE / name)
    shutil.copy2(SOURCE / 'SOURCE_MANIFEST.json', original / 'SOURCE_MANIFEST.json')
    for name in ['embedding.json', 'embedding_convergence_repair_20260920_v2.json']:
        shutil.copy2(OLD / 'protocol' / name, protocol / name)
        originals[str((OLD / 'protocol' / name).relative_to(ROOT))] = sha(OLD / 'protocol' / name)
    markers = ROOT / 'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/markers/libraries.json'
    shutil.copy2(markers, protocol / 'libraries.json')
    originals[str(markers.relative_to(ROOT))] = sha(markers)
    shutil.copy2(original / 'adaptive_ica.py', CODE / 'adaptive_ica.py')
    isolation = ".libPaths(.Library,include.site=FALSE)"
    replacements = {
        'export_geometry.R': [
            ("suppressPackageStartupMessages(library(Seurat))", isolation + "\nsuppressPackageStartupMessages(library(Seurat))"),
            ("cells<-read.csv(file.path(prep,'cells.csv'),stringsAsFactors=FALSE)",
             "cells<-read.csv(file.path(prep,'cells.csv'),stringsAsFactors=FALSE,colClasses='character',na.strings=NULL)"),
        ],
        'native_geometry.R': [
            ("root<-'/fs/scratch/PCON0080/yimin/dgscrna'\n.libPaths(c(file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments/vendor_R_dbscan64'),.libPaths()))", isolation),
            ("cells<-read.csv(file.path(cfg$geometry,'cells.csv'),stringsAsFactors=FALSE)$cell_id",
             "cells<-read.csv(file.path(cfg$geometry,'cells.csv'),stringsAsFactors=FALSE,colClasses='character',na.strings=NULL)$cell_id"),
            ("e<-read.csv(cfg$embedding,stringsAsFactors=FALSE)",
             "e<-read.csv(cfg$embedding,stringsAsFactors=FALSE,colClasses=c('character','numeric','numeric'),na.strings=NULL)"),
        ],
        'score_candidates.R': [
            ("root <- '/fs/scratch/PCON0080/yimin/dgscrna'\nbase <- file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')\n.libPaths(c(file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments/vendor_R_dbscan64'),.libPaths()))", isolation),
            ("cells <- read.csv(file.path(prep,'cells.csv'),stringsAsFactors=FALSE)",
             "cells <- read.csv(file.path(prep,'cells.csv'),stringsAsFactors=FALSE,colClasses='character',na.strings=NULL)"),
            ("marker_path <- file.path(base,'markers','libraries.json')", "marker_path <- cfg$markers"),
        ],
    }
    for name, changes in replacements.items():
        text = (original / name).read_text()
        for before, after in changes:
            assert text.count(before) == 1, (name, before)
            text = text.replace(before, after)
        (CODE / name).write_text(text)
        # Every difference is one of the declared non-numerical substitutions.
        restored = text
        for before, after in reversed(changes):
            assert restored.count(after) == 1
            restored = restored.replace(after, before)
        assert restored == (original / name).read_text()
    (CODE / 'SOURCE_PROTOCOL_MANIFEST.json').write_text(json.dumps(dict(
        status='source_adaptation_frozen_before_pilot', original_sources=originals,
        preserved_python_functions=['get_geometry', 'model_metadata', 'fit_model', 'make_embedding', 'make_partitions'],
        R_substitutions={name: [{'before': before, 'after': after} for before, after in changes]
                         for name, changes in replacements.items()},
        adapted_files={name: sha(CODE / name) for name in replacements} |
                      {'adaptive_ica.py': sha(CODE / 'adaptive_ica.py')},
        frozen_protocols={p.name: sha(p) for p in protocol.iterdir() if p.is_file()},
        numerical_body_policy='Python AST functions are executed unchanged from original/run.py; R changes reverse exactly to original bytes',
        DL_policy='Invoke installed dgscrna.reference.backend.terminal; fresh output-local cache only',
        output_policy='Separate package campaign, old result reads allowed only in verification',
        scope=dict(samples=121, budgets=2, spaces=7, representation_units=1694,
                   partitions_per_unit=13, terminal_conditions=22022),
    ), indent=2) + '\n')


if __name__ == '__main__':
    main()
