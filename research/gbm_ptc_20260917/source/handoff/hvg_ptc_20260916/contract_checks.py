#!/usr/bin/env python3
"""Scientific contracts tested in a SLURM allocation, before full fit submission."""
from common import ROOT, OUT, require_slurm, write_json, utc
require_slurm()
import sys
sys.path[:0] = [str(ROOT/'handoff/g274_table4'), str(ROOT/'handoff/gbm')]
import unittest
import tempfile
from pathlib import Path
import anndata as ad
import numpy as np
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score
from test_darmanis_region_final import RegionContract
from metrics_core import partition_metrics
from fit import annotate, representation, clusters
from darmanis_region_final import refine_final


class HVGContract(unittest.TestCase):
    def test_hdbscan_configuration_serializes(self):
        from common import clean
        import json
        x = np.random.default_rng(3).normal(size=(80, 2)).astype(np.float32)
        arm = dict(clusterer='HDBSCAN', min_cluster_size=15, min_samples=15, families=[])
        cl, info = clusters(x, arm, 42)
        json.dumps(clean(info), allow_nan=False)
        self.assertEqual(info['params']['memory']['location'], None)
        self.assertEqual(len(cl), 80)

    def test_noise_terminal(self):
        a = ad.AnnData(sp.csr_matrix(np.ones((30, 9), dtype=np.float32)))
        a.obsm['X_scaled'] = np.ones(a.shape, dtype=np.float32)
        arm = {'scoring_features': 'all', 'cutoff': 'none'}
        with tempfile.TemporaryDirectory() as d:
            m = annotate(a, np.full(30, -1), np.arange(3), arm, 42, Path(d))
            p = np.load(Path(d)/'predictions.npz', allow_pickle=False)
            self.assertTrue(m['terminal_valid'])
            self.assertFalse(m['final_valid'])
            self.assertFalse(m['training_executed'])
            self.assertTrue((p['final'] == 'Unknown').all())
            self.assertTrue((p['lineage'] == 'noise_terminal').all())

    def test_geometry_really_subsets_and_no_hidden_pca(self):
        x = np.random.default_rng(17).normal(size=(60, 50)).astype(np.float32)
        idx = np.array([0, 3, 5, 7, 11])
        g = dict(dr='none', dim=None, input_space='genes', seed=42, neighbors=15, min_dist=.1)
        z, meta = representation(x[:, idx], g)
        np.testing.assert_array_equal(z, x[:, idx])
        self.assertEqual(meta['geometry_feature_width'], 5)
        self.assertNotIn('pre_pca', meta)
        g.update(dr='PCA', dim=2)
        z, meta = representation(x[:, idx], g)
        self.assertEqual(z.shape, (60, 2))
        self.assertEqual(meta['reducer_input_shape'], [60, 5])

    def test_pair_metric_independent_pair_count(self):
        truth = np.array(['A', 'A', 'B', 'B', 'C'])
        cluster = np.array([0, 0, 1, -1, -1])
        pairs = [(i, j) for i in range(5) for j in range(i)]
        tp = sum(truth[i] == truth[j] and cluster[i] == cluster[j] for i,j in pairs)
        goldpos = sum(truth[i] == truth[j] for i,j in pairs)
        predpos = sum(cluster[i] == cluster[j] for i,j in pairs)
        m = partition_metrics(truth, cluster)
        self.assertAlmostEqual(m['pair_f1'], 2*tp/(goldpos+predpos))
        self.assertAlmostEqual(m['ari'], adjusted_rand_score(truth, cluster))

    def test_real_dl_execution_uses_full_width(self):
        # Exercise the real classifier with two learnable seed classes plus an unlabelled pool.
        import torch
        torch.set_num_threads(2)
        n, p = 130, 12
        rng = np.random.default_rng(5)
        x = rng.normal(0, .1, (n, p)).astype(np.float32)
        x[:60, :6] += 2
        x[60:120, 6:] += 2
        a = ad.AnnData(sp.csr_matrix((n, p), dtype=np.float32))
        a.obs['cluster'] = ['0']*60 + ['1']*60 + ['2']*9 + ['Noise']
        a.obs['seed'] = ['A']*60 + ['B']*60 + ['Undecided']*10
        a.var['highly_variable'] = [True]*3 + [False]*9
        a.obsm['X_scaled'] = x
        final, lineage, state = refine_final(a, 'seed')
        self.assertEqual(state['dl_status'], 'dl_completed', state)
        self.assertTrue(state['training_executed'] and state['final_valid'])
        self.assertEqual(final.iloc[-1], 'Unknown')
        self.assertEqual(state['n_training'], 120)
        self.assertEqual(state['n_pool'], 9)
        self.assertEqual(final.iloc[:120].tolist(), ['A']*60 + ['B']*60)


if __name__ == '__main__':
    loader = unittest.TestLoader()
    suite = unittest.TestSuite([loader.loadTestsFromTestCase(c) for c in [RegionContract, HVGContract]])
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    write_json(OUT/'verification/contracts.json', dict(timestamp=utc(), tests=result.testsRun,
        failures=len(result.failures), errors=len(result.errors), passed=result.wasSuccessful()))
    sys.exit(0 if result.wasSuccessful() else 1)
