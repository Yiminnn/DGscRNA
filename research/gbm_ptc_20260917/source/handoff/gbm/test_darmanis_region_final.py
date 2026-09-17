"""分region终端注释的最小契约测试；仅在SLURM执行。"""
import os
import unittest

import anndata as ad
import numpy as np
import pandas as pd

from darmanis_region_final import oracle_bounds, refine_final, select_unit


class RegionContract(unittest.TestCase):
    def data(self, labels, clusters=None):
        a = ad.AnnData(np.ones((len(labels), 3), dtype=np.float32))
        a.obs_names = [f'cell{i}' for i in range(len(labels))]
        a.obs['cluster'] = clusters or ['0'] * len(labels)
        a.obs['initial'] = labels
        a.obsm['X_scaled'] = a.X.copy()
        return a

    def test_absent_class_oracle(self):
        bounds = oracle_bounds(['Neuron', 'Neuron'], {'Neuron', 'Myeloid'})
        self.assertEqual(bounds['oracle_present_upper_bound'], 1)
        self.assertEqual(bounds['oracle_fixed7_sample_upper_bound'], 1 / 7)
        self.assertEqual(bounds['oracle_fixed7_vocab_upper_bound'], 2 / 7)

    def test_noise_excluded_and_unknown_is_pool(self):
        a = self.data(['Unknown', 'Unknown', 'A', 'B'], ['Noise', '0', '1', '2'])
        def train(sub, key, **kwargs):
            self.assertNotIn('cell0', sub.obs_names)
            self.assertEqual(sub.obs[key].tolist(), ['Undecided', 'A', 'B'])
            self.assertFalse(kwargs['use_highly_variable'])
            np.testing.assert_array_equal(sub.obsm['X_scaled'], a.X[1:])
            return {'num_classes': 2, 'input_dim': 3}
        def predict(model, sub, key, **kwargs):
            self.assertEqual(kwargs['probability_threshold'], .9)
            return pd.Series(['A', 'A', 'B'], index=sub.obs_names)
        final, lineage, status = refine_final(a, 'initial', train, predict)
        self.assertEqual(final.index.tolist(), a.obs_names.tolist())
        self.assertEqual(final.tolist(), ['Unknown', 'A', 'A', 'B'])
        self.assertEqual(lineage.tolist()[0:2], ['noise_terminal', 'dl_confident'])
        self.assertEqual(status['n_training_classes'], 2)
        self.assertEqual(status['n_pool'], 1)
        self.assertTrue(status['training_executed'])
        self.assertEqual(status['dl_status'], 'dl_completed')

    def test_noop_and_fallback(self):
        def forbidden(*args, **kwargs):
            self.fail('本分支不应训练')
        for labels, expected in [(['A', 'B'], 'empty_pool'),
                                 (['A', 'Unknown'], 'lt2_training_classes')]:
            _, _, status = refine_final(self.data(labels), 'initial', forbidden)
            self.assertEqual(status['dl_status'], expected)
            self.assertFalse(status['training_executed'])
            self.assertEqual(status['final_valid'], expected == 'empty_pool')
        def fail(*args, **kwargs):
            raise ValueError('模拟训练失败')
        final, _, status = refine_final(self.data(['A', 'B', 'Unknown']), 'initial', fail)
        self.assertEqual(status['dl_status'], 'dl_error')
        self.assertFalse(status['final_valid'])
        self.assertEqual(final.iloc[-1], 'Unknown')

    def test_unit_isolation(self):
        a = self.data(['A'] * 4)
        a.obs['donor_id'] = ['BT_S1', 'BT_S1', 'BT_S2', 'BT_S2']
        a.obs['Location'] = ['Tumor', 'Periphery', 'Tumor', 'Periphery']
        left = select_unit(a, 'BT_S1', 'Tumor')
        right = select_unit(a, 'BT_S2', 'Periphery')
        self.assertEqual(left.obs_names.tolist(), ['cell0'])
        self.assertEqual(right.obs_names.tolist(), ['cell3'])
        self.assertFalse(set(left.obs_names) & set(right.obs_names))
        with self.assertRaises(ValueError):
            select_unit(a, 'BT_S6', 'Distant')


if __name__ == '__main__':
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('测试只允许在SLURM作业内运行')
    unittest.main()
