"""Table4纯评价函数；不读真实数据、不拟合模型、不删细胞，数值调用仅限SLURM。

主表Pair-F1追溯原notebook；NPE明确为TV，KS另报。保留方差与距离RV分离。
所有公共函数返回可直接JSON序列化的数值、列表或字典；不可评返回None及reason。
"""
import hashlib
import json
import math

import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import (adjusted_rand_score, fowlkes_mallows_score,
                             pairwise_distances, v_measure_score)
from sklearn.metrics.cluster import contingency_matrix


K_NEIGHBORS = 20


def _truth(labels):
    labels = np.asarray(labels, dtype=object)
    if labels.ndim != 1 or not labels.size:
        raise ValueError('真值必须为非空一维标签，禁止删除无效细胞')
    if any(not isinstance(v, str) or not v.strip() for v in labels):
        raise ValueError('真值必须是非空字符串；缺失/NaN标签须在上游明确处理')
    return labels


def _clusters(labels, n_cells):
    labels = np.asarray(labels)
    if labels.ndim != 1 or labels.size != n_cells:
        raise ValueError('预测标签长度必须与全部真值细胞一致')
    if not np.issubdtype(labels.dtype, np.number) or np.issubdtype(labels.dtype, np.complexfloating):
        raise ValueError('cluster ID必须为整数，-1仅表示noise')
    if (not np.isfinite(labels).all() or not np.equal(labels, np.rint(labels)).all()
            or (labels < -1).any() or (labels >= 2 ** 63).any()):
        raise ValueError('cluster ID必须为-1或非负int64整数')
    return labels.astype(np.int64, copy=False)


def _matrix(x):
    if sparse.issparse(x):
        x = sparse.csr_matrix(x)
        values = x.data
    else:
        x = np.asarray(x)
        values = x
    if x.ndim != 2 or min(x.shape) < 1:
        raise ValueError('评价矩阵必须为非空二维矩阵')
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise ValueError('评价矩阵必须为实数')
    if not np.isfinite(values).all():
        raise ValueError('评价矩阵含非finite值，禁止静默删除')
    return x


def _positive_int(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f'{name}必须为正整数')
    return int(value)


def partition_metrics(y_true, clusters):
    """全体细胞的Pair-F1/ARI/FMI/V；noise作为独立partition类别，不参与映射。"""
    truth = _truth(y_true)
    pred = _clusters(clusters, len(truth))
    n = len(truth)
    result = {'n_cells': n, 'n_clusters': int(np.unique(pred[pred != -1]).size),
              'noise_rate': float(np.mean(pred == -1)), 'pair_f1': None,
              'ari': None, 'fmi': None, 'v_measure': None,
              'pair_f1_convention': 'reference_pair_count', 'undefined_reasons': {}}
    if n < 2:
        result['undefined_reasons'] = {key: 'fewer_than_two_cells'
                                      for key in ('pair_f1', 'ari', 'fmi', 'v_measure')}
        return result
    c = contingency_matrix(truth, pred, sparse=True)
    # Python整数防止pair计数平方在大样本中溢出；稀疏列联表不创建N×N数组。
    tk = sum(int(v) * (int(v) - 1) for v in c.data)
    pk = sum(int(v) * (int(v) - 1) for v in np.asarray(c.sum(axis=0)).ravel())
    qk = sum(int(v) * (int(v) - 1) for v in np.asarray(c.sum(axis=1)).ravel())
    if pk + qk == 0:
        result['pair_f1'] = 1.0
        result['pair_f1_convention'] = 'both_singleton_partitions'
    else:
        result['pair_f1'] = float(2 * tk / (pk + qk))
    result.update(ari=float(adjusted_rand_score(truth, pred)),
                  fmi=float(fowlkes_mallows_score(truth, pred)),
                  v_measure=float(v_measure_score(truth, pred)))
    return result


def macro_f1_22(y_true, clusters, fine_labels, other_label='Other'):
    """事后Hungarian-overlap对齐的固定22类macro-F1；不是原论文Pair-F1。"""
    truth = _truth(y_true)
    pred = _clusters(clusters, len(truth))
    fine = list(_truth(fine_labels))
    if len(fine) != 22 or len(set(fine)) != 22 or other_label in fine:
        raise ValueError('fine_labels必须是固定22个互异标签，且不含Other')
    if not isinstance(other_label, str) or not other_label.strip():
        raise ValueError('Other标签必须是非空字符串')
    vocabulary = fine + [other_label]
    positions = {label: i for i, label in enumerate(vocabulary)}
    if any(label not in positions for label in truth):
        raise ValueError('真值包含22类加Other以外的标签，禁止过滤')
    active = pred != -1
    unique_clusters, columns = np.unique(pred[active], return_inverse=True)
    table = np.zeros((23, len(unique_clusters)), dtype=np.int64)
    rows = np.array([positions[label] for label in truth[active]], dtype=np.int64)
    np.add.at(table, (rows, columns), 1)
    mapping = []
    if unique_clusters.size:
        # 行序固定为22类名单+Other；列为排序后的cluster ID，noise从构表前排除。
        row_idx, col_idx = linear_sum_assignment(table, maximize=True)
        mapping = [{'cluster': int(unique_clusters[j]), 'label': vocabulary[i],
                    'overlap': int(table[i, j])} for i, j in zip(row_idx, col_idx)]
        mapping.sort(key=lambda row: row['cluster'])
    lookup = {row['cluster']: row['label'] for row in mapping}
    aligned = np.array([lookup.get(int(cluster)) for cluster in pred], dtype=object)
    per_class = []
    for label in fine:
        true_mask, pred_mask = truth == label, aligned == label
        tp = int(np.count_nonzero(true_mask & pred_mask))
        fp = int(np.count_nonzero(~true_mask & pred_mask))
        fn = int(np.count_nonzero(true_mask & ~pred_mask))
        denominator = 2 * tp + fp + fn
        per_class.append({'label': label, 'support': int(true_mask.sum()),
                          'tp': tp, 'fp': fp, 'fn': fn,
                          'f1': float(2 * tp / denominator) if denominator else 0.0})
    return {'macro_f1_22': float(sum(row['f1'] for row in per_class) / 22),
            'n_cells': len(truth), 'n_classes_averaged': 22,
            'n_noise_cells': int(np.count_nonzero(~active)),
            'n_unmatched_cells': int(np.count_nonzero(active & (aligned == None))),
            'mapping': mapping, 'per_class': per_class, 'undefined_reasons': {},
            'definition': 'posthoc_hungarian_overlap_fixed22_all_cells'}


def _fixed_k(k):
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or k != K_NEIGHBORS:
        raise ValueError('此契约固定k=20，禁止静默改变邻居数')


def same_label_neighbor_counts(x, labels, k=20, chunk_size=128):
    """分块精确Euclidean邻居；等距按原cell index稳定排序，显式排除自身。"""
    _fixed_k(k)
    chunk_size = _positive_int(chunk_size, 'chunk_size')
    x, truth = _matrix(x), _truth(labels)
    n = len(truth)
    if x.shape[0] != n:
        raise ValueError('矩阵行数和全部真值标签长度不一致')
    digest = hashlib.sha256(json.dumps(truth.tolist(), ensure_ascii=False,
                                      separators=(',', ':')).encode('utf-8')).hexdigest()
    result = {'counts': None, 'k': int(k), 'n_cells': n,
              'labels_sha256': digest, 'reason': None,
              'distance': 'euclidean', 'tie_break': 'ascending_cell_index'}
    if n <= k:
        result['reason'] = 'n_cells_not_greater_than_k'
        return result
    counts = np.empty(n, dtype=np.int64)
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        # 最大中间距离数组为chunk_size×N，而不是N×N×features。
        distances = pairwise_distances(x[start:stop], x, metric='euclidean', n_jobs=1)
        if not np.isfinite(distances).all():
            result['reason'] = 'nonfinite_pairwise_distances'
            return result
        distances[np.arange(stop - start), np.arange(start, stop)] = np.inf
        nearest = np.argsort(distances, axis=1, kind='stable')[:, :k]
        counts[start:stop] = np.sum(truth[nearest] == truth[start:stop, None], axis=1)
    result['counts'] = counts.tolist()
    return result


def _neighbor_counts(counts, n, k):
    counts = np.asarray(counts)
    if counts.ndim != 1 or len(counts) != n:
        raise ValueError('邻域同类计数长度必须与全部标签一致')
    if not np.issubdtype(counts.dtype, np.number) or np.iscomplexobj(counts):
        raise ValueError('邻域同类计数必须为整数')
    if (not np.isfinite(counts).all() or not np.equal(counts, np.rint(counts)).all()
            or (counts < 0).any() or (counts > k).any()):
        raise ValueError('邻域同类计数必须为0..20中的整数')
    return counts.astype(np.int64, copy=False)


def npe_from_counts(original_counts, embedded_counts, labels, k=20):
    """复用原X的同类邻居计数；按实际出现类别(含Other)等权计算TV及KS。"""
    _fixed_k(k)
    truth = _truth(labels)
    classes = np.unique(truth)
    result = {'npe_tv': None, 'npe_ks': None, 'per_class': [],
              'n_cells': len(truth), 'n_present_classes': len(classes), 'k': int(k),
              'reason': None, 'averaging': 'equal_weight_present_true_classes_including_Other'}
    if original_counts is None or embedded_counts is None:
        result['reason'] = 'neighbor_counts_unavailable'
        return result
    x = _neighbor_counts(original_counts, len(truth), k)
    y = _neighbor_counts(embedded_counts, len(truth), k)
    if len(truth) <= k:
        result['reason'] = 'n_cells_not_greater_than_k'
        return result
    for label in classes:
        mask = truth == label
        support = int(mask.sum())
        p = np.bincount(x[mask], minlength=k + 1) / support
        q = np.bincount(y[mask], minlength=k + 1) / support
        result['per_class'].append({'label': label, 'support': support,
                                    'tv': float(0.5 * np.abs(p - q).sum()),
                                    'ks': float(np.max(np.abs(np.cumsum(p) - np.cumsum(q)))),
                                    'pmf_x': p.tolist(), 'pmf_y': q.tolist()})
    result['npe_tv'] = float(np.mean([row['tv'] for row in result['per_class']]))
    result['npe_ks'] = float(np.mean([row['ks'] for row in result['per_class']]))
    return result


def sample_index_pairs(n_cells, max_pairs=100000, seed=42):
    """均匀无放回抽取unordered pairs；Floyd抽样仅用O(max_pairs)存储。"""
    n = _positive_int(n_cells, 'n_cells')
    limit = _positive_int(max_pairs, 'max_pairs')
    total = n * (n - 1) // 2
    if total > np.iinfo(np.int64).max:
        raise ValueError('pair总数超过int64随机抽样支持范围')
    count = min(limit, total)
    rng = np.random.default_rng(seed)
    if count == total:
        ranks = range(total)
    else:
        selected = set()
        for j in range(total - count, total):
            candidate = int(rng.integers(j + 1))
            selected.add(j if candidate in selected else candidate)
        ranks = sorted(selected)
    pairs = []
    for rank in ranks:
        # 将上三角lexicographic rank逆映射为(i,j)，整数平方根避免大N浮点偏差。
        i = (2 * n - 1 - math.isqrt((2 * n - 1) ** 2 - 8 * rank)) // 2
        start = i * (2 * n - i - 1) // 2
        if start > rank:
            i -= 1
            start = i * (2 * n - i - 1) // 2
        pairs.append([int(i), int(i + 1 + rank - start)])
    return pairs


def _pairs(pairs, n):
    if len(pairs) == 0:
        return np.empty((0, 2), dtype=np.int64)
    pairs = np.asarray(pairs)
    if pairs.ndim != 2 or pairs.shape[1] != 2 or not np.issubdtype(pairs.dtype, np.integer):
        raise ValueError('pairs必须是M×2整数indices')
    if (pairs[:, 0] < 0).any() or (pairs[:, 0] >= pairs[:, 1]).any() or (pairs[:, 1] >= n).any():
        raise ValueError('pairs必须满足0<=i<j<N，无对角线')
    if len(set(map(tuple, pairs))) != len(pairs):
        raise ValueError('pairs不可重复，必须无放回')
    return pairs.astype(np.int64, copy=False)


def pair_distances(x, pairs, batch_size=1024):
    """仅计算指定pair的Euclidean距离，按batch处理feature差值，支持CSR输入。"""
    x = _matrix(x)
    pairs = _pairs(pairs, x.shape[0])
    batch_size = _positive_int(batch_size, 'batch_size')
    distances = np.empty(len(pairs), dtype=np.float64)
    for start in range(0, len(pairs), batch_size):
        batch = pairs[start:start + batch_size]
        difference = x[batch[:, 0]].astype(np.float64) - x[batch[:, 1]].astype(np.float64)
        if sparse.issparse(difference):
            squared = np.asarray(difference.multiply(difference).sum(axis=1)).ravel()
        else:
            squared = np.einsum('ij,ij->i', difference, difference)
        distances[start:start + len(batch)] = np.sqrt(squared)
    if not np.isfinite(distances).all():
        raise ValueError('计算pair距离产生非finite值，禁止过滤无效pairs')
    return distances.tolist()


def distance_residual_variance(d_x, d_y):
    """新增统一Euclidean距离补项1-Pearson²；输入pair顺序必须由调用者保持一致。"""
    x, y = np.asarray(d_x, dtype=np.float64), np.asarray(d_y, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError('d_x和d_y必须是一一对应、同长度的pair距离向量')
    result = {'rv_e': None, 'n_pairs': len(x), 'reason': None,
              'definition': '1_minus_squared_pearson_euclidean_pair_distances'}
    if len(x) < 2:
        result['reason'] = 'fewer_than_two_pairs'
    elif not np.isfinite(x).all() or not np.isfinite(y).all():
        result['reason'] = 'nonfinite_distances'
    elif (x < 0).any() or (y < 0).any():
        result['reason'] = 'negative_distances'
    elif np.ptp(x) == 0 or np.ptp(y) == 0:
        result['reason'] = 'zero_distance_variance'
    else:
        # Pearson对正比例缩放不变；先缩放，避免距离平方的上溢或下溢。
        scaled_x, scaled_y = x / np.max(x), y / np.max(y)
        centered_x = scaled_x - np.mean(scaled_x)
        centered_y = scaled_y - np.mean(scaled_y)
        correlation = np.dot(centered_x, centered_y) / (np.linalg.norm(centered_x) * np.linalg.norm(centered_y))
        if not np.isfinite(correlation):
            result['reason'] = 'nonfinite_correlation'
        else:
            correlation = float(np.clip(correlation, -1, 1))
            result['rv_e'] = float(max(0.0, 1 - correlation ** 2))
    return result


def retained_variance(method, metadata):
    """仅从fit元数据读PCA/FA方差；不推断其他方法、不替换成距离RV。"""
    result = {'retained_variance_pct': None, 'definition': None, 'reason': None}
    if method not in ('PCA', 'FA'):
        result['reason'] = 'not_defined_for_method'
        return result
    if method == 'PCA':
        result['definition'] = 'pca_explained_variance_ratio_percent'
        if 'explained_variance_ratio' not in metadata:
            result['reason'] = 'missing_fit_metadata'
            return result
        ratios = np.asarray(metadata['explained_variance_ratio'], dtype=float)
        if (ratios.ndim != 1 or not ratios.size or not np.isfinite(ratios).all()
                or (ratios < 0).any() or ratios.sum() > 1 + 1e-8):
            result['reason'] = 'invalid_fit_variance'
            return result
        result['retained_variance_pct'] = float(100 * ratios.sum())
    else:
        result['definition'] = 'fa_shared_loading_variance_fraction_percent'
        if any(key not in metadata for key in ('fa_signal_variance', 'fa_noise_variance')):
            result['reason'] = 'missing_fit_metadata'
            return result
        try:
            signal, noise = float(metadata['fa_signal_variance']), float(metadata['fa_noise_variance'])
        except (TypeError, ValueError):
            result['reason'] = 'invalid_fit_variance'
            return result
        if not np.isfinite([signal, noise]).all() or signal < 0 or noise < 0:
            result['reason'] = 'invalid_fit_variance'
        elif signal == noise == 0:
            result['reason'] = 'zero_total_model_variance'
        else:
            scale = max(signal, noise)
            result['retained_variance_pct'] = float(100 * (signal / scale) / (signal / scale + noise / scale))
    return result


def visualization_metadata(method):
    """原文未给ED/DD自动判据；未人工判读时不猜测。"""
    if method in ('none', 'NO-reduction'):
        return {'value': 'N/A', 'reason': 'no_2d_embedding'}
    return {'value': 'Not assessed', 'reason': 'requires_visual_assessment'}
