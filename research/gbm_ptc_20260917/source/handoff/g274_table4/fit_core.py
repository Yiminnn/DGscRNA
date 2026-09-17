"""Table4的无标签数值核心；导入不读写文件，不访问gold或历史partition。"""
from dataclasses import dataclass
import time
import warnings

import numpy as np
from scipy import sparse

DR_ORDER = ('PCA', 'FA', 'ICA', 'Isomap', 'UMAP', 'TSNE', 'none')
CLUSTER_ORDER = ('KMeans', 'GMM', 'HDBSCAN')


@dataclass(frozen=True)
class FitConfig:
    gmm_covariance: str
    seed: int = 42
    n_components: int = 2
    k: int = 23
    min_cluster_size: int = 15
    min_samples: int = 15
    n_neighbors: int = 15
    target_sum: float = 10000.0
    max_value: float = 10.0

    def __post_init__(self):
        if self.gmm_covariance not in ('full', 'diag', 'full_2d_diag_none'):
            raise ValueError('GMM协方差必须明确指定full、diag或full_2d_diag_none')
        if self.n_components != 2 or self.k != 23:
            raise ValueError('冻结主实验为2维和预先指定K23，不按gold类别数重算')
        if self.min_cluster_size != 15 or self.min_samples != 15:
            raise ValueError('冻结HDBSCAN min_cluster_size=min_samples=15')
        if self.n_neighbors < 2 or self.target_sum <= 0 or self.max_value <= 0:
            raise ValueError('邻居数、归一化目标和截断值必须有效')


def prepare_matrix(counts, target_sum=10000.0, max_value=10.0):
    """完整输入→normalize_total→log1p→scale；不做HVG或额外细胞/基因筛选。"""
    import anndata as ad
    import scanpy as sc

    if len(counts.shape) != 2 or min(counts.shape) < 2:
        raise ValueError('至少需要两个细胞和两个基因')
    data = counts.data if sparse.issparse(counts) else np.asarray(counts)
    if not np.isfinite(data).all() or (data < 0).any():
        raise ValueError('counts必须非负且finite')
    if not np.equal(data, np.rint(data)).all():
        raise ValueError('输入必须为原始整数counts，不能重复归一化')
    if target_sum <= 0 or max_value <= 0:
        raise ValueError('归一化目标与截断值必须为正数')
    x = sparse.csr_matrix(counts, dtype=np.float32, copy=True)
    x.sum_duplicates()
    x.eliminate_zeros()
    if (np.asarray(x.sum(axis=1)).ravel() <= 0).any():
        raise ValueError('存在全零细胞，禁止静默删除')
    a = ad.AnnData(x)
    sc.pp.normalize_total(a, target_sum=target_sum)
    sc.pp.log1p(a)
    # 显式稠密化，避免依赖稀疏scale版本间的隐式行为及提示。
    a.X = a.X.toarray()
    sc.pp.scale(a, max_value=max_value)
    result = np.asarray(a.X, dtype=np.float32)
    if not np.isfinite(result).all():
        raise ValueError('归一化/缩放后出现非finite值')
    return result


def _validate_matrix(x):
    x = np.asarray(x)
    if x.ndim != 2 or min(x.shape) < 2 or not np.isfinite(x).all():
        raise ValueError('拟合输入必须为至少2×2的finite矩阵')
    return x


def _fit_metadata(model, caught):
    result = {'params': model.get_params(deep=False),
              'warnings': [f'{w.category.__name__}: {w.message}' for w in caught]}
    for key in ('n_iter_', 'converged_', 'kl_divergence_', 'reconstruction_error_'):
        value = getattr(model, key, None)
        if value is not None and not callable(value):
            result[key.rstrip('_')] = value.item() if isinstance(value, np.generic) else value
    return result


def reduce_matrix(x, method, config):
    """每个reducer直接接收同一X；none保留所有列，绝不经过PCA/HVG。"""
    from sklearn.decomposition import PCA, FactorAnalysis, FastICA
    from sklearn.manifold import Isomap, TSNE

    x = _validate_matrix(x)
    if method == 'none':
        return x, {'params': {'method': 'none'}, 'warnings': []}
    if method not in DR_ORDER:
        raise ValueError(f'未知降维方法：{method}')
    kwargs = {'n_components': config.n_components}
    if method == 'PCA':
        model = PCA(**kwargs, random_state=config.seed)
    elif method == 'FA':
        model = FactorAnalysis(**kwargs, max_iter=5000, random_state=config.seed)
    elif method == 'ICA':
        model = FastICA(**kwargs, max_iter=5000, tol=1e-4,
                        random_state=config.seed, whiten='unit-variance')
    elif method == 'Isomap':
        model = Isomap(**kwargs, n_neighbors=config.n_neighbors, n_jobs=1)
    elif method == 'UMAP':
        import umap
        model = umap.UMAP(**kwargs, n_neighbors=config.n_neighbors,
                          random_state=config.seed, n_jobs=1)
    else:
        model = TSNE(**kwargs, random_state=config.seed, init='pca', n_jobs=1)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit_start = time.perf_counter()
        z = model.fit_transform(x)
        fit_transform_seconds = time.perf_counter() - fit_start
    z = np.asarray(z, dtype=np.float32)
    if z.shape != (x.shape[0], config.n_components) or not np.isfinite(z).all():
        raise ValueError('reducer返回的坐标形状或数值无效')
    metadata = _fit_metadata(model, caught)
    metadata.update(fit_transform_seconds=fit_transform_seconds,
                    timing_boundary='model.fit_transform(X) only')
    if method == 'PCA':
        metadata['model_variance'] = {'explained_variance_ratio': model.explained_variance_ratio_.tolist()}
    elif method == 'FA':
        metadata['model_variance'] = {
            'fa_signal_variance': float(np.square(model.components_).sum()),
            'fa_noise_variance': float(model.noise_variance_.sum())}
    return z, metadata


def cluster_parameters(method, dr_method, config):
    """策略按预先指定DR名称决定，禁止按实际维数静默切换GMM模型。"""
    if dr_method not in DR_ORDER:
        raise ValueError(f'未知降维方法：{dr_method}')
    if method == 'KMeans':
        return {'n_clusters': config.k, 'n_init': 10, 'random_state': config.seed}
    if method == 'GMM':
        covariance = config.gmm_covariance
        if covariance == 'full_2d_diag_none':
            covariance = 'diag' if dr_method == 'none' else 'full'
        return {'n_components': config.k, 'reg_covar': 1e-4,
                'random_state': config.seed, 'covariance_type': covariance}
    if method == 'HDBSCAN':
        return {'min_cluster_size': config.min_cluster_size,
                'min_samples': config.min_samples, 'core_dist_n_jobs': 1}
    raise ValueError(f'未知聚类方法：{method}')


def summarize_labels(labels):
    labels = np.asarray(labels)
    if labels.ndim != 1 or not labels.size or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError('聚类标签必须是非空一维整数数组')
    if (labels < -1).any():
        raise ValueError('只允许-1表示noise，其余标签必须非负')
    n_noise = int((labels == -1).sum())
    n_clusters = len(np.unique(labels[labels != -1]))
    return {'n_cells': len(labels), 'n_clusters': n_clusters, 'n_noise': n_noise,
            'noise_frac': float(n_noise / len(labels)),
            'outcome': 'all_noise' if n_noise == len(labels) else 'clustered'}


def cluster_matrix(z, method, dr_method, config):
    """不接收gold；K23来自用户预先指定，HDBSCAN保持k-free。"""
    from sklearn.cluster import KMeans
    from sklearn.mixture import GaussianMixture

    z = _validate_matrix(z)
    params = cluster_parameters(method, dr_method, config)
    if method == 'KMeans':
        model = KMeans(**params)
    elif method == 'GMM':
        model = GaussianMixture(**params)
    else:
        import hdbscan
        model = hdbscan.HDBSCAN(**params)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        labels = np.asarray(model.fit_predict(z))
    if labels.shape != (z.shape[0],):
        raise ValueError('聚类标签长度与输入细胞数不匹配')
    meta = _fit_metadata(model, caught)
    meta.update(summarize_labels(labels))
    return labels.astype(np.int64, copy=False), meta
