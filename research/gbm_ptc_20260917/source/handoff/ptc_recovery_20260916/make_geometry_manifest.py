"""Create predetermined configurations; no data/labels are read."""
from ptc_common import BASE, GROUPS, write_json

tasks = []
for group, samples in GROUPS.items():
    for sample in samples:
        specs = []
        for feature in ['all', 'hvg500', 'hvg1000', 'hvg2000', 'hvg3000', 'hvg5000']:
            for dr in ['PCA', 'FA', 'ICA', 'Isomap', 'UMAP', 'TSNE', 'none']:
                specs.append((feature, 'direct', dr, None if dr=='none' else 2, 42, 'factorial'))
        for feature in ['all', 'hvg2000']:
            specs.extend([(feature, 'pca30', 'UMAP', 2, 42, 'PCA_control'),
                          (feature, 'direct', 'PCA', 30, 42, 'PCA_control')])
            for seed in [7, 17, 29, 101]:
                for space in ['direct', 'pca30']:
                    specs.append((feature, space, 'UMAP', 2, seed, 'seed_sensitivity'))
        for feature, space, dr, dim, seed, family in specs:
            gid = f'{feature}__{space}__{dr}{dim or 0}__s{seed}'
            full_roster = seed==42 and feature in ['all', 'hvg2000'] and (dr=='UMAP' or dim==30)
            tasks.append(dict(task_id=len(tasks), group=group, sample=sample, geometry_id=gid,
                feature=feature, input_space=space, dr=dr, dim=dim, seed=seed,
                neighbors=15, min_dist=0.1, family=family, full_roster=full_roster,
                clusterers=['KMeans', 'GMM', 'HDBSCAN'] if seed==42 else ['HDBSCAN']))
write_json(BASE/'protocol/single_sample_geometries.json', tasks)
print(f'{len(tasks)} frozen geometry tasks; no reference annotations used')
