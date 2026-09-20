"""Saved-prediction cross-tool agreement and frozen NL022 expression evidence.

No fitting or annotation changes. Every data operation requires SLURM.
"""
from pathlib import Path
import json
import os
import sys
import traceback
from itertools import combinations

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE = ROOT/'handoff/reviewer_completion_20260920/evidence'
DEST = ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/evidence'
sys.path.insert(0, str(ROOT/'handoff/paper_claim_validation_20260917'))
from common import OUT, L1, require_slurm, checked, sha, write_json, utc, complete
from prediction_helpers import read_native, map_labels, UNKNOWN

METHODS = ['DG-scRNA', 'scType', 'scCATCH', 'SCINA', 'SingleR', 'scDeepSort']
SOURCES = {}


def source(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def status(state, completed, remaining, **extra):
    write_json(DEST/'status.json', dict(stage='D_F_GBM', status=state, updated_at=utc(),
        jobs=[os.environ.get('SLURM_JOB_ID', '')], completed=completed, remaining=remaining,
        evidence=[str(DEST.relative_to(ROOT))], **extra))


def savefig(fig, stem):
    import matplotlib.pyplot as plt
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(DEST/f'{stem}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)


def summarize_patients(frame, keys, values, stem):
    """Sample means within patient, then equally weighted patients; conditional CI."""
    import numpy as np
    import pandas as pd
    patient = frame.groupby(keys+['patient'], dropna=False)[values].mean().reset_index()
    patient.to_csv(DEST/f'{stem}_patient.csv', index=False)
    rng = np.random.default_rng(20260920)
    rows = []
    for config, part in patient.groupby(keys, dropna=False, sort=True):
        if not isinstance(config, tuple): config = (config,)
        base = dict(zip(keys, config))
        for value in values:
            x = part[value].dropna().to_numpy()
            draws = x[rng.integers(0, len(x), size=(2000, len(x)))].mean(1) if len(x) else []
            rows.append(dict(**base, metric=value, n_patients=len(x),
                mean=float(x.mean()) if len(x) else float('nan'),
                CI95_low=float(np.quantile(draws, .025)) if len(x) else float('nan'),
                CI95_high=float(np.quantile(draws, .975)) if len(x) else float('nan')))
    summary = pd.DataFrame(rows)
    summary.to_csv(DEST/f'{stem}_summary.csv', index=False)
    return patient, summary


def agreement():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, fowlkes_mallows_score
    from sklearn.metrics import precision_recall_fscore_support
    assert checked(OUT/'comparison_summary')
    cohort = pd.read_csv(source(OUT/'protocol/cohort.csv'))
    cv = pd.read_csv(source(OUT/'comparison_summary/patient_heldout_results.csv'), dtype={'cutoff': str})
    source(OUT/'comparison_summary/manifest.json')
    source(OUT/'markers/panel_L1_mapping.csv')
    records, pairrows, grouprows, classrows, provenance, parity = [], [], [], [], [], []
    n_unique = set()
    for sample in cohort.itertuples():
        truthfile = source(OUT/'evaluation_inputs'/sample.sample/'truth.csv.gz')
        truth = pd.read_csv(truthfile, dtype=str, keep_default_na=False)
        assert truth.cell_id.is_unique and len(truth) == sample.n_cells
        assert truth.L1.isin(L1).all()
        n_unique.add(sample.sample)
        for c in ['all121', 'primary97']:
            if c == 'primary97' and not sample.primary: continue
            chosen = cv[(cv.cohort == c) & (cv.patient == sample.patient)]
            assert set(chosen.method) == set(METHODS) and len(chosen) == len(METHODS)
            mapped, native_called = [], []
            cells = truth[['cell_id', 'L1']].copy()
            cells.insert(0, 'sample', sample.sample)
            cells.insert(1, 'patient', sample.patient)
            cells.insert(2, 'cohort', c)
            for method in METHODS:
                cfg = chosen[chosen.method == method].iloc[0]
                pred, path = read_native(method, sample.sample, cfg.budget, cfg.route, cfg.library, cfg.cutoff)
                source(path)
                assert np.array_equal(pred.cell_id, truth.cell_id), (sample.sample, method)
                p = map_labels(method, cfg.library, pred.prediction)
                native = ~pred.prediction.isin(UNKNOWN).to_numpy()
                supported = np.isin(p, L1)
                mapped.append(p); native_called.append(native)
                cells[method] = p
                _, _, f1, support = precision_recall_fscore_support(truth.L1, p, labels=L1, zero_division=0)
                metrics = dict(accuracy=float((p == truth.L1).mean()),
                    macroF1_present=float(f1[support > 0].mean()), coverage=float(native.mean()),
                    unknown_rate=float((~native).mean()), mapped_coverage=float(supported.mean()),
                    off_vocabulary_rate=float((native & ~supported).mean()))
                evaluation = OUT/'GBM'/sample.sample/'hvg2000/evaluation' if method == 'DG-scRNA' else OUT/'comparators'/method/sample.sample/'evaluation'
                saved = pd.read_csv(evaluation/'metrics.csv', dtype={'cutoff': str})
                mask = np.ones(len(saved), dtype=bool)
                for key in ['budget', 'route', 'library', 'cutoff']: mask &= saved[key].eq(cfg[key]).to_numpy()
                if method == 'DG-scRNA': mask &= saved.stage.eq('terminal090').to_numpy()
                selected = saved[mask]; assert len(selected) == 1
                for key, value in metrics.items():
                    np.testing.assert_allclose(value, selected.iloc[0][key], rtol=0, atol=1e-12,
                        err_msg=str((sample.sample, method, key)))
                parity.append(dict(cohort=c, sample=sample.sample, method=method, n_cells=len(truth),
                    six_metrics_match_saved=True))
                records.append(dict(cohort=c, sample=sample.sample, patient=sample.patient, method=method,
                    n_cells=len(truth), **metrics,
                    ARI_author_L1=adjusted_rand_score(truth.L1, p),
                    NMI_author_L1=normalized_mutual_info_score(truth.L1, p),
                    FMI_author_L1=fowlkes_mallows_score(truth.L1, p)))
                provenance.append(dict(cohort=c, sample=sample.sample, patient=sample.patient, method=method,
                    budget=cfg.budget, route=cfg.route, library=cfg.library, cutoff=cfg.cutoff,
                    fold=int(cfg.fold), endpoint='terminal090' if method == 'DG-scRNA' else cfg.cutoff,
                    source=str(path.relative_to(ROOT)), sha256=sha(path)))
            matrix = np.asarray(mapped).T
            called = np.asarray(native_called).T
            supported = np.isin(matrix, L1)
            y = truth.L1.to_numpy()
            for i, j in combinations(range(len(METHODS)), 2):
                a, b = matrix[:, i], matrix[:, j]
                for subset, mask in [('all_cells', np.ones(len(y), dtype=bool)),
                    ('both_native_called', called[:, i] & called[:, j]),
                    ('both_supported_L1', supported[:, i] & supported[:, j])]:
                    n = int(mask.sum()); same = a == b
                    good = same & supported[:, i] & supported[:, j]
                    pairrows.append(dict(cohort=c, sample=sample.sample, patient=sample.patient,
                        method_a=METHODS[i], method_b=METHODS[j], subset=subset, n_all_cells=len(y), n_subset=n,
                        subset_coverage=n/len(y), n_same_label=int((same & mask).sum()),
                        n_same_supported_label=int((good & mask).sum()),
                        exact_agreement=float(same[mask].mean()) if n else np.nan,
                        supported_agreement=float(good[mask].mean()) if n else np.nan,
                        joint_accuracy=float(((a == y) & (b == y))[mask].mean()) if n else np.nan,
                        ARI=adjusted_rand_score(a[mask], b[mask]) if n >= 2 else np.nan,
                        NMI=normalized_mutual_info_score(a[mask], b[mask]) if n >= 2 else np.nan,
                        FMI=fowlkes_mallows_score(a[mask], b[mask]) if n >= 2 else np.nan))
            for group, cols in [('all_six', list(range(6))), ('four_marker_methods', list(range(4)))]:
                pm = matrix[:, cols]
                all_supported = supported[:, cols].all(1)
                exact = (pm == pm[:, :1]).all(1)
                unanimous = exact & all_supported
                counts = np.stack([(pm == label).sum(1) for label in L1], axis=1)
                max_count = counts.max(1)
                # Strict majority over all group tools, not only those making a supported call.
                majority = np.asarray(L1)[counts.argmax(1)]
                majority[max_count <= len(cols)/2] = 'No_majority'
                cells[group+'_unanimous_supported'] = unanimous
                cells[group+'_majority'] = majority
                for label in ['ALL']+L1:
                    labelmask = np.ones(len(y), dtype=bool) if label == 'ALL' else y == label
                    for category, mask in [('all_cells', labelmask),
                        ('unanimous_supported', labelmask & unanimous),
                        ('nonunanimous_or_unsupported', labelmask & ~unanimous),
                        ('all_supported_but_disagree', labelmask & all_supported & ~exact),
                        ('any_unsupported', labelmask & ~all_supported)]:
                        n = int(mask.sum())
                        values = dict(cohort=c, sample=sample.sample, patient=sample.patient, group=group,
                            author_L1=label, category=category, n_all_cells=len(y), n_class=int(labelmask.sum()),
                            n_subset=n, fraction_all_cells=n/len(y),
                            fraction_within_class=n/int(labelmask.sum()) if labelmask.any() else np.nan,
                            mean_tool_accuracy=float((pm == y[:, None])[mask].mean()) if n else np.nan,
                            strict_majority_accuracy=float((majority == y)[mask].mean()) if n else np.nan,
                            strict_majority_coverage=float((majority != 'No_majority')[mask].mean()) if n else np.nan,
                            unanimous_accuracy=float((pm[:, 0] == y)[mask].mean()) if n and category == 'unanimous_supported' else np.nan,
                            raw_exact_fraction=float(exact[mask].mean()) if n else np.nan)
                        (grouprows if label == 'ALL' else classrows).append(values)
            # Frozen labels retained per cell; no consensus labels replace method outputs.
            target = DEST/'predictions'/c
            target.mkdir(parents=True, exist_ok=True)
            cells.to_csv(target/f'{sample.sample}.csv.gz', index=False, compression='gzip')
        if len(n_unique) % 10 == 0:
            print('AGREEMENT_SAMPLE', len(n_unique), sample.sample, flush=True)
            status('running', [f'cross-tool evaluation: {len(n_unique)}/121 samples'],
                ['remaining agreement samples', 'NL022 expression figures'])
    assert len(n_unique) == 121
    pairs = pd.DataFrame(pairrows); metrics = pd.DataFrame(records); groups = pd.DataFrame(grouprows)
    pairs.to_csv(DEST/'pairwise_agreement_sample.csv.gz', index=False)
    metrics.to_csv(DEST/'method_author_concordance_sample.csv', index=False)
    groups.to_csv(DEST/'consensus_strata_sample.csv', index=False)
    pd.DataFrame(classrows).to_csv(DEST/'consensus_author_L1_strata_sample.csv.gz', index=False)
    pd.DataFrame(provenance).to_csv(DEST/'selected_prediction_sources.csv', index=False)
    pd.DataFrame(parity).to_csv(DEST/'saved_metric_parity.csv', index=False)
    _, pairsummary = summarize_patients(pairs, ['cohort','method_a','method_b','subset'],
        ['subset_coverage','exact_agreement','supported_agreement','joint_accuracy','ARI','NMI','FMI'], 'pairwise_agreement')
    summarize_patients(metrics, ['cohort','method'],
        ['accuracy','macroF1_present','coverage','unknown_rate','mapped_coverage','off_vocabulary_rate',
         'ARI_author_L1','NMI_author_L1','FMI_author_L1'], 'method_author_concordance')
    _, groupsummary = summarize_patients(groups, ['cohort','group','category'],
        ['fraction_all_cells','mean_tool_accuracy','strict_majority_accuracy','strict_majority_coverage',
         'unanimous_accuracy','raw_exact_fraction'], 'consensus_strata')
    summarize_patients(pd.DataFrame(classrows), ['cohort','group','category','author_L1'],
        ['fraction_within_class','mean_tool_accuracy','strict_majority_accuracy'], 'consensus_author_L1_strata')
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), layout='constrained')
    for ax, metric, title in zip(axes, ['supported_agreement','FMI','ARI'],
        ['Same supported L1 / all cells', 'FMI of annotation partitions', 'ARI of annotation partitions']):
        matrix = np.full((6, 6), np.nan)
        sub = pairsummary[(pairsummary.cohort == 'primary97') & (pairsummary.subset == 'all_cells') & (pairsummary.metric == metric)]
        for row in sub.itertuples():
            i, j = METHODS.index(row.method_a), METHODS.index(row.method_b)
            matrix[i, j] = matrix[j, i] = row.mean
        image = ax.imshow(np.ma.masked_invalid(matrix), cmap='viridis', vmin=0, vmax=1)
        ax.set_xticks(range(6), METHODS, rotation=40, ha='right', fontsize=8)
        ax.set_yticks(range(6), METHODS, fontsize=8)
        for i in range(6):
            for j in range(6):
                if i != j: ax.text(j, i, f'{matrix[i,j]:.2f}', ha='center', va='center', fontsize=8,
                    color='white' if matrix[i,j] < .5 else 'black')
        ax.set_title(title, fontsize=10)
    fig.colorbar(image, ax=axes, shrink=.6)
    fig.suptitle('Primary97 / 55 patients: true cross-tool comparisons\nPatient means of sample metrics; shared abstention is not supported agreement', fontsize=12)
    savefig(fig, 'GBM_cross_tool_agreement')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    colors = ['#0072B2', '#D55E00']
    labels = ['Unanimous supported', 'Disagreement / unsupported']
    for idx, group in enumerate(['four_marker_methods', 'all_six']):
        sub = groupsummary[(groupsummary.cohort == 'primary97') & (groupsummary.group == group)]
        for axis, metric in zip(axes, ['fraction_all_cells', 'mean_tool_accuracy']):
            for j, cat in enumerate(['unanimous_supported', 'nonunanimous_or_unsupported']):
                row = sub[(sub.category == cat) & (sub.metric == metric)].iloc[0]
                x = idx+(j-.5)*.34
                axis.bar(x, row['mean'], width=.32, color=colors[j], label=labels[j] if idx == 0 else None)
                axis.errorbar(x, row['mean'], yerr=[[row['mean']-row.CI95_low],[row.CI95_high-row['mean']]],
                    color='black', capsize=3, lw=1)
            axis.set_xticks([0,1], ['Four marker tools','All six tools'])
            axis.set_ylim(0,1)
    axes[0].set_ylabel('Fraction of all cells'); axes[1].set_ylabel('Mean tool accuracy versus author L1')
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle('Agreement is not a truth standard\nPatient means and conditional 95% bootstrap intervals', fontsize=12)
    savefig(fig, 'GBM_consensus_vs_author_L1')
    return dict(n_samples=121, n_patients=cohort.patient.nunique(), n_cohort_sample_units=218,
        n_pairwise_rows=len(pairs), n_metric_parity_rows=len(parity), methods=METHODS)


def expression():
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sample = json.loads(source(OUT/'protocol/input_audit.json').read_text())['display_sample']
    assert sample == 'NL022'
    cohort = pd.read_csv(OUT/'protocol/cohort.csv'); patient = cohort[cohort['sample'] == sample].iloc[0].patient
    cv = pd.read_csv(source(OUT/'DG_fixed_partition_selection/patient_heldout_results.csv'), dtype={'cutoff': str})
    cfg = cv[(cv.cohort == 'primary97') & (cv.patient == patient) & (cv.method == 'DG-scRNA')].iloc[0]
    prep = OUT/'GBM'/sample/'hvg2000'
    src = OUT/'inputs'/sample
    im = json.loads(source(src/'input_manifest.json').read_text())
    assert checked(src, 'input_manifest.json', 'INPUT_COMPLETE')
    for file in ['x.bin','i.bin','p.bin','genes.csv','cells_fit.csv']:
        assert sha(source(src/file)) == im['fitting_files'][file]
    X = sp.csr_matrix((np.fromfile(src/'x.bin', dtype='<f8'), np.fromfile(src/'i.bin', dtype='<i4'),
        np.fromfile(src/'p.bin', dtype='<i4')), shape=(im['n_cells'], im['n_genes']))
    X.eliminate_zeros(); assert (X.data > 0).all()
    truth = pd.read_csv(source(OUT/'evaluation_inputs'/sample/'truth.csv.gz'), dtype=str, keep_default_na=False)
    fit = pd.read_csv(src/'cells_fit.csv', dtype=str)
    assert list(truth.cell_id) == list(fit.cell_id)
    genes = pd.read_csv(source(prep/'Seurat_gene_names.csv'), dtype=str)
    assert list(genes.source) == list(pd.read_csv(src/'genes.csv').gene)
    assert genes.Seurat.is_unique
    pos = {g:i for i,g in enumerate(genes.Seurat)}
    totals = np.asarray(X.sum(1)).ravel()
    X.data = np.log1p(X.data*np.repeat(10000/totals, np.diff(X.indptr)))
    features = source(prep/'DL_features.txt').read_text().splitlines()
    dl = np.memmap(source(prep/'DL.float32.bin'), dtype='<f4', mode='r', shape=(len(truth),len(features)))
    # Full NL022 HVG parity, not a sampled normalization check.
    diff = float(np.max(np.abs(X[:,[pos[g] for g in features]].toarray().astype(np.float32)-dl)))
    assert diff < 1e-6, diff
    coords = pd.read_csv(source(prep/'UMAP2.csv'), index_col=0)
    assert list(coords.index) == list(truth.cell_id)
    libraries = json.loads(source(OUT/'markers/libraries.json').read_text())
    mapping = pd.read_csv(OUT/'markers/panel_L1_mapping.csv', dtype=str)
    conditions = [('fixed_glioma','CM2_glioma_other','mean'), ('training_patient_selected',cfg.library,cfg.cutoff)]
    labelsets = {'author_L1':truth.L1.to_numpy()}
    labels = truth[['cell_id','L1']].copy()
    provenance, markerrows = [], []
    detected = X.getnnz(0)/len(truth)
    mean_all = np.asarray(X.mean(0)).ravel()
    plot_genes = []
    for condition, library, cutoff in conditions:
        pred, path = read_native('DG-scRNA', sample, 'hvg2000', 'UMAP2_HDBSCAN_R', library, cutoff)
        source(path)
        assert list(pred.cell_id) == list(truth.cell_id)
        source(path.parent/'terminal_manifest.json')
        terminal = json.loads((path.parent/'terminal_manifest.json').read_text())
        provenance.append(dict(condition=condition, sample=sample, patient=patient, library=library, cutoff=cutoff,
            budget='hvg2000', route='UMAP2_HDBSCAN_R', endpoint='final090',
            terminal_manifest=terminal, prediction_source=str(path.relative_to(ROOT))))
        for endpoint, values in [('initial',pred.initial), ('terminal090',pred.prediction)]:
            key=condition+'_'+endpoint
            labelsets[key] = map_labels('DG-scRNA', library, values)
            labels[key+'_native'] = values
            labels[key+'_L1'] = labelsets[key]
        localmap = mapping[mapping.library == library].set_index('panel').L1.to_dict()
        for panel, markers in libraries[library].items():
            for gene in markers:
                # Original Seurat names are preserved; no marker renaming here.
                ix = pos.get(gene)
                markerrows.append(dict(condition=condition, library=library, panel=panel, mapped_L1=localmap[panel],
                    gene=gene, in_scoring_RNA=ix is not None, in_DL_HVG=gene in features,
                    detection_fraction_all_cells=float(detected[ix]) if ix is not None else np.nan,
                    mean_log1p_all_cells=float(mean_all[ix]) if ix is not None else np.nan))
    markerdf = pd.DataFrame(markerrows)
    markerdf.to_csv(DEST/'NL022_selected_library_marker_inventory.csv', index=False)
    labels.to_csv(DEST/'NL022_frozen_endpoint_labels.csv.gz', index=False)
    marker_detail = pd.read_csv(source(OUT/'marker_evidence_summary/gene_panel_source_records.csv.gz'), dtype=str, keep_default_na=False)
    marker_detail[marker_detail.library.isin([a[1] for a in conditions])].to_csv(
        DEST/'NL022_selected_library_gene_source_records.csv.gz', index=False)
    # Three display genes per supported L1 class, chosen on overall detection only.
    # This ranks presentation genes, not marker libraries or annotation parameters.
    display = markerdf[(markerdf.condition == 'training_patient_selected') & markerdf.in_scoring_RNA
        & markerdf.mapped_L1.isin(L1)].drop_duplicates(['mapped_L1','gene'])
    display = display.sort_values(['mapped_L1','detection_fraction_all_cells','gene'], ascending=[True,False,True])
    display = display.groupby('mapped_L1', sort=False).head(3)
    display.to_csv(DEST/'NL022_display_gene_selection.csv', index=False)
    plot_genes = display.gene.drop_duplicates().tolist()
    assert len(plot_genes) > 0
    all_markers = markerdf[markerdf.in_scoring_RNA].gene.drop_duplicates().tolist()
    expression_rows = []
    for endpoint, values in labelsets.items():
        for label in sorted(set(values)):
            mask = values == label; sub = X[mask][:,[pos[g] for g in all_markers]]
            means = np.asarray(sub.mean(0)).ravel(); fractions = sub.getnnz(0)/mask.sum()
            expression_rows.extend(dict(endpoint=endpoint, label=label, n_cells=int(mask.sum()), gene=gene,
                mean_log1p=float(avg), detection_fraction=float(frac))
                for gene, avg, frac in zip(all_markers,means,fractions))
    expression_summary = pd.DataFrame(expression_rows)
    expression_summary.to_csv(DEST/'NL022_all_selected_marker_expression_by_endpoint.csv.gz', index=False)
    pd.DataFrame(X[:,[pos[g] for g in plot_genes]].toarray(), index=truth.cell_id, columns=plot_genes).rename_axis('cell_id').to_csv(
        DEST/'NL022_display_gene_cell_expression.csv.gz')
    dots = [('author_L1','Author L1'),('training_patient_selected_initial','Selected initial marker labels'),
        ('training_patient_selected_terminal090','Selected terminal DL / 0.90')]
    fig, axes = plt.subplots(3,1,figsize=(max(12,.38*len(plot_genes)),10),layout='constrained')
    max_mean = expression_summary[expression_summary.gene.isin(plot_genes)].mean_log1p.max()
    for ax, (endpoint,title) in zip(axes,dots):
        rows = [x for x in L1+['Unknown','UNMAPPABLE','AMBIGUOUS_NEURON','NO_L1_COUNTERPART'] if x in set(labelsets[endpoint])]
        rows += sorted(set(labelsets[endpoint])-set(rows))
        table = expression_summary[expression_summary.endpoint == endpoint]
        for i,label in enumerate(rows):
            frame = table[table.label == label].set_index('gene').reindex(plot_genes)
            image = ax.scatter(range(len(plot_genes)),np.full(len(plot_genes),i),
                s=frame.detection_fraction*110,c=frame.mean_log1p,cmap='viridis',vmin=0,vmax=max_mean,linewidths=.25,edgecolors='#888')
        ax.set_xticks(range(len(plot_genes)),plot_genes,rotation=60,ha='right',fontsize=8)
        ax.set_yticks(range(len(rows)),[f'{lab} (n={(labelsets[endpoint] == lab).sum()})' for lab in rows],fontsize=8)
        ax.set_title(title,fontsize=10); ax.invert_yaxis(); ax.grid(alpha=.15)
        ax.set_xlim(-.6,len(plot_genes)-.4)
    fig.colorbar(image,ax=axes,shrink=.6,label='Mean log1p(count / cell total × 10,000)')
    fig.legend(handles=[Line2D([],[],marker='o',ls='',markerfacecolor='#666',markeredgecolor='#888',
        markersize=np.sqrt(frac*110),label=f'{frac:.0%} detected') for frac in [.1,.5,1]],
        loc='outside lower center',ncol=3,frameon=False,fontsize=8)
    fig.suptitle(f'NL022: markers from {cfg.library}\nUp to 3 genes per supported class, ranked by label-free overall detection; full marker table retained',fontsize=11)
    savefig(fig,'NL022_selected_marker_dotplot')
    # Violin selection: first detection-ranked gene per supported class, <=11 panels.
    violin_genes = display.groupby('mapped_L1',sort=False).head(1).gene.drop_duplicates().tolist()
    endpoint='training_patient_selected_terminal090'; values=labelsets[endpoint]
    rows=[lab for lab in L1+['Unknown','UNMAPPABLE','AMBIGUOUS_NEURON','NO_L1_COUNTERPART'] if lab in set(values)]
    ncols=3; nrows=(len(violin_genes)+ncols-1)//ncols
    fig, axes=plt.subplots(nrows,ncols,figsize=(14,3.4*nrows),squeeze=False,layout='constrained')
    for ax,gene in zip(axes.flat,violin_genes):
        vector=X[:,pos[gene]].toarray().ravel()
        data=[vector[values==lab] for lab in rows]
        for i,a in enumerate(data):
            if len(a)>1 and np.ptp(a)>0:
                parts=ax.violinplot([a],positions=[i],showmedians=True,showextrema=False)
                for body in parts['bodies']:body.set_facecolor('#0072B2');body.set_alpha(.55)
            else:ax.plot([i-.15,i+.15],[a[0],a[0]],color='#0072B2')
        ax.set_xticks(range(len(rows)),rows,rotation=45,ha='right',fontsize=7)
        ax.set_title(gene);ax.set_ylabel('log1p normalized RNA',fontsize=8)
    for ax in list(axes.flat)[len(violin_genes):]:ax.axis('off')
    fig.suptitle('NL022: terminal DL labels, including abstentions\nOne detected marker per supported class; densities are descriptive, not independent validation',fontsize=12)
    savefig(fig,'NL022_selected_marker_violin')
    # Same frozen UMAP display coordinates for every endpoint and both marker choices.
    palette=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#332288','#88CCEE','#44AA99','#AA4499','#999933']
    colors=dict(zip(L1,palette));colors.update(Unknown='#bdbdbd',UNMAPPABLE='#333333',AMBIGUOUS_NEURON='#777777',NO_L1_COUNTERPART='#666666')
    panels=[('author_L1','Author L1'),('fixed_glioma_initial','Fixed glioma / initial'),
        ('fixed_glioma_terminal090','Fixed glioma / terminal DL'),('training_patient_selected_initial','Selected library / initial'),
        ('training_patient_selected_terminal090','Selected library / terminal DL')]
    fig, axes=plt.subplots(2,3,figsize=(15,10),layout='constrained')
    for ax,(endpoint,title) in zip(axes.flat,panels):
        values=labelsets[endpoint]
        for label in sorted(set(values)):
            mask=values==label
            ax.scatter(coords.iloc[mask,0],coords.iloc[mask,1],s=3,c=[colors.get(label,'#555555')],linewidths=0,rasterized=True)
        ax.set_title(title,fontsize=11);ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel('Fixed display UMAP1');ax.set_ylabel('Fixed display UMAP2')
    unknown=labelsets['training_patient_selected_terminal090']=='Unknown'
    axes.flat[-1].scatter(coords.iloc[:,0],coords.iloc[:,1],s=3,c=np.where(unknown,'#D55E00','#d9d9d9'),linewidths=0,rasterized=True)
    axes.flat[-1].set_title(f'Selected terminal Unknown: {unknown.sum()} / {len(unknown)}',fontsize=11)
    axes.flat[-1].set_xticks([]);axes.flat[-1].set_yticks([])
    present=set(np.concatenate(list(labelsets.values())))
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=colors.get(lab,'#555555'),label=lab,markersize=5)
        for lab in L1+['Unknown','UNMAPPABLE','AMBIGUOUS_NEURON','NO_L1_COUNTERPART'] if lab in present],
        loc='outside lower center',ncol=6,frameon=False,fontsize=8)
    fig.suptitle('NL022 / HVG2000 / native-R UMAP→HDBSCAN / terminal 0.90\nFixed glioma mean cutoff versus training-patient-selected '+cfg.library+' / '+cfg.cutoff,fontsize=12)
    savefig(fig,'NL022_markers_initial_terminal_same_coordinates')
    write_json(DEST/'NL022_expression_provenance.json',dict(sample=sample,n_cells=len(truth),
        normalization='log1p(count / retained-gene cell total * 10000)',
        normalization_all_HVG_max_abs_R_difference=diff,normalization_checked_cells=len(truth),
        normalization_checked_genes=len(features),conditions=provenance,
        n_plot_genes=len(plot_genes),n_all_marker_genes_present=len(all_markers),
        no_new_fitting=True,display_sample_chosen_by_median_size_before_prediction_inspection=True,
        gene_display_rule='Top3 available genes per supported L1 by overall detection, gene-name tie break; no truth labels used'))
    return dict(display_sample=sample,n_cells=len(truth),normalization_max_abs_difference=diff,
        n_marker_gene_panel_pairs=len(markerdf),n_available_marker_genes=len(all_markers),n_display_genes=len(plot_genes))


def run():
    require_slurm()
    DEST.mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    status('running',[],['cross-tool agreement','NL022 marker expression'])
    a=agreement()
    status('running',['121-sample cross-tool agreement and author-L1 parity'],['NL022 marker expression'])
    e=expression()
    source(__file__)
    note='''# GBM selected-marker expression and true cross-tool agreement

This supplement evaluates saved predictions; it neither refits a model nor changes
any original labels. DG-scRNA is always terminal DL/refinement at confidence 0.90.
Six tools are DG-scRNA, scType, scCATCH, SCINA, SingleR and scDeepSort. The first
four are marker methods; SingleR uses labelled training-patient references and
scDeepSort uses its published human brain atlas. Their information conditions
remain distinct. Existing five-fold patient-label-heldout library/threshold
choices are reused as retrospective analyses, not newly blinded evaluations.

## Denominators and interpretation

All121 comprises 121 samples / 59 patients and primary97 comprises 97 samples /
55 patients. No cell is removed from the all-cell denominator. Pairwise tables
also report both-native-called and both-supported-L1 subsets and their coverage.
Native-call coverage differs from supported-L1 coverage: off-vocabulary calls
remain separate from abstention and count as errors against author L1. Exact
agreement includes matching Unknown/unsupported strings and is explicitly
distinct from supported agreement, which requires the same genuine L1 label.

All-tool unanimity requires every method to give the same supported L1 class.
Shared Unknown or unmappable labels never create a biological consensus. Strict
majority needs more than half of all tools in that group; abstainers remain in
the voting denominator. Consensus predictions are diagnostics and do not replace
the methods' saved labels. Author L1 is the existing annotation reference, not
independent biological gold truth. Consensus does not establish correctness.

ARI/NMI/FMI compare partitions induced by final mapped annotations. They are
not accuracy scores and are not the upstream cluster partitions; every mapped
status is a distinct category in all-cell partition metrics. FMI can be high
under class imbalance and label permutations. Conditional supported-subset
metrics are supplied alongside coverage; they cannot replace all-cell results.
Undefined empty strata are NA, not zero. Per-class consensus files expose all
11 author categories including unsupported or empty strata.

Samples are averaged within patient before cohort averaging. Confidence
intervals use 2000 patient bootstrap draws conditional on the fixed selected
predictions; they do not refit the cross-validation selection and are not an
independent significance test of methods sharing training folds. The earlier
paired method tests are reused, not duplicated or reinterpreted here.

## NL022 expression evidence

The previously frozen median-size display sample NL022 is used, with exactly
the original HVG2000 native-R UMAP coordinates. Dotplots and violins use genuine
log-normalized RNA expression reconstructed from the frozen eligible-gene counts.
Every cell and all 2000 HVGs are checked numerically against R's exported matrix.
All selected-library marker genes, missing-gene states and original gene-specific
source/assay records are retained. Assay metadata are not recoded as RNA-only.

The fixed CM2_glioma_other / mean context and the training-patient-selected
library / cutoff are shown on the same partition. Figures distinguish initial
marker calls, terminal DL calls and Unknown. Dotplots show the fraction with
nonzero expression (size) and mean log1p expression over all cells in that row
(color). Up to three display genes per supported class are ranked by overall
detection independent of labels, with lexicographic ties. Violin panels use one
such gene per supported class. Full gene summaries remain available, preventing
the limited presentation panel from concealing the remainder of the marker set.

Expression separation is partially induced by using these same genes to annotate
cells and is descriptive consistency evidence, not independent validation of
biological cell types, disease mechanisms, or pipeline optimality. No PTC result,
SignacX historical zero-T row, old notebook or earlier result is modified.
'''
    (DEST/'README.md').write_text(note)
    write_json(DEST/'manifest.json',dict(status='completed',agreement=a,expression=e,
        files={str(p.relative_to(DEST)):sha(p) for p in DEST.rglob('*') if p.is_file() and p.name not in ['manifest.json','status.json','COMPLETE']},
        sources=SOURCES,job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(DEST)
    status('completed',['121 samples / 6 tools / patient-aggregated agreement, ARI/NMI/FMI',
        'Unknown and coverage denominators; true consensus versus author L1',
        'NL022 selected-library RNA dotplot, violin and same-coordinate endpoints; full-R normalization parity'],[],
        summary=dict(agreement=a,expression=e))
    print('EVIDENCE_COMPLETE',json.dumps(dict(agreement=a,expression=e)),flush=True)


if __name__=='__main__':
    try: run()
    except Exception as exc:
        status('failed',[],['repair evidenced failure and resume'],error=str(exc),traceback=traceback.format_exc())
        raise
