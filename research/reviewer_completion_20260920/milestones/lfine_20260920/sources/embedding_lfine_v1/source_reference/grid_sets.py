"""Cross clustering granularity with scoring formula, on the same cells.

The two candidate levers for raising macro-F1 are (a) finer clustering, which raises the
attainable ceiling, and (b) replacing the panel-length divisor, which recovers score the
current formula throws away. Each has been measured alone; this measures the grid, plus
the ceiling at every granularity so gains are always reported against what is reachable.
"""
import sys, os, json
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.io as sio
import hdbscan

B = "/fs/scratch/PCON0080/yimin/dgscrna"
OUT = f"{B}/results/grid_sets"; os.makedirs(OUT, exist_ok=True)
SEED = 42
# All 8 sets of the markers_v3 roster, so the cohort figure has the same rows as the
# single-sample figure (handoff/FROZEN_single_sample_TKU3186). CARE_TME and BrainAtlas112
# built the golden label and must be marked with a dagger wherever they are plotted.
SETS = ["CM2_glioma_other","CM2_glioblastoma","UNION_CellMarker2_brain","UNION_all",
        "BrainAtlas112","CARE_TME","Liu_devbrain","CM2_brain_normal"]
MCS_GRID = [15]   # frozen single-sample setting; cohort must match it or the two figures
                  # report different numbers for the same sample (TKU3186: 0.65 at mcs25 vs 0.78 at mcs15)
NOISE = 'Noise'


def macro_f1(pred, gold, classes):
    fs = []
    for c in classes:
        tp = int(((pred == c) & (gold == c)).sum())
        fp = int(((pred == c) & (gold != c)).sum())
        fn = int(((pred != c) & (gold == c)).sum())
        fs.append(2*tp/(2*tp+fp+fn) if (2*tp+fp+fn) else 0.0)
    return float(np.mean(fs))


def ceiling_at(cl, gold, classes):
    """Max macro-F1 over cluster->label assignments, greedy from the majority start.

    A lower bound on the true maximum (single-coordinate moves only); exact enumeration is
    intractable at fine granularity. Reported as such.
    """
    ks = list(pd.unique(cl))
    lab = {k: pd.Series(gold[cl == k]).value_counts().index[0] for k in ks}
    cur = macro_f1(np.array([lab[k] for k in cl], dtype=object), gold, classes)
    for _ in range(6):
        moved = False
        for k in ks:
            keep = lab[k]
            for cand in classes:
                if cand == lab[k]:
                    continue
                lab[k] = cand
                v = macro_f1(np.array([lab[x] for x in cl], dtype=object), gold, classes)
                if v > cur + 1e-12:
                    cur, keep, moved = v, cand, True
                lab[k] = keep
            lab[k] = keep
        if not moved:
            break
    return cur


def compose_lfine(obs):
    """Lfine = L3 for TME cells, Malignant_<MalState> for malignant cells (report Table 2).

    Where BOTH are empty the frozen run falls back to L1, so a cell is never unlabelled.
    Verified against handoff/FROZEN_single_sample_TKU3186: 4384/4384 cells identical.
    """
    bad = {'nan', '', 'None', 'NA'}
    n = len(obs)
    get = lambda k: obs[k].astype(str).values if k in obs else np.array(['nan'] * n)
    L1, L3, MS = get('L1'), get('L3'), get('MalState')
    return np.array([f'Malignant_{m}' if m not in bad else (l3 if l3 not in bad else l1)
                     for l1, l3, m in zip(L1, L3, MS)], dtype=object)


# One-to-many credit map, generalised from handoff/FROZEN_single_sample_TKU3186/report/
# v3_lfine_targets.csv. The pipeline only ever emits the 13-term L1 vocabulary, so it CANNOT
# name an Lfine class exactly; a prediction is credited for every Lfine class it can
# legitimately denote. Cross-family entries are deliberate: a malignant cell in the OPC-like
# state genuinely runs an OPC programme, so predicting OPC on it is not an error.
LFINE_PREFIX = {
    'Malignant':         ('Malignant_',),
    'TAM':               ('TAM_',),
    'Lymphocyte':        ('Lymphocyte', 'TIL_'),
    'Oligodendrocyte':   ('Oligodendrocyte',),
    'OPC':               ('OPC', 'Malignant_OPC'),
    'Astrocyte':         ('Astrocyte', 'Malignant_AC'),
    'Excitatory neuron': ('ExN_', 'Malignant_Neuron'),
    'Inhibitory neuron': ('InN_', 'Malignant_Neuron'),
    'AMBIGUOUS_NEURON':  ('ExN_', 'InN_', 'Malignant_Neuron'),
    'Endothel':          ('Endothelial',),
    'Pericyte':          ('Mural_',),
    'Other':             ('Other',),
}


def lfine_targets(lf_classes):
    """L1 prediction -> set of Lfine classes present in THIS sample that it can denote."""
    return {k: {c for c in lf_classes if any(c.startswith(p) for p in pref)}
            for k, pref in LFINE_PREFIX.items()}


def macro_f1_lfine(pred, lfine, lf_classes, TGT):
    """Report rule, verified to reproduce Table 10 of REPORT_TKU3186 exactly.

    A cell is CORRECT iff its Lfine class is in the target set of its prediction. An INCORRECT
    cell counts as a false positive for EVERY class in that prediction's target set -- so a
    wrong call is charged once per class it wrongly claimed, not once overall.
    """
    ok = np.array([g in TGT.get(q, ()) for g, q in zip(lfine, pred)])
    fs = []
    for c in lf_classes:
        g = (lfine == c)
        has = np.array([c in TGT.get(q, ()) for q in pred])
        tp = int((g & ok).sum()); fn = int((g & ~ok).sum()); fp = int(((~ok) & has & ~g).sum())
        fs.append(2*tp/(2*tp+fp+fn) if (2*tp+fp+fn) else 0.0)
    return float(np.mean(fs)) if fs else np.nan


def score(deg, panels, order, mode):
    names = list(panels)
    S = np.zeros((len(names), len(order)), dtype=np.float32)
    for j, c in enumerate(order):
        d = deg[c]
        for i, p in enumerate(names):
            gs = panels[p]
            if not gs:
                continue
            tot = float(sum(d[g] for g in gs if g in d))
            n = len(gs)
            S[i, j] = (tot / n if mode == 'shipped_div_len' else
                       tot / np.sqrt(n) if mode == 'div_sqrt_len' else
                       tot)                                   # sum_no_norm
    return S, names


def main(samp):
    MAPALL = pd.read_csv(f"{B}/handoff/markers_v3/mapping_L1_v3.csv")
    LIB = {}
    for st in SETS:
        mk = pd.read_csv(f"{B}/handoff/markers_v3/{st}.csv", index_col=0)
        LIB[st] = ({c: [g for g in mk[c].dropna().tolist() if g] for c in mk.columns},
                   MAPALL[MAPALL.marker_set == st].set_index('panel').gold_class.to_dict())

    M = f"{B}/data_bench/GSE274546/mtx/{samp}"
    X = sio.mmread(f"{M}/matrix.mtx").T.tocsr()
    genes = [l.strip() for l in open(f"{M}/genes.tsv")]
    obs = pd.read_csv(f"{M}/obs.csv")
    a = ad.AnnData(X.astype(np.float32)); a.var_names = genes
    a.obs = obs.set_index('CellID'); a.var_names_make_unique()
    sc.pp.filter_genes(a, min_cells=3)
    a.layers['counts'] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    a.layers['lognorm'] = a.X.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    # scale-before-PCA is load-bearing: annotate_v3.py does it, and omitting it changes the
    # embedding enough to alter the cluster count (TKU3186: 4 vs the frozen 6 at mcs=15).
    # The mcs=15 arm must reproduce the frozen run or the grid is not a controlled comparison.
    sc.pp.scale(a, max_value=10); sc.tl.pca(a, n_comps=30, random_state=SEED)
    a.X = a.layers['lognorm'].copy()
    sc.pp.neighbors(a, n_neighbors=15, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)

    gold = a.obs['L1'].astype(str).values
    lfine = compose_lfine(a.obs)
    # Three class-inclusion rules scored from the SAME predictions. Raising the threshold
    # changes what is averaged, not how well anything is annotated; both are reported so the
    # distinction stays visible.
    CLASS_RULES = {'k>=15': 15, 'k>=50': 50, 'k>=100': 100}
    LF_THR = 20          # chosen on the TKU3186 sweep: a plateau (identical to 15), 97% of cells
    classes = [c for c in pd.unique(gold) if (gold == c).sum() >= 15 and c != 'Other']
    lf_vc = pd.Series(lfine).value_counts()
    lf_classes = [c for c in lf_vc.index if lf_vc[c] >= LF_THR and c not in ('Other', 'nan')]
    # Target sets are built over EVERY Lfine class present, not just the scored ones: a cell
    # whose true class fell below LF_THR must still be able to count as CORRECT, otherwise it
    # is charged as a false positive against every scored class the prediction can denote.
    TGT_LF = lfine_targets(list(lf_vc.index))
    if len(classes) < 2:
        pd.DataFrame([dict(sample=samp, note='fewer than 2 scored classes')]).to_csv(
            f"{OUT}/grid_{samp}.csv", index=False)
        return

    rows = []
    for mcs in MCS_GRID:
        lab = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=mcs).fit_predict(a.obsm['X_umap'])
        cl = np.array([NOISE if c == -1 else str(c) for c in lab], dtype=object)
        order = [c for c in pd.unique(cl) if c != NOISE]
        if len(order) < 2:
            rows.append(dict(sample=samp, mcs=mcs, n_clusters=len(order), note='<2 clusters'))
            continue
        a.obs['_cl'] = pd.Categorical(cl)
        sc.tl.rank_genes_groups(a, '_cl', groups=order, method='wilcoxon',
                                n_genes=100, layer='lognorm', use_raw=False)
        deg = {c: dict(zip(a.uns['rank_genes_groups']['names'][c],
                           a.uns['rank_genes_groups']['logfoldchanges'][c])) for c in order}
        ceil = ceiling_at(cl, gold, classes)
        for setname, (panels, mp) in LIB.items():
          for mode in ('shipped_div_len',):   # match the frozen single-sample run
            S, names = score(deg, panels, order, mode)
            call = {}
            for j, c in enumerate(order):
                col = S[:, j]
                call[c] = mp.get(names[int(np.argmax(col))], 'Unknown') if col.max() > 0 else 'Unknown'
            pred = np.array([call.get(k, 'Unknown') for k in cl], dtype=object)
            malg = gold == 'Malignant'
            malp = pred == 'Malignant'
            tp = int((malp & malg).sum()); fp = int((malp & ~malg).sum()); fn = int((~malp & malg).sum())
            for rname, thr in CLASS_RULES.items():
                cls_r = [c for c in pd.unique(gold) if (gold == c).sum() >= thr and c != 'Other']
                if len(cls_r) < 2:
                    continue
                rows.append(dict(sample=samp, mcs=mcs, mode=mode, marker_set=setname, class_rule=rname,
                                 k=len(cls_r), n_clusters=len(order),
                                 macroF1_rule=macro_f1(pred, gold, cls_r),
                                 ceiling_rule=ceiling_at(cl, gold, cls_r)))
            rows.append(dict(
                sample=samp, mcs=mcs, mode=mode, marker_set=setname, n_clusters=len(order),
                noise_frac=float((cl == NOISE).mean()), ceiling=ceil,
                macroF1=macro_f1(pred, gold, classes),
                coverage=float((pred != 'Unknown').mean()),
                acc_on_called=float((pred[pred != 'Unknown'] == gold[pred != 'Unknown']).mean())
                if (pred != 'Unknown').any() else np.nan,
                mal_f1=2*tp/(2*tp+fp+fn) if (2*tp+fp+fn) else np.nan,
                n_classes=len(classes), n_cells=int(a.n_obs),
                # --- Lfine, one-to-many credit (the granularity Dr. Cheng asked for) ---
                lfine_macroF1=macro_f1_lfine(pred, lfine, lf_classes, TGT_LF),
                lfine_n_classes=len(lf_classes),
                lfine_pct_cells=float(np.isin(lfine, lf_classes).mean())))
    pd.DataFrame(rows).to_csv(f"{OUT}/grid_{samp}.csv", index=False)
    print(f"OK {samp}", flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
