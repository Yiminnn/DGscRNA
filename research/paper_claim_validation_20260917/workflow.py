"""Editable English decision tree with explicit completed/pending evidence nodes."""
import os
from common import OUT, require_slurm, checked, write_json, sha, utc

def run():
    require_slurm()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    d=OUT/'summary';d.mkdir(exist_ok=True)
    core=checked(d,'aggregate_manifest.json','AGGREGATE_COMPLETE')
    state='completed' if core else 'running'
    controls='completed' if checked(OUT/'controls_summary') else 'running'
    comparators='completed' if checked(OUT/'comparison_summary') else 'running'
    resources='completed' if checked(OUT/'scalability_summary') else 'running'
    ptc='completed' if checked(OUT/'PTC_summary') else 'pending GBM'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(11,12));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    def node(x,y,w,h,title,body,color='#0072B2'):
        ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=.005',facecolor='#f3f7fa',edgecolor=color,lw=1.2,zorder=3))
        ax.text(x,y+h*.23,title,ha='center',va='center',fontsize=8.4,fontweight='bold',color=color,zorder=4)
        ax.text(x,y-h*.15,body,ha='center',va='center',fontsize=7.3,linespacing=1.25,zorder=4)
    def arrow(x,y,X,Y):ax.add_patch(FancyArrowPatch((x,y),(X,Y),arrowstyle='-|>',mutation_scale=9,color='#687782',lw=1,zorder=1))
    ax.text(.5,.985,'DG-scRNA: native R workflow and evidence map',ha='center',fontsize=15,fontweight='bold')
    ax.text(.5,.965,'Blue = ablation choice  |  Gray = fixed rule  |  Orange = evaluation / evidence scope',ha='center',fontsize=9)
    x=.35;w=.62
    node(x,.917,w,.065,'Inputs and evaluation boundary','GSE274546: 121 samples / 59 patients; primary 97 / 55\nAuthor cells; genes detected in >=3 cells; labels stored separately','#525a61')
    arrow(x,.880,x,.862)
    node(x,.83,w,.06,'A0 Batch scope [existing PTC comparison]','GBM per sample: RNA, no batch integration\nPTC: NMT and TTU CCA; archived all-8 baseline kept separate')
    arrow(x,.795,x,.775)
    node(x,.738,w,.064,f'A1  Feature budget [GBM {state}]','Native VST: 500 / 1,000 / 2,000 / 3,000 / 5,000 / all\nLogNormalize 10,000; center/scale; record every stage gene list')
    arrow(x,.701,x,.684)
    node(x,.66,w,.04,'PCA30 [original default]','HVG geometry; seed 42','#525a61')
    arrow(x,.635,.17,.613);arrow(x,.635,.52,.613)
    node(.17,.586,.29,.047,'A2  PCA30','Retain 30 PCs for clustering')
    node(.52,.586,.29,.047,'A2  PCA30 -> UMAP2','uwot cosine; neighbors 30; min_dist .3')
    for u,v in [(.17,.09),(.17,.26),(.52,.44),(.52,.61)]:arrow(u,.557,v,.541)
    for xx,title,body in [(.09,'A3 SNN','Louvain r=.5'),(.26,'A3 HDBSCAN','R minPts=50'),(.44,'A3 SNN','Louvain r=.5'),(.61,'A3 HDBSCAN','R minPts=50')]:
        node(xx,.515,.15,.047,title,body);arrow(xx,.488,x,.47)
    node(x,.44,w,.053,'Cluster DEG and original density score','Wilcoxon; density sums DEG log2FC>1, full-panel denominator\nSingleton factor .8; ties -> Undecided; HDBSCAN noise scored as cluster','#525a61')
    arrow(x,.41,x,.394)
    node(x,.363,w,.055,f'A4  Marker context [GBM {state}]','16 frozen libraries; 3 cutoffs: none / mean / 0.5\nBrain / glioma / immune / vascular / unions / AllHuman / curated references')
    arrow(x,.33,x,.313)
    node(x,.277,w,.063,f'A5  Marker-only vs terminal DL [GBM {state}]','Native panel labels -> normalized selected-gene MLP input\n256/128; LeakyReLU; legacy Softmax+CE; Adamax; 10 epochs\nFill only Undecided; preserve trained / no-op / unavailable state')
    arrow(x,.239,x,.225)
    node(x,.202,w,.04,'A6  Confidence threshold','0.90 primary; 0.70 from the same saved probabilities')
    arrow(x,.178,x,.160)
    node(x,.123,w,.064,'Frozen patient-based evaluation','All cells scored; Unknown / unmapped count as errors; coverage separate\nAuthor L1 macro-F1; per-class confusion; ARI/NMI only for clustering\nPatient-heldout configuration selection; paired patient uncertainty','#C46518')
    ax.axvline(.7,ymin=.07,ymax=.95,color='#d4dbe0',lw=1)
    notes=[(.94,'Scope of the core comparison',
        '6 budgets x 4 clustering routes\nx 16 libraries x 3 cutoffs.\nOriginal 2,000/UMAP/HDBSCAN\nremains a fixed anchor.\nA high score is not assumed.'),
        (.79,'What HVG changes',
        'Single-sample GBM RNA:\ngeometry + DL genes vary;\nDEG scoring keeps all RNA genes.\nPTC CCA budget: anchors,\ngeometry, scoring and DL vary.\nA geometry-only arm must fix DL\nand the scoring gene universe.'),
        (.60,'Ablation completion',
        f'A7 Geometry-only: {controls}\nA8 Dimensions / clustering /\nembedding seeds: {controls}\nA9 MLP width / epochs /\nmodel seeds: {controls}\nPTC retention: {ptc}'),
        (.43,'Reference interpretation',
        'CARE_TME / BrainAtlas112\ncontributed to author labels.\nTheir results are concordance,\nnot independent validation.\nNo generic neuron is assigned\nto an excitatory/inhibitory class\nusing evaluation outcomes.'),
        (.245,'Publication claim boundary',
        f'Best among evaluated choices,\nwithin specified cohorts/metrics.\nHistorical NMT selected SNN.\nInternal ablation does not prove\nsuperiority to other tools.\nFair comparators: {comparators}\nResource repeats (3): {resources}')]
    for y,title,body in notes:
        ax.text(.725,y,title,fontsize=9,color='#C46518',weight='bold',va='top')
        ax.text(.725,y-.025,body,fontsize=8,va='top',linespacing=1.4)
    ax.text(.04,.047,'Retrospective validation: frozen patient folds isolate test labels during configuration selection. GBM fits each sample separately;\nPTC integrates full-group expression, with patient labels held out for selection. No global-optimum or never-seen-cohort claim.',fontsize=8,va='top')
    fig.subplots_adjust(left=.01,right=.99,top=.99,bottom=.01)
    for ext in ['png','pdf','svg']:fig.savefig(d/f'workflow_decision_tree.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    write_json(d/'workflow_manifest.json',dict(GBM_core_status=state,job=os.environ['SLURM_JOB_ID'],
        source_sha256=sha(__file__),completed_at=utc(),followups_separately_tracked=True,
        controls=controls,comparators=comparators,scalability=resources,PTC=ptc))

if __name__=='__main__':run()
