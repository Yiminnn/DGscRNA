"""Editable full-workflow decision tree; run rendering inside SLURM."""
import os,json
from pathlib import Path

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
    root=Path('/fs/scratch/PCON0080/yimin/dgscrna')
    out=root/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary'
    out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(8.3,10.7));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    blue='#0072B2';orange='#D55E00';gray='#50565B'
    def node(x,y,w,h,title,body='',kind='ablation',size=8):
        color=blue if kind=='ablation' else gray
        face='#EEF7FC' if kind=='ablation' else '#F3F4F5'
        ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.006,rounding_size=0.008',
                     facecolor=face,edgecolor=color,lw=1.15,zorder=3))
        ax.text(x,y+(h*.19 if body else 0),title,ha='center',va='center',fontsize=size,fontweight='bold',color=color,zorder=4)
        if body:ax.text(x,y-h*.16,body,ha='center',va='center',fontsize=size-1,color='#252525',zorder=4,linespacing=1.25)
    def arrow(x1,y1,x2,y2):
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=9,color='#65717A',lw=1,zorder=1))
    ax.text(.5,.98,'DG-scRNA: R-reference workflow and ablation map',ha='center',va='top',fontsize=13,fontweight='bold')
    ax.text(.5,.954,'Blue: varied in the experiment grid   |   Gray: fixed algorithm or evaluation rule',ha='center',fontsize=8)
    x=.345;w=.61
    node(x,.909,w,.056,'Input and normalization','Fixed published cells; genes detected in ≥3 cells; LogNormalize 10,000',kind='fixed')
    arrow(x,.875,x,.860)
    node(x,.829,w,.060,'A1  Integration scope','PTC: archived / fresh all 8 versus NMT and TTU separately')
    arrow(x,.793,.18,.774);arrow(x,.793,.51,.774)
    node(.18,.737,.28,.066,'A2  Feature budget','500 / 1k / 2k / 3k / 5k / all\nRecord anchor, geometry, score, DL genes',size=7.9)
    node(.51,.737,.28,.066,'A3  Batch correction','CCA expression / RNA / Harmony PCs\nRNA vs Harmony: matched score + DL',size=7.9)
    arrow(.18,.698,x,.694);arrow(.51,.698,x,.694)
    node(x,.679,w,.021,'Paired choices → expression → fixed scaling → PCA (30 PCs)',kind='fixed',size=7.2)
    arrow(x,.662,.18,.650)
    node(.18,.622,.28,.044,'A4  PCA','Use the 30 saved PCs')
    node(.51,.622,.28,.044,'A4  UMAP','2D from those PCs; seed fixed')
    arrow(.329,.622,.361,.622)
    for start,end in [(.18,.09),(.18,.265),(.51,.425),(.51,.60)]:arrow(start,.594,end,.590)
    for xx,title,body in [(.09,'A5  SNN','PCA graph'),(.265,'A5  HDBSCAN','PCA distance'),(.425,'A5  SNN','UMAP graph'),(.60,'A5  HDBSCAN','UMAP distance')]:
        node(xx,.562,.15,.044,title,body,size=7.4)
        arrow(xx,.534,x,.506)
    node(x,.478,w,.054,'Cluster DEG statistics','Seurat v4-compatible Wilcoxon; logFC >1 for density scoring',kind='fixed')
    arrow(x,.445,x,.430)
    node(x,.399,w,.060,'A6  Marker context','PTC: 17 original libraries; reviewer data: prespecified CellMarker tissues\nPrimary / disease / immune / stromal / relevant organs / unions / all')
    arrow(x,.362,x,.347)
    node(x,.321,w,.049,'A7  Density cutoff','none / mean / 0.5; full marker denominators and tie rule fixed')
    arrow(x,.290,x,.275)
    node(x,.247,w,.053,'A8  DL refinement','Marker-only ablation → original-style MLP; fill Undecided only')
    arrow(x,.215,x,.200)
    node(x,.174,w,.048,'A9  Confidence threshold','0.90 primary; 0.70 sensitivity; record no-op / untrainable states')
    arrow(x,.144,x,.129)
    node(x,.093,w,.063,'Evaluation after terminal DL','Reviewer curated labels; PTC TCR agreement and archived-label concordance\nClustering ARI/NMI is a separate endpoint; retain Unknowns',kind='fixed')
    # Side notes specify the scope of each ablation rather than imply a full factorial.
    ax.axvline(.698,ymin=.055,ymax=.935,color='#CDD2D6',lw=.8)
    notes=[(.913,'Cohort scope',
      'PTC: 92,404 historical cells.\nNMT: MT-1/2 + N-1/2.\nTTU: TU-1/2 + T-1/2.\n\n11 reviewer datasets.\nHCL: 59 original tissue groups;\nall 599,926 cells retained.\nHCL is a conditional tissue analysis.'),
      (.721,'Two distinct HVG contrasts',
      '1. Full CCA feature budget:\nanchor + geometry + score + DL\nchange together.\n\n2. Geometry-only control:\nfixed all-gene CCA expression;\nfixed 2,000 score and DL genes.\nOnly geometry genes change.'),
      (.514,'Algorithm parameters held fixed',
      'SNN: Louvain, resolution 0.5.\nHDBSCAN: R minPts = 50.\nUMAP: cosine, 30 neighbors,\nmin_dist = 0.3, seed = 42.\n\nDEG: original R-compatible rule.\nDensity: singleton factor 0.8;\nfull panel denominator;\nties → Undecided.'),
      (.303,'DL parameters held fixed',
      'MLP: 256 / 128, LeakyReLU.\nLegacy Softmax → CrossEntropy.\nAdamax 0.001; 10 epochs.\n90/10 split; batch 256.\nOnly unresolved cells updated.\n\nNo evaluation labels enter fitting.'),
      (.126,'Interpretation',
      'Best among measured conditions\nis dataset- and endpoint-specific.\nThe grid is not fully factorial.\nNo global-optimum claim.\nSaved labels are concordance,\nnot independent ground truth.')]
    for y,title,body in notes:
        ax.text(.72,y,title,ha='left',va='top',fontsize=8.4,fontweight='bold',color=orange)
        ax.text(.72,y-.023,body,ha='left',va='top',fontsize=7.1,linespacing=1.3)
    ax.text(.04,.026,'PTC varies A1–A9; reviewer cohorts use the reference 2,000-gene preparation and vary A4–A9.\n“All genes” means genes shared after the fixed detection filter; no variance ranking. Adaptations for tiny batches are logged.',
            fontsize=7,ha='left',va='center',linespacing=1.3)
    fig.subplots_adjust(left=.025,right=.99,top=.995,bottom=.015)
    for ext in ['png','pdf','svg']:fig.savefig(out/f'workflow_decision_tree.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    (out/'workflow_figure_manifest.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'formats':['SVG editable text','PDF','PNG 300 dpi'],'ablation_nodes':list(range(1,10))},indent=2)+'\n')

if __name__=='__main__':run()
