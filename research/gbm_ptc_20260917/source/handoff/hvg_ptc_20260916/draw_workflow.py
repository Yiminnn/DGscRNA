#!/usr/bin/env python3
"""Publication vector workflow with explicit ablation and terminal-state nodes."""
from common import OUT,require_slurm
require_slurm()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Polygon,FancyArrowPatch

fig,ax=plt.subplots(figsize=(16,10))
ax.set_xlim(0,20.6);ax.set_ylim(-.25,13.9);ax.axis('off')
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
palette={'normal':('#F1F4F6','#52616B'),'ablation':('#FFF0D5','#BC7414'),
         'final':('#DFF1E8','#187D5B'),'stop':('#F9E7E3','#AA6255'),'reference':('#E7F1FA','#3676A5')}

def box(x,y,w,h,title,body,kind='normal',size=10):
    fill,edge=palette[kind]
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.06,rounding_size=0.14',
                              linewidth=1.4,edgecolor=edge,facecolor=fill,zorder=3))
    ax.text(x,y+h/2-.23,title,ha='center',va='top',fontsize=size+1,weight='bold',color='#172D3B',zorder=4)
    ax.text(x,y+h/2-.62,body,ha='center',va='top',fontsize=size,linespacing=1.35,color='#263942',zorder=4)

def diamond(x,y,w,h,text,kind='normal'):
    fill,edge=palette[kind]
    ax.add_patch(Polygon([(x,y+h/2),(x+w/2,y),(x,y-h/2),(x-w/2,y)],closed=True,
                         edgecolor=edge,facecolor=fill,lw=1.4,zorder=3))
    ax.text(x,y,text,ha='center',va='center',fontsize=10.4,color='#263942',linespacing=1.3,zorder=4)

def arrow(points,label=None,at=None,color='#637780',style='solid'):
    for a,b in zip(points[:-2],points[1:-1]):
        ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=1.25,ls=style,zorder=1)
    ax.add_patch(FancyArrowPatch(points[-2],points[-1],arrowstyle='-|>',mutation_scale=12,
                                lw=1.25,color=color,linestyle=style,zorder=2))
    if label:
        ax.text(*at,label,ha='center',va='center',fontsize=9.4,color=color,
                bbox={'facecolor':'white','edgecolor':'none','pad':1.4},zorder=5)

ax.text(.25,13.6,'DG-scRNA: from controlled feature ablations to terminal annotation',fontsize=19,weight='bold',color='#172D3B')
ax.text(.25,13.12,'A  |  Frozen cells, controlled geometry and clustering',fontsize=12,weight='bold',color='#172D3B')
box(2.45,11.55,4.15,2.05,'Input and preprocessing','121 GBM samples / 59 patients\nAuthor-QC cells; genes in ≥3 cells\nNormalize 10k → log1p → scale',kind='reference')
box(7.45,11.55,4.15,2.05,'A1  |  Feature selection','All / 500 / 1k / 2k / 3k / 5k HVG\nVST on counts; legacy flavor bridge\nA5: HVG2000 ∪ marker genes',kind='ablation')
box(12.45,11.55,4.15,2.05,'A2–A3  |  Representation','Direct input vs PCA30 → UMAP\nSix reducers (2D), or no reduction\nUMAP 2 / 10 / 30D; PCA30 control',kind='ablation')
box(17.45,11.55,4.15,2.05,'A4  |  Clustering','KMeans / diagonal GMM (K=23)\nHDBSCAN (15/15)\nSeparate K and density sensitivity',kind='ablation')
for x in [2.45,7.45,12.45]:arrow([(x+2.14,11.55),(x+2.86,11.55)])
box(2.45,8.95,4.15,1.35,'Fixed annotation features','Full log-normalized + scaled genes\nA5: explicit HVG-only scoring contrast',kind='ablation',size=9.5)
diamond(12.45,8.95,3.75,1.45,'≥2 non-noise\nclusters?')
box(17.45,8.95,3.85,1.35,'Structural abstention','Terminal Unknown\nNo fitted DL claim',kind='stop',size=9.5)
arrow([(2.45,10.45),(2.45,9.68)])
arrow([(17.45,10.45),(17.45,10.15),(12.45,10.15),(12.45,9.72)])
arrow([(14.4,8.95),(15.45,8.95)],'No',(14.95,9.18))
ax.text(5.4,7.62,'B  |  Terminal marker and DL/refinement decisions',fontsize=12,weight='bold',color='#172D3B')
box(2.9,6.1,4.55,1.95,'Cluster marker annotation','Top 100 Wilcoxon DEGs; legacy scores\nCM2_glioma_other; fixed L1 mapping\nPositive unique winner or Undecided\nUnsupported singleton DEG → unavailable',size=9.1)
diamond(7.55,6.1,3.3,1.6,'Non-noise\nundecided pool?')
diamond(12.45,6.1,3.4,1.6,'≥2 seed\ntraining classes?')
box(17.45,6.1,3.85,1.45,'Structurally untrainable','Known calls retained\nPool Unknown; no fitted DL',kind='stop',size=9.5)
arrow([(12.45,8.2),(12.45,7.94),(2.9,7.94),(2.9,7.14)],'Yes',(10,8.13))
arrow([(2.45,8.2),(2.45,7.94),(2.9,7.94)],color='#BC7414')
arrow([(5.25,6.1),(5.82,6.1)])
arrow([(9.26,6.1),(10.68,6.1)],'Yes',(9.98,6.35))
arrow([(14.22,6.1),(15.45,6.1)],'No',(14.86,6.36))
box(7.55,3.65,3.55,1.48,'Empty-pool no-op','Known calls retained\nNoise remains Unknown',kind='final',size=9.6)
box(12.45,3.65,4.25,1.7,'A6  |  DL/refinement','Fixed full-gene input; 15 epochs\nFill pool if probability ≥0.9\nNoise excluded; known calls unchanged',kind='ablation',size=9.6)
arrow([(7.55,5.23),(7.55,4.47)],'No',(7.81,4.88))
arrow([(12.45,5.23),(12.45,4.58)],'Yes',(12.77,4.9))
box(17.45,1.05,4.25,1.82,'Terminal DG-scRNA output','Saved per-cell final calls + lineage\nExecuted / no-op / untrainable status\nTechnical or structural missingness explicit',kind='final',size=9.1)
arrow([(17.45,8.19),(19.94,8.19),(19.94,1.05),(19.65,1.05)])
arrow([(17.45,5.31),(17.45,2.03)])
arrow([(12.45,2.73),(12.45,2.37),(16.1,2.37),(16.1,2.03)])
arrow([(7.55,2.84),(7.55,2.49),(15.63,2.49),(15.63,2.03)])
box(10.05,.75,6.45,1.48,'Evaluation after fitting','Strict L1 + legacy set-valued concordance; separate partition metrics\nPatient-paired effects, availability, coverage, noise and bounds',kind='reference',size=9.2)
box(3.5,.75,4.75,1.48,'Author reference labels','97 samples / 55 patients: primary analysis\nAll 121 samples: reported separately',kind='reference',size=9.3)
arrow([(15.26,.75),(13.35,.75)])
arrow([(5.95,.75),(6.75,.75)],color='#3676A5',style='dashed')
ax.text(.4,3.88,'Orange nodes: named ablations',fontsize=10.8,weight='bold',color='#96611E')
ax.text(.4,3.48,'Five fixed seeds and UMAP parameter sensitivities.\nMarker-only calls appear only in the A6 ablation.\nHVG changes geometry; full genes stay available\nfor scoring and DL. Reference labels never enter fitting.',fontsize=9.6,va='top',linespacing=1.45,color='#52616B')
fig.subplots_adjust(left=.025,right=.995,top=.98,bottom=.04)
out=OUT/'figures';out.mkdir(parents=True,exist_ok=True)
for ext in ['pdf','svg','png']:
    fig.savefig(out/f'decision_tree_main.{ext}',dpi=180,bbox_inches='tight',facecolor='white')
plt.close(fig)
