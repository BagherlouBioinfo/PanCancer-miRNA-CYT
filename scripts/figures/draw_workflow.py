from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Polygon,FancyArrowPatch
import argparse
parser=argparse.ArgumentParser(description="Render the original study workflow without changing its layout")
parser.add_argument("--output-dir",required=True)
args=parser.parse_args()
O=Path(args.output_dir).resolve()
O.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.8,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(6,6.0));fig.subplots_adjust(0,0,1,1);ax.set(xlim=(-4,164),ylim=(10,186));ax.axis('off');ink='#223443';edge='#647787';pale='#F5F7F9';blue='#E6EFF7';green='#EDF4F0'
def box(x,y,w,h,head,body,fill=pale,fs=8.8):
 ax.add_patch(Rectangle((x,y),w,h,facecolor=fill,edgecolor=edge,lw=.85,zorder=3))
 ax.text(x+w/2,y+h-2.2,head,ha='center',va='top',fontsize=9.5,fontweight='bold',color=ink,zorder=4)
 ax.text(x+w/2,y+h-7,body,ha='center',va='top',fontsize=fs,color=ink,linespacing=1.28,zorder=4)
def arrow(a,b,label=None,xy=None):
 ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=8,lw=.85,color=edge,zorder=2,shrinkA=0,shrinkB=0))
 if label:ax.text(*xy,label,fontsize=8.8,ha='left',va='center',color=ink,zorder=5,bbox=dict(facecolor='white',edgecolor='none',pad=.4))
def line(xs,ys):ax.plot(xs,ys,color=edge,lw=.85,zorder=1)
box(3,163,115,16,'Matched TCGA cohort','9,286 primary tumours · 31 cancer types\nMatched miRNA, mRNA and clinical data')
arrow((60.5,163),(60.5,158))
box(3,142,115,16,'Within-cancer discovery and TargetScan integration','miRNA–CYT · miRNA–gene · gene–CYT\nPrimary CPE analysis: 7,523 tumours · 20 cancers')
arrow((60.5,142),(60.5,138))
ax.add_patch(Polygon([(60.5,138),(118,126),(60.5,114),(3,126)],facecolor='white',edgecolor=edge,lw=.9,zorder=3))
ax.text(60.5,126,'All three arms: |ρ| ≥ 0.30,\nBH-FDR < 0.05 in the same cancer,\nwith recurrence in ≥3 cancers?',ha='center',va='center',fontsize=8.8,color=ink,linespacing=1.22,zorder=4)
arrow((118,126),(128,126),'No',(120,129.2))
ax.add_patch(Rectangle((128,119),29,14,facecolor=pale,edgecolor=edge,lw=.85,zorder=3));ax.text(142.5,126,'Not in the\nprimary set',ha='center',va='center',fontsize=8.8,color=ink,zorder=4)
arrow((60.5,114),(60.5,108),'Yes',(63,111))
box(3,92,115,16,'Fixed primary association set','47 pairs · 6 miRNAs · 40 genes\n315 qualifying pair–cancer records',fill=blue)
# Four direct branches from the fixed set; no downstream analysis redefines selection.
line([60.5,60.5,80],[92,88,88]);line([40.5,119.5],[88,88]);line([80,80],[88,48])
for x in [40.5,119.5]:arrow((x,88),(x,83))
line([40.5,119.5],[48,48]);arrow((40.5,48),(40.5,43));arrow((119.5,48),(119.5,43))
box(3,58,75,25,'Sensitivity and specificity testing','LUMP purity sensitivity\nImmune-population adjustment\nCorrelation-matched non-target controls',fs=8.7)
box(82,58,75,25,'Pooled and exploratory analyses','Random-effects meta-analysis\nFeature and cluster survival\nMediation · Clustering · miR-146a/b\nEnrichment and drug annotation',fs=8.7)
box(3,19,75,24,'Cellular-context assessment','Rank-based interactions and HC3 inference\nConjunction and sequential checks\n31 contexts representing 19 pairs',fs=8.7)
box(82,19,75,24,'Independent CPTAC assessment','Fixed pairs and thresholds; no reselection\n101 matched breast tumours\n96 with WES-derived purity',fill=green,fs=8.7)
# Arrows represent workflow dependencies, not causal paths.
for ext in ['pdf','svg']:fig.savefig(O/f'Figure_1.{ext}',facecolor='white')
fig.savefig(O/'Figure_1.png',dpi=600,facecolor='white');fig.savefig(O/'preview.png',dpi=180,facecolor='white');plt.close(fig)
