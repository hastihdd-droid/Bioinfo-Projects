from pathlib import Path
import json,csv,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
BASE=Path(__file__).resolve().parent;O=BASE/'figures';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
rows=list(csv.DictReader((BASE/'results/training_quality_summary.csv').open()))
fig,ax=plt.subplots(figsize=(7.1,3.7))
for method,color,marker,label in [('X-RAY DIFFRACTION','#277f8e','o','X-ray (10 fragments)'),('SOLUTION NMR','#b46a35','^','NMR (3 fragments)')]:
 r=[r for r in rows if r['method']==method];ax.scatter([float(x['mean_pairwise_rmsd_angstrom']) for x in r],[float(x['mean_pairwise_Q']) for x in r],s=55,c=color,marker=marker,label=label,edgecolor='white',linewidth=.7,zorder=3)
for name,offset in [('6NAN',(7,5)),('1IRH',(6,5)),('1KTH',(-42,-8)),('4U32',(8,2))]:
 r=next(r for r in rows if r['pdb']==name);ax.annotate(name,(float(r['mean_pairwise_rmsd_angstrom']),float(r['mean_pairwise_Q'])),xytext=offset,textcoords='offset points',fontsize=9)
ax.set(xlabel='Mean pairwise fitted Cα RMSD (Å)',ylabel='Mean pairwise Q score',xlim=(.65,1.8),ylim=(.70,.96));ax.grid(alpha=.18);ax.legend(frameon=False,loc='upper right',fontsize=9);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(O/f'training_structural_quality.{ext}',dpi=320,bbox_inches='tight')
plt.close(fig)
records=json.loads((BASE/'results/chimerax_results.json').read_text());assert all(r['status']=='success' for r in records)
ref=np.array([x['coord'] for x in records[0]['reference_trace']]);center=ref.mean(0);_,_,vt=np.linalg.svd(ref-center,full_matrices=False);rot=vt.T
fig=plt.figure(figsize=(12,4.7));axes=[]
for j,r in enumerate(records):
 ax=fig.add_subplot(1,3,j+1,projection='3d');axes.append(ax)
 reference={x['residue']:np.array(x['coord']) for x in r['reference_trace']};target={x['residue']:np.array(x['coord']) for x in r['target_trace']};core={x['target_residue'] for x in r['mapping'] if x['is_pruned_core']};mapped={x['target_residue'] for x in r['mapping']};lo,hi=min(mapped),max(mapped)
 def transform(x):return (np.array(x)-center)@rot
 def segments(trace,ids):
  ids=sorted(ids);return np.array([[transform(trace[a]),transform(trace[b])] for a,b in zip(ids[:-1],ids[1:]) if b==a+1])
 refseg=segments(reference,set(reference));tarseg=segments(target,{n for n in target if lo<=n<=hi});coreseg=segments(target,core)
 ax.add_collection3d(Line3DCollection(refseg,colors='#bf8b2e',linewidths=2.5,alpha=.95))
 ax.add_collection3d(Line3DCollection(tarseg,colors='#a6c5d2',linewidths=1.2,alpha=.85))
 if len(coreseg):ax.add_collection3d(Line3DCollection(coreseg,colors='#286a8c',linewidths=2.6))
 cxyz=np.array([transform(target[n]) for n in core]);ax.scatter(*cxyz.T,s=5,c='#286a8c',depthshade=False)
 points=np.vstack([transform(list(reference.values())),transform([target[n] for n in target if lo<=n<=hi])]);mid=(points.min(0)+points.max(0))/2;span=max(np.ptp(points,axis=0))*1.06
 for dim,setlim in enumerate([ax.set_xlim,ax.set_ylim,ax.set_zlim]):setlim(mid[dim]-span/2,mid[dim]+span/2)
 ax.set_box_aspect((1,1,1));ax.view_init(elev=18,azim=-68);ax.set_axis_off();ax.set_title(r['Entry'],fontsize=18,weight='bold',pad=1)
 ax.text2D(.5,-.01,f"Core: {r['core_rmsd_angstrom']:.3f} Å / {r['core_pairs']} pairs\nAll pairs: {r['full_rmsd_angstrom']:.3f} Å (n={r['full_pairs']})\nMean core pLDDT: {r['mean_core_pLDDT']:.1f}",transform=ax.transAxes,ha='center',va='top',fontsize=15)
fig.legend(handles=[Line2D([0],[0],color='#bf8b2e',lw=2.5,label='Experimental BPTI (3TGI:I)'),Line2D([0],[0],color='#286a8c',lw=2.5,label='AlphaFold matched core'),Line2D([0],[0],color='#a6c5d2',lw=2,label='Other residues in the matched region')],loc='lower center',ncol=3,frameon=False,fontsize=14,bbox_to_anchor=(.5,-.04))
fig.subplots_adjust(left=.01,right=.99,top=.92,bottom=.23,wspace=.02)
for ext in ['png','pdf']:fig.savefig(O/f'false_negative_superpositions.{ext}',dpi=320,bbox_inches='tight')
plt.close(fig)
# Independent coordinate and sequence checks on the actual transformations.
checks=[]
for r in records:
 ref={x['residue']:np.array(x['coord']) for x in r['reference_trace']};tar={x['residue']:np.array(x['coord']) for x in r['target_trace']}
 full=[np.sum((ref[x['reference_residue']]-tar[x['target_residue']])**2) for x in r['mapping']];core=[v for v,x in zip(full,r['mapping']) if x['is_pruned_core']]
 fr=float(np.sqrt(np.mean(full)));cr=float(np.sqrt(np.mean(core)));assert abs(fr-r['full_rmsd_angstrom'])<1e-4;assert abs(cr-r['core_rmsd_angstrom'])<1e-4
 checks.append(dict(Entry=r['Entry'],independent_full_rmsd=fr,independent_core_rmsd=cr,core_target_residue_min=min(x['target_residue'] for x in r['mapping'] if x['is_pruned_core']),core_target_residue_max=max(x['target_residue'] for x in r['mapping'] if x['is_pruned_core'])))
(BASE/'results/coordinate_verification.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
