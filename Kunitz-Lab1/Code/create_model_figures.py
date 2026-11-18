from pathlib import Path
import csv,json,math,hashlib,sys,os
sys.path.insert(0,str(Path('work/report_packages').resolve()))
os.environ['MPLCONFIGDIR']=str(Path('work/mplconfig').resolve())
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,PathPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D
R=Path('work/raw_audit');B=R/'results/hmm/build_20261107T212501_925161Z';O=Path('work/model_evidence_figures');O.mkdir(exist_ok=True)
record=json.loads((B/'build_record.json').read_text());hmm=B/'kunitz_initial.hmm';sto=B/'kunitz_seed.sto'
assert hashlib.sha256(hmm.read_bytes()).hexdigest()==record['model_sha256']
assert hashlib.sha256(sto.read_bytes()).hexdigest()==record['input_sha256']
seqs={}
for line in sto.read_text().splitlines():
 if not line.strip() or line.startswith(('#','//')):continue
 key,seq=line.split()[:2];seqs[key]=seqs.get(key,'')+seq
assert len(seqs)==13 and {len(x) for x in seqs.values()}=={53}
seed={};key=None
for l in (R/'data/processed/PF00014_structural_seed.fasta').read_text().splitlines():
 if l.startswith('>'):key=l[1:].split()[0];seed[key]=''
 else:seed[key]+=l.strip()
assert all(s.replace('-','').replace('.','')==seed[k] for k,s in seqs.items())
with (R/'data/processed/PF00014_fragment_manifest.csv').open() as f:manifest={r['seed_id']:r for r in csv.DictReader(f)}
conserved_c=[i+1 for i in range(53) if all(s[i]=='C' for s in seqs.values())]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
# Display actual structural-alignment-derived amino-acid rows.
fig,ax=plt.subplots(figsize=(12,3.75));ax.set_xlim(-18,53.6);ax.set_ylim(-1.1,13.6)
colors={'C':'#f5d779','D':'#efbbbb','E':'#efbbbb','K':'#c1d8ee','R':'#c1d8ee','H':'#c1d8ee','F':'#d4e3ce','Y':'#d4e3ce','W':'#d4e3ce','I':'#d4e3ce','L':'#d4e3ce','V':'#d4e3ce','M':'#d4e3ce','A':'#e5ecdf','G':'#eadff0','P':'#eadff0','S':'#dce9ee','T':'#dce9ee','N':'#dce9ee','Q':'#dce9ee','-':'#eeeeee'}
for row,(k,s) in enumerate(seqs.items()):
 y=12-row;m=manifest[k];label=f"{m['entity_identifier'].split('_')[0]}:{m['source_label_chain']} {m['start']}-{m['end']}"
 ax.text(-.7,y+.47,label,ha='right',va='center',fontsize=11.5,fontfamily='DejaVu Sans Mono')
 for i,a in enumerate(s):
  ax.add_patch(Rectangle((i,y),.98,.94,facecolor=colors.get(a,'#f6f6f6'),edgecolor='white',lw=.2))
  ax.text(i+.49,y+.47,a,ha='center',va='center',fontsize=12,fontfamily='DejaVu Sans Mono',weight='bold' if a=='C' else 'normal')
for i in [1,5,10,15,20,25,30,35,40,45,50,53]:ax.text(i-.5,13.25,str(i),ha='center',fontsize=10)
ax.text(-.7,13.25,'PDB chain / residues',ha='right',fontsize=10)
ax.text(26,-.8,'Structural-alignment column',ha='center',fontsize=11);ax.axis('off');fig.subplots_adjust(left=.01,right=.995,top=.98,bottom=.03)
for ext in ['png','pdf']:fig.savefig(O/f'structural_alignment.{ext}',dpi=350,bbox_inches='tight',pad_inches=.05)
plt.close(fig)
# Parse model match emissions: HMMER stores negative natural log probabilities.
alpha=None;probs=[];states=[]
for line in hmm.read_text().splitlines():
 fields=line.split()
 if line.startswith('HMM '):alpha=fields[1:]
 if fields and fields[0].isdigit() and len(fields)>=21:
  states.append(int(fields[0]));probs.append([0 if v=='*' else math.exp(-float(v)) for v in fields[1:21]])
assert states==list(range(1,54)) and len(alpha)==20
P=np.array(probs);assert np.max(abs(P.sum(axis=1)-1))<1e-4
P/=P.sum(axis=1,keepdims=True)
H=-(P*np.log2(np.clip(P,1e-300,None))).sum(axis=1);I=np.log2(20)-H;heights=P*I[:,None]
aa_colors={a:'#376d91' for a in alpha}
for a in 'DE':aa_colors[a]='#bd4d44'
for a in 'KRH':aa_colors[a]='#3b69ad'
for a in 'STNQ':aa_colors[a]='#438b75'
for a in 'AGP':aa_colors[a]='#777777'
aa_colors['C']='#bf9224'
font=FontProperties(family='DejaVu Sans',weight='bold');glyphs={a:TextPath((0,0),a,size=1,prop=font) for a in alpha}
fig,ax=plt.subplots(figsize=(12,2.9))
for i,row in enumerate(heights):
 bottom=0
 for j in np.argsort(row):
  ht=float(row[j]);a=alpha[j];path=glyphs[a];box=path.get_extents()
  trans=Affine2D().translate(-box.x0,-box.y0).scale(.88/box.width,ht/box.height).translate(i+1-.44,bottom)
  ax.add_patch(PathPatch(path,transform=trans+ax.transData,color=aa_colors[a],lw=0));bottom+=ht
ax.set(xlim=(.3,53.7),ylim=(0,max(3.2,float(I.max())+.15)),xlabel='HMM match state',ylabel='Information (bits)')
ax.set_xticks([1,5,10,15,20,25,30,35,40,45,50,53]);ax.tick_params(labelsize=10);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(O/f'hmm_logo.{ext}',dpi=350,bbox_inches='tight',pad_inches=.05)
plt.close(fig)
# Counts of actual recalibrated predictions; no normalization by background size.
f=json.loads(Path('work/reanalysis/full_scan_reanalysis.json').read_text());counts=f['taxonomic_counts'];order=['Vertebrates','Arthropods','Cnidarians','Other metazoans','Nematodes'];values=[counts[k] for k in order];assert sum(values)==398
fig,ax=plt.subplots(figsize=(5,3));bars=ax.barh(order[::-1],values[::-1],color='#376d91',height=.65)
for bar,v in zip(bars,values[::-1]):ax.text(v+3,bar.get_y()+bar.get_height()/2,str(v),va='center',fontsize=11)
ax.set(xlim=(0,260),xlabel='Predicted proteins');ax.grid(axis='x',alpha=.17);ax.set_axisbelow(True);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(O/f'taxonomic_distribution.{ext}',dpi=350,bbox_inches='tight',pad_inches=.05)
plt.close(fig)
with (O/'hmm_match_emissions.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['state']+alpha+['information_bits_uniform_background'])
 for state,row,info in zip(states,P,I):w.writerow([state]+list(row)+[info])
provenance={'model_sha256':hashlib.sha256(hmm.read_bytes()).hexdigest(),'alignment_sha256':hashlib.sha256(sto.read_bytes()).hexdigest(),'alignment_rows':13,'alignment_columns':53,'universally_conserved_cysteine_columns':conserved_c,'logo':'Match-emission probabilities exp(-stored_cost), row normalized for file rounding; stack height log2(20)-Shannon entropy relative to uniform amino-acid background; letter height probability times stack height; no finite-sample correction','taxonomic_total':sum(values),'threshold_bits':17.2,'model_threshold_results_changed':False}
(O/'figure_provenance.json').write_text(json.dumps(provenance,indent=2));print(json.dumps(provenance,indent=2))
