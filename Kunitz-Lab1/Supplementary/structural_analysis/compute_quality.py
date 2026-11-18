from pathlib import Path
import csv,json,itertools,hashlib
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
BASE=Path(__file__).resolve().parent;I=BASE/'inputs';O=BASE/'results';O.mkdir(exist_ok=True)
seqs={};key=None
for l in (I/'seed_alignment.fasta').read_text().splitlines():
 if l.startswith('>'):key=l[1:].split()[0];seqs[key]=''
 else:seqs[key]+=l.strip()
assert len(seqs)==13 and {len(s) for s in seqs.values()}=={53}
columns=[i for i in range(53) if all(seqs[k][i] not in '-.' for k in seqs)];assert len(columns)==52
parser=PDBParser(QUIET=True);xyz={};length={}
for k,s in seqs.items():
 residues=[r for r in parser.get_structure(k,I/'seed_fragments'/(k+'.pdb')).get_residues() if 'CA' in r]
 pdbseq=''.join(seq1(r.resname) for r in residues);assert pdbseq==s.replace('-','').replace('.',''),k
 col_to_res={};j=0
 for i,a in enumerate(s):
  if a not in '-.':col_to_res[i]=j;j+=1
 xyz[k]=np.array([residues[col_to_res[i]]['CA'].coord.astype(float) for i in columns]);length[k]=len(residues)
old={tuple(sorted([r['seed_a'],r['seed_b']])):float(r['pairwise_fitted_rmsd_angstrom']) for r in csv.DictReader((I/'archived_pairwise_rmsd.csv').open())}
rows=[]
for a,b in itertools.combinations(sorted(seqs),2):
 x,y=xyz[a],xyz[b];xc=x-x.mean(0);yc=y-y.mean(0);u,sv,vt=np.linalg.svd(xc.T@yc);rot=u@vt
 if np.linalg.det(rot)<0:u[:,-1]*=-1;rot=u@vt
 rms=float(np.sqrt(np.mean(np.sum((xc@rot-yc)**2,axis=1))));q=52**2/((1+(rms/3)**2)*length[a]*length[b])
 assert abs(rms-old[(a,b)])<1e-5,(a,b,rms,old[(a,b)])
 rows.append(dict(seed_a=a,seed_b=b,aligned_CA_pairs=52,N1=length[a],N2=length[b],rmsd_angstrom=rms,Q_score=q,archived_rmsd_angstrom=old[(a,b)]))
assert len(rows)==78
with (O/'training_pairwise_quality.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
meta={r['seed_id']:r for r in csv.DictReader((I/'fragment_manifest.csv').open())};summary=[]
for a in sorted(seqs):
 r=[r for r in rows if a in [r['seed_a'],r['seed_b']]];m=meta[a]
 summary.append(dict(seed_id=a,pdb=m['entity_identifier'].split('_')[0],chain=m['source_label_chain'],method=m['experimental_method'],length=length[a],mean_pairwise_rmsd_angstrom=float(np.mean([x['rmsd_angstrom'] for x in r])),maximum_pairwise_rmsd_angstrom=max(x['rmsd_angstrom'] for x in r),mean_pairwise_Q=float(np.mean([x['Q_score'] for x in r]))))
with (O/'training_quality_summary.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
result=dict(training_fragments=13,alignment_columns=53,common_gap_free_columns=52,pairs_recomputed=78,maximum_archive_difference=max(abs(r['rmsd_angstrom']-r['archived_rmsd_angstrom']) for r in rows),mean_pairwise_rmsd_angstrom=float(np.mean([r['rmsd_angstrom'] for r in rows])),minimum_pairwise_rmsd_angstrom=min(r['rmsd_angstrom'] for r in rows),maximum_pairwise_rmsd_angstrom=max(r['rmsd_angstrom'] for r in rows),mean_pairwise_Q=float(np.mean([r['Q_score'] for r in rows])),minimum_pairwise_Q=min(r['Q_score'] for r in rows),maximum_pairwise_Q=max(r['Q_score'] for r in rows),Q_definition='Nalign^2 / ((1+(RMSD/3A)^2)*N1*N2), using fixed MUSTANG correspondence; not a new PDBeFold alignment',source_coordinates_modified=False,training_membership_modified=False)
(O/'training_quality_results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
