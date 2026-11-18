from pathlib import Path
import json,math
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
BASE=Path(__file__).resolve().parent;I=BASE/'inputs';O=BASE/'results'
records=json.loads((O/'chimerax_results.json').read_text());cases=json.loads((I/'case_reference.json').read_text());checks=[]
for r in records:
 assert r['status']=='success'
 case=cases[r['Entry']];metadata=json.loads((I/(r['Entry']+'_metadata.json')).read_text())[0]
 assert metadata['sequence']==case['original_sequence']
 residues=[res for res in PDBParser(QUIET=True).get_structure(r['Entry'],I/(r['Entry']+'.pdb')).get_residues() if 'CA' in res]
 assert ''.join(seq1(x.resname) for x in residues)==case['original_sequence']
 assert case['reference_label']==1 and case['score_bits']<17.2 and case['outcome']=='FN'
 ref={x['residue']:np.array(x['coord']) for x in r['reference_trace']};tar={x['residue']:np.array(x['coord']) for x in r['target_trace']}
 squared=[float(np.sum((ref[x['reference_residue']]-tar[x['target_residue']])**2)) for x in r['mapping']];core=[x for x,m in zip(squared,r['mapping']) if m['is_pruned_core']]
 full=math.sqrt(sum(squared)/len(squared));pruned=math.sqrt(sum(core)/len(core));assert abs(full-r['full_rmsd_angstrom'])<1e-4 and abs(pruned-r['core_rmsd_angstrom'])<1e-4
 assert max(core)<=4.0001
 pl=[m['target_pLDDT'] for m in r['mapping'] if m['is_pruned_core']];assert abs(sum(pl)/len(pl)-r['mean_core_pLDDT'])<1e-4
 checks.append(dict(Entry=r['Entry'],split=case['split'],sequence_exactly_matches_project=True,full_pairs=len(squared),pruned_core_pairs=len(core),full_rmsd_independently_verified=full,core_rmsd_independently_verified=pruned,mean_core_pLDDT=r['mean_core_pLDDT']))
assert len(checks)==3
(O/'structure_verification.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
