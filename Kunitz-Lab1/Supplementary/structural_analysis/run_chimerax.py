# Run with ChimeraX --nogui --exit --script /absolute/path/run_chimerax.py
from pathlib import Path
import json,traceback
from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
BASE=Path(__file__).resolve().parent;I=BASE/'inputs';O=BASE/'results';O.mkdir(exist_ok=True)
records=[]
for acc in ['D3GGZ8','O62247','A0A1Q1NL17']:
 try:
  run(session,'close all')
  run(session,'open '+str(I/'3TGI_chain_I.pdb'))
  run(session,'open '+str(I/(acc+'.pdb')))
  cmd='matchmaker #2/A to #1/I pairing ss matrix BLOSUM-62 alg nw ssFraction 0.3 gapExtend 1 hgap 18 sgap 18 ogap 6 computeSS true cutoffDistance 2.0 showAlignment false logParameters true'
  results=run(session,cmd)
  assert len(results)==1
  result=results[0]
  full_ref=result['full ref atoms'];full_match=result['full match atoms'];core_ref=result['final ref atoms'];core_match=result['final match atoms']
  mapping=[dict(reference_residue=int(r.residue.number),target_residue=int(t.residue.number),target_pLDDT=float(t.bfactor),is_pruned_core=any(t==x for x in core_match)) for r,t in zip(full_ref,full_match)]
  target=next(m for m in session.models.list() if isinstance(m,AtomicStructure) and m.id==(2,))
  ref=next(m for m in session.models.list() if isinstance(m,AtomicStructure) and m.id==(1,))
  # Scene coordinates reflect the actual Matchmaker transform.
  def ca_trace(model):
   return [dict(residue=int(a.residue.number),resname=a.residue.name,coord=[float(x) for x in a.scene_coord],bfactor=float(a.bfactor)) for a in model.atoms if a.name=='CA']
  record=dict(Entry=acc,status='success',command=cmd,full_pairs=len(full_ref),core_pairs=len(core_ref),full_rmsd_angstrom=float(result['full RMSD']),core_rmsd_angstrom=float(result['final RMSD']),mean_core_pLDDT=sum(float(a.bfactor) for a in core_match)/len(core_match),minimum_core_pLDDT=min(float(a.bfactor) for a in core_match),mapping=mapping,reference_trace=ca_trace(ref),target_trace=ca_trace(target),transform=result['transformation matrix'].matrix.tolist())
  records.append(record)
  run(session,'save '+str(O/(acc+'_superposed.pdb'))+' models #1,2')
  print('STRUCTURAL_RESULT',acc,record['full_pairs'],record['core_pairs'],record['full_rmsd_angstrom'],record['core_rmsd_angstrom'],record['mean_core_pLDDT'],flush=True)
 except Exception as e:records.append(dict(Entry=acc,status='error',error=str(e),traceback=traceback.format_exc()));print(traceback.format_exc(),flush=True)
(O/'chimerax_results.json').write_text(json.dumps(records,indent=2))
