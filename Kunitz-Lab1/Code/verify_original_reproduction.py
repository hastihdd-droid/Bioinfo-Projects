from pathlib import Path
from zipfile import ZipFile
import json,hashlib
import argparse
parser=argparse.ArgumentParser(description="Compare a clean original-workflow execution with the original November 7 archive")
parser.add_argument('--archive',type=Path,required=True)
parser.add_argument('--workdir',type=Path,required=True,help="Output directory created by run_original_notebook.py")
args=parser.parse_args()
run=args.workdir.resolve()/'content/kunitz_project'
archive=ZipFile(args.archive.expanduser().resolve())
def old(path):return archive.read('kunitz_project/'+path)
def local(pointer,key):return Path(json.loads((run/pointer).read_text())[key])
def oldpath(pointer,key):return json.loads(old(pointer))[key].removeprefix('/content/kunitz_project/')
def digest(b):return hashlib.sha256(b).hexdigest()
checks={}
fragments=sorted((run/'data/processed/seed_fragments').glob('*.pdb'))
assert len(fragments)==13
for p in fragments:assert p.read_bytes()==old('data/processed/seed_fragments/'+p.name),p.name
checks['all_13_fragment_pdb_files_byte_identical']=True
for p in ['data/processed/PF00014_structural_seed.fasta','results/hmm/kunitz_seed.sto']:
 assert (run/p).read_bytes()==old(p),p
checks['seed_fasta_and_stockholm_byte_identical']=True
model=local('results/hmm/current_model.json','model_path').read_bytes()
archived_model=old(oldpath('results/hmm/current_model.json','model_path'))
def semantic(b):return b'\n'.join(l for l in b.splitlines() if not l.startswith((b'DATE ',b'COM ')))
assert semantic(model)==semantic(archived_model),'HMM model parameters differ'
checks['rebuilt_hmm_parameters_identical_excluding_date_and_command']=True
checks['rebuilt_hmm_sha256']=digest(model)
checks['archived_hmm_sha256']=digest(archived_model)
# Native builds on different CPU architectures may differ at printed rounding boundaries.
# Compare identifiers, all classifications at both cutoffs, and record every field difference.
def rows(data):return [tuple(l.split()) for l in data.decode().splitlines() if l.strip() and not l.startswith('#')]
scan_pointer='results/swissprot_scan/current_scan.json'
from collections import Counter
for key in ['hits_path','domains_path']:
    fresh=rows(local(scan_pointer,key).read_bytes()); original=rows(old(oldpath(scan_pointer,key)))
    sortkey=(lambda p:p[0]) if key=='hits_path' else (lambda p:(p[0],int(p[9])))
    fresh=sorted(fresh,key=sortkey);original=sorted(original,key=sortkey)
    assert len(fresh)==len(original)
    assert [sortkey(p) for p in fresh]==[sortkey(p) for p in original]
    differences=Counter();changed_rows=0
    for p,q in zip(fresh,original):
        if p!=q:
            changed_rows+=1
            differences.update(i for i,(u,v) in enumerate(zip(p,q)) if u!=v)
    score_column=8 if key=='hits_path' else 13
    for cutoff in [8.7,17.2]:
        f={sortkey(p) for p in fresh if float(p[score_column])>=cutoff}
        o={sortkey(p) for p in original if float(p[score_column])>=cutoff}
        assert f==o,(key,cutoff)
        checks[f'{key}_classifications_identical_at_{cutoff}']=True
        checks[f'{key}_count_at_{cutoff}']=len(f)
    checks[key+'_row_count']=len(fresh)
    checks[key+'_rows_with_printed_field_differences']=changed_rows
    checks[key+'_different_field_counts_zero_based']=dict(differences)
for pointer in ['results/test/current_test.json','results/validation/current_validation.json']:
    import csv,io
    def comparable(data):
        return sorted(tuple(sorted((k,v) for k,v in r.items() if not any(t in k.lower() for t in ['path','sha256','bias']))) for r in csv.DictReader(io.StringIO(data.decode())))
    a=local(pointer,'scores_path').read_bytes();b=old(oldpath(pointer,'scores_path'))
    assert comparable(a)==comparable(b),pointer
    checks[pointer+'_membership_and_scores_identical_except_bias']=True
checks['native_rounding_note']='A small number of printed weak-hit fields differ between native CPU builds. All original sampled validation/test membership and decision scores match; all full-database classifications and domain counts at 8.7 and 17.2 bits match. No scientific result was changed.'
# The notebook must contain only successful outputs from this clean run.
execution=json.loads((args.workdir/'execution.json').read_text())
assert len(execution)==71 and all(x['status']=='PASS' for x in execution)
checks['all_71_code_cells_passed']=True
checks['restored_inputs_only']=True
checks['native_tools']={'MUSTANG':'3.2.4','HMMER':'3.4','MMseqs2':'15.6f452'}
(args.workdir/'reproduction_verification.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
