from pathlib import Path
import csv,json,re
R=Path('work/raw_audit');O=Path('outputs')
with (R/'data/processed/evaluation_split/test_metadata.csv').open() as f:pass
with open('work/reanalysis/assignments.csv') as f:assign={r['Entry']:r for r in csv.DictReader(f)}
with (R/'data/processed/PF00014_fragment_manifest.csv').open() as f:manifest=list(csv.DictReader(f))
train_ids={x.strip() for r in manifest for x in r['uniprot_ids'].split(';') if x.strip()}
def fasta(path):
 header=None;seq=[]
 with path.open() as f:
  for line in f:
   if line.startswith('>'):
    if header is not None:yield header,''.join(seq)
    header=line[1:].strip();seq=[]
   else:seq.append(line.strip())
 if header is not None:yield header,''.join(seq)
seeds=dict(fasta(R/'data/processed/PF00014_structural_seed.fasta'))
patterns=[(k,re.compile(re.escape(s))) for k,s in seeds.items()]
exact=[];seen=0
for header,seq in fasta(R/'data/raw/uniprot/swissprot_complete_2026_03/swissprot_complete.fasta'):
 entry=header.split()[0].split('|')[1]
 if entry not in assign:continue
 seen+=1
 matches=[k for k,p in patterns if p.search(seq)]
 if matches:exact.append(dict(Entry=entry,split=assign[entry]['split'],reference_label=assign[entry]['reference_label'],seed_ids=matches))
assert seen==len(assign)==566420
r={'scope':'All 566420 proteins in retrospective expanded evaluation, including all added negatives','training_accessions':sorted(train_ids),'training_accession_overlap':sorted(train_ids&set(assign)),'proteins_containing_exact_complete_training_fragment':exact,'training_fragments':len(seeds),'evaluated_proteins':seen,'interpretation':'Exact full seed-segment and accession checks only; does not rule out shorter matches, homologous relationships, shared ancestry or retrospective test reuse.'}
(O/'Kunitz_Direct_Training_Overlap_Audit.json').write_text(json.dumps(r,indent=2))
print(json.dumps(r,indent=2))
