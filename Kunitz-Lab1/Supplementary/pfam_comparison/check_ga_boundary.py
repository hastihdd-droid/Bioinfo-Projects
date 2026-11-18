from pathlib import Path
import sys,json,csv
sys.path.insert(0,str(Path('work/pfam_packages').resolve()))
import pyhmmer
P=Path('work/pfam_comparison');targets={}
for l in (P/'pfam.tbl').read_text().splitlines():
 if not l or l.startswith('#'):continue
 a=l.split();seq=float(a[5]);dom=float(a[8])
 if seq==21.0 or dom==21.0:targets[a[0]]=(seq,dom)
seqs=[];head=None;parts=[];alpha=pyhmmer.easel.Alphabet.amino()
def save():
 if head in targets:seqs.append(pyhmmer.easel.TextSequence(name=head,sequence=''.join(parts)).digitize(alpha))
with (P/'evaluation.fasta').open() as f:
 for line in f:
  if line.startswith('>'):save();head=line[1:].strip();parts=[]
  else:parts.append(line.strip())
 save()
with pyhmmer.plan7.HMMFile(P/'PF00014.hmm') as f:model=f.read()
rows=[]
if seqs:
 hits=next(pyhmmer.hmmer.hmmsearch(model,pyhmmer.easel.DigitalSequenceBlock(alpha,seqs),cpus=2,bias_filter=False,F1=1.,F2=1.,F3=1.,E=1000000,domE=1000000,seed=42))
 for h in hits:rows.append(dict(Entry=h.name,sequence_score_native=h.score,domain_score_native=h.best_domain.score,GA_positive=h.score>=21.0 and h.best_domain.score>=21.0))
(P/'GA_boundary_native_scores.json').write_text(json.dumps({'scope':'Native-score check for hits rounded to the published21.0bit boundary; primary comparison remains one-decimal','targets':len(targets),'rows':rows},indent=2));print(json.dumps(rows,indent=2))
