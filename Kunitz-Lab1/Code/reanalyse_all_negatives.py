import csv,json,math,hashlib
from collections import defaultdict,Counter
from pathlib import Path
import sys
sys.path.insert(0,str(Path('work/report_packages').resolve()))
import numpy as np
R=Path('work/raw_audit');O=Path('work/reanalysis');O.mkdir(exist_ok=True)
design=json.loads(Path('work/expanded_design.json').read_text())
(O/'design.json').write_text(json.dumps(design,indent=2))
def read(p):
 with (R/p).open() as f:return list(csv.DictReader(f))
assignment=read('data/processed/evaluation_split/split_assignments.csv');original={r['Entry']:r for r in assignment}
neg=read('data/processed/SwissProt_negative_candidate_metadata.csv');pos=read('data/processed/SwissProt_eligible_positives.csv')
labels={r['Entry']:0 for r in neg};labels.update({r['Entry']:1 for r in pos});assert len(labels)==566420
parent={x:x for x in labels}
def find(x):
 while x!=parent[x]:parent[x]=parent[parent[x]];x=parent[x]
 return x
def union(a,b):
 a,b=find(a),find(b)
 if a!=b:parent[max(a,b)]=min(a,b)
oldgroups={}
for r in assignment:
 g=r['similarity_group'];x=r['Entry']
 if g in oldgroups:union(x,oldgroups[g])
 else:oldgroups[g]=x
seq_first={};seen=set();seqhash={}
fasta=R/'data/raw/uniprot/swissprot_complete_2026_03/swissprot_complete.fasta'
def consume(header,lines):
 if not header:return
 x=header.split()[0].split('|')[1];
 if x not in labels:return
 sequence=''.join(lines);h=hashlib.sha256(sequence.encode()).hexdigest();seqhash[x]=h;seen.add(x)
 if h in seq_first:union(x,seq_first[h])
 else:seq_first[h]=x
header=None;lines=[]
with fasta.open() as f:
 for line in f:
  if line.startswith('>'):consume(header,lines);header=line[1:].strip();lines=[]
  else:lines.append(line.strip())
consume(header,lines);assert seen==set(labels)
groups=defaultdict(list)
for x in labels:groups[find(x)].append(x)
alloc={};conflicts=[]
for root,members in groups.items():
 ys={labels[x] for x in members}
 if len(ys)>1:conflicts.extend(members);continue
 old={original[x]['split'] for x in members if x in original}
 assert len(old)<=1
 side=next(iter(old)) if old else ('validation' if int(hashlib.sha256(('42:'+min(members)).encode()).hexdigest(),16)%2==0 else 'test')
 for x in members:alloc[x]=side
# Persist the split before reading scores.
with (O/'assignments.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['Entry','reference_label','split','group','sequence_sha256'])
 for x in sorted(alloc):w.writerow([x,labels[x],alloc[x],find(x),seqhash[x]])
print('Frozen new split:',dict(Counter((alloc[x],labels[x]) for x in alloc)),'conflicts',len(conflicts),flush=True)
# Read only fixed scoring outputs; retain all no-hit proteins.
scores={}
scan='results/swissprot_scan/20261107T212956_904384Z/swissprot_hits.tbl'
with (R/scan).open() as f:
 for line in f:
  if line.startswith('#') or not line.strip():continue
  z=line.split();x=z[0].split('|')[1] if z[0].startswith(('sp|','tr|')) else z[0];scores[x]=float(z[8])
diffs=[]
for split,stamp in [('validation','20261107T212941_588852Z'),('test','20261107T212942_643707Z')]:
 old=read(f'results/{split}/{stamp}/{split}_scores.csv')
 for r in old:
  x=r['Entry'];a=float(r['best_domain_score_bits']) if r['best_domain_score_bits'] else -math.inf;b=scores.get(x,-math.inf)
  if a!=b:diffs.append({'Entry':x,'label':int(r['reference_label']),'old':a,'full':b})
print('Score discrepancies',len(diffs),diffs[:10],flush=True)
# Database E-value reporting may add extremely weak hits. All positive scores
# must agree; any changed negative must remain below the original cutoff.
assert all(d['label']==0 and max(d['old'],d['full'])<8.7 for d in diffs)
ids={side:sorted(x for x in alloc if alloc[x]==side) for side in ['validation','test']}
def metric(side,threshold):
 count=Counter()
 for x in ids[side]:
  p=scores.get(x,-math.inf)>=threshold;y=labels[x]
  count['TP' if y and p else 'FN' if y else 'FP' if p else 'TN']+=1
 tp,fn,fp,tn=[count[x] for x in ['TP','FN','FP','TN']]
 den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
 return dict(TP=tp,FN=fn,FP=fp,TN=tn,sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),precision=tp/(tp+fp) if tp+fp else 0,MCC=(tp*tn-fp*fn)/den if den else 0,F1=2*tp/(2*tp+fp+fn))
# Efficient cumulative validation sweep over unique score bins.
bins=defaultdict(lambda:[0,0]);total=Counter(labels[x] for x in ids['validation'])
for x in ids['validation']:
 s=scores.get(x,-math.inf)
 if math.isfinite(s):bins[s][labels[x]]+=1
values=sorted(bins);tp=total[1];fp=total[0]
# Remove missing hits, then move each finite bin below cutoff.
missing=Counter(labels[x] for x in ids['validation'] if x not in scores);tp-=missing[1];fp-=missing[0]
sweep=[]
for lo,hi in zip(values[:-1],values[1:]):
 fp-=bins[lo][0];tp-=bins[lo][1];fn=total[1]-tp;tn=total[0]-fp
 den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));m=(tp*tn-fp*fn)/den if den else 0
 sweep.append(dict(lower=lo,upper=hi,MCC=m,TP=tp,FP=fp,FN=fn,TN=tn))
best=max(r['MCC'] for r in sweep);optimal=[r for r in sweep if abs(r['MCC']-best)<1e-12]
# Merge adjoining tied intervals, as specified before scoring.
blocks=[]
for r in optimal:
 if blocks and blocks[-1][1]==r['lower']:blocks[-1][1]=r['upper']
 else:blocks.append([r['lower'],r['upper']])
lo,hi=max(blocks,key=lambda v:(v[1]-v[0],(v[0]+v[1])/2));threshold=(lo+hi)/2
selected={'threshold_bits':threshold,'optimal_interval_lower_exclusive':lo,'optimal_interval_upper_inclusive':hi,'metric':'MCC','design_sha256':hashlib.sha256(Path('work/expanded_design.json').read_bytes()).hexdigest(),'split_sha256':hashlib.sha256((O/'assignments.csv').read_bytes()).hexdigest(),'validation_metrics':metric('validation',threshold),'test_metrics_not_used':True}
(O/'selected_threshold.json').write_text(json.dumps(selected,indent=2))
# Evaluate test only after the threshold record has been written.
result=dict(status='Retrospective expanded-negative reanalysis of archived scores; fixed HMM',selected_threshold=selected,test_metrics=metric('test',threshold),counts={side:{str(y):sum(labels[x]==y for x in ids[side]) for y in [0,1]} for side in ids},label_conflicts_excluded=conflicts,weak_negative_score_discrepancies=diffs,original_cutoff_metrics={side:metric(side,8.7) for side in ids})
with (O/'test_predictions.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['Entry','reference_label','score_bits','predicted_label','outcome'])
 for x in ids['test']:
  s=scores.get(x,-math.inf);p=s>=threshold;y=labels[x]
  w.writerow([x,y,s if math.isfinite(s) else '',int(p),'TP' if y and p else 'FN' if y else 'FP' if p else 'TN'])
with (O/'threshold_sweep.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(sweep[0]));w.writeheader();w.writerows(sweep)
(O/'results.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='weak_negative_score_discrepancies'},indent=2))
