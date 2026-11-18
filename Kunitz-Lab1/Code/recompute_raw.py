import csv,json,math,hashlib
from pathlib import Path
root=Path('work/raw_audit')
out=Path('work/raw_recomputed');out.mkdir(exist_ok=True)
def read(rel):
 with (root/rel).open() as f:return list(csv.DictReader(f))
def j(rel):return json.loads((root/rel).read_text())
def rel(p):return p.removeprefix('/content/kunitz_project/')
def metrics(rows,t):
 tp=tn=fp=fn=0
 for row in rows:
  y=int(row['reference_label']);s=row['raw_best_domain_score'];pred=s is not None and s>=t
  if y and pred:tp+=1
  elif y:fn+=1
  elif pred:fp+=1
  else:tn+=1
 den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
 return dict(TP=tp,TN=tn,FP=fp,FN=fn,sensitivity=tp/(tp+fn) if tp+fn else None,specificity=tn/(tn+fp) if tn+fp else None,precision=tp/(tp+fp) if tp+fp else None,MCC=(tp*tn-fp*fn)/den if den else None)
allrows={};audit={}; hashes={}
for split,stamp in [('validation','20261107T212941_588852Z'),('test','20261107T212942_643707Z')]:
 dr=f'results/{split}/{stamp}'
 meta=read(f'data/processed/evaluation_split/{split}_metadata.csv')
 saved=read(f'{dr}/{split}_scores.csv')
 sd={r['Entry']:r for r in saved}
 hits={}
 for line in (root/f'{dr}/{split}_hits.tbl').read_text().splitlines():
  if not line or line.startswith('#'):continue
  f=line.split(); assert f[0] not in hits;hits[f[0]]=float(f[8])
 assert set(hits)<=set(sd)
 assert len(meta)==len(saved)==1914
 assert {r['Entry'] for r in meta}==set(sd)
 for row in meta:
  entry=row['Entry'];score=hits.get(entry);row['raw_best_domain_score']=score
  assert int(row['reference_label'])==int(sd[entry]['reference_label'])
  assert row['Sequence']==sd[entry].get('Sequence',row['Sequence'])
  sv=sd[entry]['best_domain_score_bits'];assert (score is None and not sv) or (score is not None and float(sv)==score)
 allrows[split]=meta
 audit[split]={'rows':len(meta),'raw_hits':len(hits),'missing_positive_hits':sum(int(r['reference_label'])==1 and r['raw_best_domain_score'] is None for r in meta)}
 record=j(f'{dr}/run_record.json');hashes[split]=record['model_sha256']
# Independently optimize over finite raw validation scores using midpoints.
v=allrows['validation'];unique=sorted({r['raw_best_domain_score'] for r in v if r['raw_best_domain_score'] is not None})
sweep=[]
for lo,hi in zip(unique[:-1],unique[1:]):
 m=metrics(v,(lo+hi)/2);sweep.append((m['MCC'],lo,hi))
best=max(x[0] for x in sweep);opt=[x for x in sweep if x[0]==best]
lo,hi=max([(x[1],x[2]) for x in opt],key=lambda x:(x[1]-x[0],(x[0]+x[1])/2))
t=(lo+hi)/2
assert math.isclose(t,8.7) and hashes['test']==hashes['validation']
for split,rows in allrows.items():
 audit[split]['recomputed_metrics']=metrics(rows,t)
 scores={label:[r['raw_best_domain_score'] for r in rows if int(r['reference_label'])==label and r['raw_best_domain_score'] is not None] for label in [0,1]}
 audit[split]['lowest_positive']=min(scores[1]);audit[split]['highest_reported_negative']=max(scores[0])
# Independently check frozen membership and overlap from complete sequences.
v=allrows['validation'];trows=allrows['test']
for column in ['Entry','Sequence','similarity_group']:
 assert not ({r[column] for r in v}&{r[column] for r in trows}),column
# Verify exclusions and membership of named difficult examples.
positives=read('data/processed/SwissProt_positive_candidates.csv')
eligible=read('data/processed/SwissProt_eligible_positives.csv')
assert len(positives)==368 and len(eligible)==348
assert {r['Entry'] for r in eligible}=={r['Entry'] for r in v+trows if int(r['reference_label'])==1}
training_audit=j('results/robustness/current_training_audit.json')
training=read(rel(training_audit['audit_path']));print('training fields',list(training[0]))
# Recompute expanded-background check from actual membership tables.
neg=read('data/processed/SwissProt_negative_candidate_metadata.csv')
negids={r['Entry'] for r in neg};assert len(neg)==len(negids)==566072
scan='results/swissprot_scan/20261107T212956_904384Z'
preds=read(scan+'/classification/SwissProt_predicted_Kunitz.csv')
meta=read(scan+'/distribution/metadata_corrected_20261107T213427_876118Z/predicted_protein_distribution_metadata.csv')
md={r['Entry']:r for r in meta}; assert len(md)==426
old_neg={r['Entry'] for r in v+trows if int(r['reference_label'])==0}
expanded=negids-old_neg
fpids={r['Entry'] for r in preds if r['Entry'] in expanded}
assert len(expanded)==562592 and len(fpids)==15
# Reconstruct independently, without using stored predicted labels.
scanraw={}
with (root/scan/'swissprot_hits.tbl').open() as f:
 for line in f:
  if line.startswith('#') or not line.strip():continue
  x=line.split();accession=x[0].split('|')[1] if x[0].startswith(('sp|','tr|')) else x[0]
  assert accession not in scanraw
  scanraw[accession]=float(x[8])
raw_fpids={x for x in expanded if scanraw.get(x,-math.inf)>=8.7}
assert raw_fpids==fpids
rawexpanded=[dict(reference_label='0',raw_best_domain_score=scanraw.get(x)) for x in expanded]
expandedmetric=metrics([r for r in trows if int(r['reference_label'])==1]+rawexpanded,8.7)
results={'status':'Recomputed from verified extracted raw HMMER tblout and metadata; HMMER itself was not rerun','threshold_recomputed_from_validation':t,'optimal_interval':{'lower_exclusive':lo,'upper_inclusive':hi},'same_model_hash_validation_test':hashes['test'],'no_shared_accessions_sequences_or_80pct_groups':True,'splits':audit,'expanded_background_audit':expandedmetric,'expanded_negative_count':len(expanded),'expanded_FP_accessions':sorted(fpids),'positive_exclusions_total':20,'positive_exclusion_criterion':'training overlap or close training similarity, not HMM scores'}
(out/'raw_recalculation.json').write_text(json.dumps(results,indent=2))
with (out/'test_positive_raw_scores.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['Entry','similarity_group','raw_best_domain_score']);w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in trows if int(r['reference_label'])==1)
print(json.dumps(results,indent=2))
