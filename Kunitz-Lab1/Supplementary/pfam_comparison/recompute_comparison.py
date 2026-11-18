from pathlib import Path
import csv,json,math,hashlib,collections,gzip
BASE=Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
P=BASE/'computed';P.mkdir(exist_ok=True)
I=BASE/'inputs'
protocol=json.loads((BASE/'protocol.json').read_text());assert protocol['structural_model_rule'].startswith('Original fixed')
with gzip.open(I/'benchmark_membership.csv.gz','rt') as f:assignment={r['Entry']:r for r in csv.DictReader(f)}
def read_tbl(p):
 hits={}
 with gzip.open(p,'rt') as f:
  for l in f:
   if l.startswith('#') or not l.strip():continue
   v=l.split();assert v[0] not in hits;hits[v[0]]={'domain':float(v[8]),'sequence':float(v[5])}
 assert set(hits)<=set(assignment)
 return hits
def metrics(rows,hits,threshold,seqthreshold=None):
 c=collections.Counter()
 for a in rows:
  r=assignment[a];h=hits.get(a);y=int(r['reference_label']);pr=h is not None and h['domain']>=threshold and (seqthreshold is None or h['sequence']>=seqthreshold)
  c['TP' if y and pr else 'FN' if y else 'FP' if pr else 'TN']+=1
 tp,fn,fp,tn=[c[x] for x in ['TP','FN','FP','TN']];den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
 return dict(TP=tp,FN=fn,FP=fp,TN=tn,sensitivity=tp/(tp+fn),precision=tp/(tp+fp) if tp+fp else None,specificity=tn/(tn+fp),F1=2*tp/(2*tp+fp+fn),MCC=(tp*tn-fp*fn)/den if den else 0)
v={a for a,r in assignment.items() if r['split']=='validation'};test={a for a in assignment}-v
pf=read_tbl(I/'pfam.tbl.gz');st=read_tbl(I/'structural.tbl.gz')
# Calibrate reference using validation membership and scores only.
counts=collections.defaultdict(collections.Counter)
for a in v:
 if a in pf:counts[pf[a]['domain']][int(assignment[a]['reference_label'])]+=1
unique=sorted(counts);npositive=sum(int(assignment[a]['reference_label']) for a in v);nnegative=len(v)-npositive
tp=sum(x[1] for x in counts.values());fp=sum(x[0] for x in counts.values());fn=npositive-tp;tn=nnegative-fp;sweep=[]
for lo,hi in zip(unique[:-1],unique[1:]):
 tp-=counts[lo][1];fn+=counts[lo][1];fp-=counts[lo][0];tn+=counts[lo][0]
 den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));m=(tp*tn-fp*fn)/den if den else 0
 sweep.append(dict(lower_exclusive=lo,upper_inclusive=hi,midpoint=(lo+hi)/2,MCC=m,TP=tp,FN=fn,FP=fp,TN=tn))
best=max(x['MCC'] for x in sweep);optimal=[x for x in sweep if x['MCC']==best];chosen=max(optimal,key=lambda x:(x['upper_inclusive']-x['lower_exclusive'],x['midpoint']));cut=chosen['midpoint']
selected=dict(model='PF00014.31',threshold_bits=cut,validation_optimal_interval=[chosen['lower_exclusive'],chosen['upper_inclusive']],criterion='validation MCC',protocol_sha256=hashlib.sha256((BASE/'protocol.json').read_bytes()).hexdigest(),test_metrics_used=False)
(P/'selected_reference_threshold.json').write_text(json.dumps(selected,indent=2))
with (P/'reference_validation_sweep.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(sweep[0]));w.writeheader();w.writerows(sweep)
# Only after persisting selected cutoff, evaluate the fixed test membership.
record=json.loads((BASE/'search_record.json').read_text());gaseq,gadom=record['gathering_cutoffs']
# Published GA applies to native scores. Recompute boundary candidates to avoid rounding a subthreshold score upward.
pfGA={a:dict(h) for a,h in pf.items()}
boundary=json.loads((BASE/'GA_boundary_native_scores.json').read_text())
for b in boundary['rows']:pfGA[b['Entry']]={'domain':b['domain_score_native'],'sequence':b['sequence_score_native']}
models=[('Structural HMM',st,17.2,None),('Pfam validation-calibrated',pf,cut,None),('Pfam curated GA',pfGA,gadom,gaseq)]
results=[]
for name,h,t,seqt in models:results.append(dict(model=name,domain_cutoff_bits=t,sequence_cutoff_bits=seqt,validation=metrics(v,h,t,seqt),test=metrics(test,h,t,seqt)))
expected=json.loads((I/'original_results.json').read_text())
assert results[0]['test']==expected['test_metrics']
for k in ['TP','FN','FP','TN','MCC']:assert results[0]['validation'][k]==expected['selected_threshold']['validation_metrics'][k]
old={}
with gzip.open(I/'original_test_predictions.csv.gz','rt') as f:
 for r in csv.DictReader(f):old[r['Entry']]=r
assert all(int(r['predicted_label'])==int(a in st and st[a]['domain']>=17.2) for a,r in old.items())
disagreements=[];predictions=[]
for a in sorted(test):
 y=int(assignment[a]['reference_label']);sr=st.get(a);pr=pf.get(a);sp=bool(sr and sr['domain']>=17.2);pp=bool(pr and pr['domain']>=cut);g=pfGA.get(a);gp=bool(g and g['domain']>=gadom and g['sequence']>=gaseq)
 row=dict(Entry=a,reference_label=y,structural_score_bits=sr['domain'] if sr else '',pfam_domain_score_bits=pr['domain'] if pr else '',pfam_sequence_score_bits=pr['sequence'] if pr else '',structural_prediction=int(sp),pfam_calibrated_prediction=int(pp),pfam_GA_prediction=int(gp))
 predictions.append(row)
 if sp!=pp or sp!=gp:disagreements.append(row)
for name,rows in [('test_predictions.csv',predictions),('test_disagreements.csv',disagreements)]:
 with (P/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(predictions[0]));w.writeheader();w.writerows(rows)
summary={'scope':'Reference-model comparison on the existing expanded benchmark; not external validation','source_model':record,'selected_reference_threshold':selected,'models':results,'structural_test_decisions_reproduced':len(old),'test_disagreements':disagreements,'limits':protocol['limits'],'GA_boundary_native_check':boundary}
(P/'comparison_results.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
