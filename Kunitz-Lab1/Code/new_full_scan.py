import csv,json,math
from pathlib import Path
R=Path('work/raw_audit');O=Path('work/reanalysis')
result=json.loads((O/'results.json').read_text());t=result['selected_threshold']['threshold_bits']
scan=R/'results/swissprot_scan/20261107T212956_904384Z'
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
oldpred=read(scan/'classification/SwissProt_predicted_Kunitz.csv');oldpos={r['Entry'] for r in read(scan/'annotation_comparison/predicted_with_PF00014.csv')}
meta=read(scan/'distribution/summary/predictions_with_taxonomic_bins.csv');md={r['Entry']:r for r in meta}
scores={}
with (scan/'swissprot_hits.tbl').open() as f:
 for line in f:
  if not line.strip() or line.startswith('#'):continue
  x=line.split();entry=x[0].split('|')[1] if x[0].startswith(('sp|','tr|')) else x[0];scores[entry]=float(x[8])
newids={x for x,s in scores.items() if s>=t};assert newids<={r['Entry'] for r in oldpred}
domcount={x:0 for x in newids};domains=[];oldcount=0
with (scan/'swissprot_domains.domtbl').open() as f:
 for line in f:
  if not line.strip() or line.startswith('#'):continue
  x=line.split();entry=x[0].split('|')[1] if x[0].startswith(('sp|','tr|')) else x[0];s=float(x[13]);oldcount+=s>=8.7
  if s>=t:
   assert entry in newids;domcount[entry]+=1;domains.append(dict(Entry=entry,score=s,hmm_start=int(x[15]),hmm_end=int(x[16]),ali_start=int(x[17]),ali_end=int(x[18])))
assert oldcount==514 and all(n>0 for n in domcount.values())
miss=oldpos-newids;extra=newids-oldpos
review=read(scan/'discordant_domain_review/interpretation/discordant_prediction_interpretation.csv');rv={r['Entry']:r for r in review}
from collections import Counter
stats={'threshold_bits':t,'database_proteins':575748,'predicted_proteins':len(newids),'PF00014_detected':len(oldpos&newids),'PF00014_missed':len(miss),'additional_predictions':len(extra),'annotation_agreement':len(oldpos&newids)/len(newids),'domains_above_threshold':len(domains),'multidomain_proteins':sum(n>1 for n in domcount.values()),'review_counts':dict(Counter(rv[x]['review_category'] for x in extra)),'taxonomic_counts':dict(Counter(md[x]['taxonomic_bin'] for x in newids)),'explicit_fragment_count':sum(md[x]['is_explicit_fragment'].lower()=='true' for x in newids),'missed_PF00014':[dict(Entry=x,score=scores.get(x),is_fragment=md[x]['is_explicit_fragment'],length=md[x].get('Length'),taxon=md[x]['taxonomic_bin']) for x in sorted(miss)],'additional_accessions':sorted(extra)}
errors=[r for r in read(O/'test_predictions.csv') if r['outcome'] in ['FP','FN']]
for r in errors:
 if r['Entry'] in rv:r['annotation_review']=rv[r['Entry']]['review_category']
stats['test_errors']=errors
(O/'full_scan_reanalysis.json').write_text(json.dumps(stats,indent=2))
for name,rows in [('full_scan_predictions.csv',[dict(Entry=x,score_bits=scores[x],PF00014=int(x in oldpos),domains=domcount[x]) for x in sorted(newids)]),('test_errors.csv',errors)]:
 with (O/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
print(json.dumps(stats,indent=2))
