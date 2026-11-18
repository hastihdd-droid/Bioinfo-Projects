import csv,json,math,hashlib
from pathlib import Path
O=Path('work/reanalysis');R=Path('work/raw_audit')
r=json.loads((O/'results.json').read_text());t=r['selected_threshold']['threshold_bits'];assert t==17.2
with (O/'test_predictions.csv').open() as f:rows=list(csv.DictReader(f))
from collections import Counter
c=Counter(x['outcome'] for x in rows);assert c=={'TP':172,'FN':2,'FP':2,'TN':283816}
for x in rows:
 score=float(x['score_bits']) if x['score_bits'] else -math.inf
 assert int(x['predicted_label'])==int(score>=t)
assert len(rows)==283992
with (O/'assignments.csv').open() as f:a=list(csv.DictReader(f))
assert len(a)==566420 and len({x['Entry'] for x in a})==566420
for column in ['Entry','group','sequence_sha256']:
 v={x[column] for x in a if x['split']=='validation'};s={x[column] for x in a if x['split']=='test'}
 assert not(v&s)
f=json.loads((O/'full_scan_reanalysis.json').read_text());assert f['predicted_proteins']==f['PF00014_detected']+f['additional_predictions']==398
assert f['PF00014_detected']+f['PF00014_missed']==398
assert sum(f['taxonomic_counts'].values())==398
with (R/'results/robustness/training_audit_20261107T213423_499191Z/positive_training_similarity_audit.csv').open() as h:tr={x['Entry']:x for x in csv.DictReader(h)}
remote=[x for x in rows if x['reference_label']=='1' and float(tr[x['Entry']]['identity_percent'])<40]
assert len(remote)==9 and sum(x['predicted_label']=='1' for x in remote)==7
print('PASS: all reanalysis predictions match the frozen threshold; counts and full-scan summaries reconcile.')
print('PASS: no shared accessions, connected groups or exact complete-sequence hashes across the new split.')
print('PASS: 7/9 test positives below 40% training identity detected at 17.2 bits.')
print('Scope: archived raw scores were reprocessed; HMMER was not rerun, and distant homology is not excluded.')
