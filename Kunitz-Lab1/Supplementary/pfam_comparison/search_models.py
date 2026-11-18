from pathlib import Path
import sys,csv,json,time,hashlib,inspect
sys.path.insert(0,str(Path('work/pfam_packages').resolve()))
import pyhmmer
P=Path('work/pfam_comparison');R=Path('work/raw_audit');start=time.time()
with open('work/reanalysis/assignments.csv') as f:ids={r['Entry'] for r in csv.DictReader(f)}
source=R/'data/raw/uniprot/swissprot_complete_2026_03/swissprot_complete.fasta';target=P/'evaluation.fasta';count=0
with source.open() as f,target.open('w') as out:
 header=None;seq=[]
 def save():
  global count
  if header is not None:
   entry=header.split()[0].split('|')[1]
   if entry in ids:out.write('>'+entry+'\n'+''.join(seq)+'\n');count+=1
 for line in f:
  if line.startswith('>'):save();header=line[1:].strip();seq=[]
  else:seq.append(line.strip())
 save()
assert count==len(ids)==566420
original=R/'results/hmm/build_20261107T212501_925161Z/kunitz_initial.hmm';assert hashlib.sha256(original.read_bytes()).hexdigest()=='3d315c6484a3ad992f3b15d52d881da790dffcb5f010139c8fb73aeca9cb0308'
files=[original,P/'PF00014.hmm'];models=[]
for p in files:
 with pyhmmer.plan7.HMMFile(p) as f:models.append(f.read())
print('Prepared',count,'identical evaluation targets; PyHMMER',pyhmmer.__version__,flush=True)
settings=dict(cpus=2,bias_filter=False,F1=1.,F2=1.,F3=1.,E=1000000,domE=1000000,seed=42)
record={'pyhmmer_version':pyhmmer.__version__,'settings':settings,'sequence_count':count,'protocol_sha256':hashlib.sha256((P/'protocol.json').read_bytes()).hexdigest(),'input_fasta_sha256':hashlib.file_digest(target.open('rb'),'sha256').hexdigest(),'model_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'reference_accession':models[1].accession,'reference_match_states':models[1].M,'reference_nseq':models[1].nseq,'gathering_cutoffs':list(models[1].cutoffs.gathering)}
(P/'search_record.json').write_text(json.dumps(record,indent=2))
for model,name in zip(models,['structural','pfam']):
 tick=time.time();print('Searching',name,flush=True)
 with pyhmmer.easel.SequenceFile(target,digital=True) as seqfile:
  hits=next(pyhmmer.hmmer.hmmsearch(model,seqfile,parallel='queries',**settings))
 with (P/(name+'.tbl')).open('wb') as f:hits.write(f,format='targets',header=True)
 print('Finished',name,'hits',len(hits),'seconds',round(time.time()-tick,1),flush=True)
record['elapsed_seconds']=round(time.time()-start,2);record['search_completed']=True
(P/'search_record.json').write_text(json.dumps(record,indent=2))
