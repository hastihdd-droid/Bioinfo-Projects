from pathlib import Path
import os,json,hashlib,zipfile,contextlib,io
base=Path.cwd();run=base/'work/submission_audit_run'; w=run/'work';w.mkdir(parents=True,exist_ok=True)
for name in ['raw_audit','report_packages']:
 p=w/name
 if not p.exists():p.symlink_to(base/'work'/name,target_is_directory=True)
(run/'outputs').mkdir(exist_ok=True)
(w/'expanded_design.json').write_bytes((base/'work/expanded_design.json').read_bytes())
os.chdir(run)
log=[]
ns={'__name__':'__main__','Path':Path}
nb=json.loads((base/'outputs/Lab1_ElhamHadadi_Reanalysis_Executed.ipynb').read_text());cells=[c for c in nb['cells'] if c['cell_type']=='code']
# Replay actual calculation cells; extraction inspected separately with archive manifest.
for i,c in enumerate(cells[1:],2):
 buf=io.StringIO()
 with contextlib.redirect_stdout(buf):exec(compile(''.join(c['source']),f'notebook_cell_{i}','exec'),ns)
 log.append(f'Cell {i}: PASS\n'+buf.getvalue())
for script in ['final_consistency_check.py','check_direct_training_overlap.py']:
 buf=io.StringIO()
 with contextlib.redirect_stdout(buf):exec(compile((base/'work'/script).read_text(),script,'exec'),ns)
 log.append(script+':\n'+buf.getvalue())
for name in ['results.json','selected_threshold.json','full_scan_reanalysis.json']:
 assert json.loads((w/'reanalysis'/name).read_text())==json.loads((base/'work/reanalysis'/name).read_text()),name
archive=Path('/Users/hasti/Downloads/kunitz_project_backup_20261107T213430_789701Z.zip')
with zipfile.ZipFile(archive) as z:
 manifest=json.loads(z.read('backup_file_manifest.json'));verified=0
 for e in manifest['files']:
  name='kunitz_project/'+e['relative_path']
  with z.open(name) as f:h=hashlib.file_digest(f,'sha256').hexdigest()
  assert h==e['sha256'] and z.getinfo(name).file_size==e['size_bytes'],name
  verified+=1
log.append(f'Full original archive: {verified} files verified against size and SHA256 manifest.')
(base/'work/submission_audit_execution.log').write_text('\n'.join(log))
print('\n'.join(x.splitlines()[0] for x in log));print('Replayed calculation results exactly equal saved results.')
