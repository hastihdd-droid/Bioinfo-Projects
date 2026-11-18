"""Execute every original-workflow cell in a new local working directory.
Usage: python run_original_notebook.py --archive /path/to/Nov7.zip --workdir /path/to/new_run
Run inside the environment defined by environment.yml. No existing processed results are loaded.
"""
import argparse, json, os, time, traceback
from pathlib import Path

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--notebook',type=Path)
    args=parser.parse_args()
    archive=args.archive.expanduser().resolve()
    if not archive.is_file():parser.error('Archive not found')
    run=args.workdir.expanduser().resolve()
    if run.exists() and any(run.iterdir()):parser.error('Use an empty working directory')
    run.mkdir(parents=True,exist_ok=True)
    notebook=args.notebook or Path(__file__).resolve().parents[1]/'Notebooks/Lab1.ElhamHadadi.snapshot-fixed.ipynb'
    os.environ['KUNITZ_CONTENT_ROOT']=str(run/'content')
    os.environ['KUNITZ_SNAPSHOT']=str(archive)
    os.environ['MPLCONFIGDIR']=str(run/'matplotlib')
    os.environ['IPYTHONDIR']=str(run/'ipython')
    os.environ['MPLBACKEND']='Agg'
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    import nbformat
    nb=nbformat.read(notebook,as_version=4)
    shell=InteractiveShell.instance()
    log=[];count=0
    for index,cell in enumerate(nb.cells):
        if cell.cell_type!='code':continue
        count+=1;start=time.time()
        print(f'Executing cell {index} ({count})',flush=True)
        with capture_output() as captured:
            result=shell.run_cell(cell.source,store_history=True)
        outputs=[]
        if captured.stdout:outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
        if captured.stderr:outputs.append(nbformat.v4.new_output('stream',name='stderr',text=captured.stderr))
        for item in captured.outputs:outputs.append(nbformat.v4.new_output('display_data',data=item.data,metadata=item.metadata))
        error=result.error_before_exec or result.error_in_exec
        if error:outputs.append(nbformat.v4.new_output('error',ename=type(error).__name__,evalue=str(error),traceback=traceback.format_exception(error)))
        cell.outputs=outputs;cell.execution_count=count
        log.append({'cell_index':index,'execution_count':count,'seconds':round(time.time()-start,3),'status':'FAIL' if error else 'PASS'})
        nbformat.write(nb,run/'Original_Workflow_Executed.ipynb')
        (run/'execution.json').write_text(json.dumps(log,indent=2)+'\n')
        if error:raise RuntimeError(f'Cell {index} failed: {error}') from error
    print(f'PASS: all {count} code cells executed. See {run}',flush=True)

if __name__=='__main__':main()
