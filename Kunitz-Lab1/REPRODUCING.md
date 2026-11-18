# Reproducing the original Kunitz workflow

Use `Notebooks/Lab1.ElhamHadadi.snapshot-fixed.ipynb` and the original input backup:

`kunitz_project_backup_20261107T213430_789701Z.zip`

SHA256: `8f253ab895334be8af6baf0334a4d748a588bdee8fc297dc8180a34606153598`

Download: https://drive.google.com/file/d/1WBFyxo18TtFOom_M_cH_BvdCJNNUbiA9/view

The complete submission ZIP includes this backup in `Supplementary/`. GitHub links it separately because of its size. The November 3 backup is not needed.

## Local execution

Create and activate the environment, then run every notebook cell in a new directory:

```sh
conda env create -f Code/environment.yml
conda activate kunitz-lab1
python Code/run_original_notebook.py --archive /absolute/path/kunitz_project_backup_20261107T213430_789701Z.zip --workdir /absolute/path/new_kunitz_run
```

To verify the completed run against the archive:

```sh
python Code/verify_original_reproduction.py --archive /absolute/path/kunitz_project_backup_20261107T213430_789701Z.zip --workdir /absolute/path/new_kunitz_run
```

The runner requires an empty working directory. It saves an executed notebook and `execution.json`. The original workflow has 71 code cells, including the added archive setup. Expected original sampled-test result: threshold 8.7 bits; TP174, TN1740, FP0, FN0. The original complete scan at this cutoff returns 426 proteins. These are historical original-workflow results; the report's 17.2-bit expanded evaluation is reproduced by the separate reanalysis notebook.

The setup verifies the complete archive SHA256, verifies individual restored files against its manifest, and restores raw inputs only. It does not restore processed selection tables, alignments, models, scores or result pointers. Archived UniProt responses are selected by recorded request and release 2026_03. Raw PDB coordinates and the original PDB candidate search keep the structural input fixed. All later calculations and native MUSTANG, hmmbuild, MMseqs2 and hmmsearch commands execute again. A local run writes the two backup archives locally instead of mounting Google Drive.

## Colab

Open the repaired notebook. Put the November 7 ZIP in `/content` or the Drive `MyDrive/kunitz_project_backups` folder, then Run all. Alternatively, set `KUNITZ_SNAPSHOT` to the archive's full path before running. Colab already supplies the Python notebook environment; the original installation cells install missing native tools with apt. Drive mounting requires your authorization. The verification record describes the tested local environment; it does not claim a fresh Colab session was tested.

## Separate final-report analyses

Run `Notebooks/Lab1_ElhamHadadi_Reanalysis_Executed.ipynb` for the expanded-negative evaluation. It uses archived scores, retains the structural HMM, selects 17.2 bits on validation, and yields test TP172/FN2/FP2/TN283816 and 398 full-scan predictions. The Pfam and structural-analysis supplementary folders include their own execution instructions.

The clean-run verification compares newly extracted fragments, alignment, model parameters and raw search rows with the original archive. HMM file hashes differ because the native builder writes the new build date; the model parameters must match after excluding that metadata line. No evaluation threshold or biological dataset is changed by the reproducibility repair.
