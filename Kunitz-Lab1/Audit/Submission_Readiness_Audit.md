# Final submission readiness audit

Updated 18 November 2026 after a clean local reproduction. This checklist supersedes the earlier warning about an unavailable November 3 backup. Evidence is saved with the submission; the live Colab notebook was not used for this local execution.

| Status | Item | Evidence |
|---|---|---|
| PASS | Original workflow from raw inputs | All 71 code cells executed successfully in a new local working directory. The extra cell sets up the verified input archive. See `original_workflow_execution.json`, `original_workflow_execution.log` and the executed notebook. |
| PASS | Input archive and old dependency | The notebook verifies the November 7 archive SHA256 `8f253ab895334be8af6baf0334a4d748a588bdee8fc297dc8180a34606153598` and checks restored files against its manifest. No separate November 3 archive is required. The full submission includes the archive; GitHub links it separately. |
| PASS | Genuine recomputation | Only `data/raw/` inputs were restored. Structural selection, coordinate extraction, MUSTANG alignment, native hmmbuild, training consistency check, MMseqs2 grouping, validation/test searches, complete Swiss-Prot search, error review, taxonomy, figures and backups executed again. Processed tables, models and search outputs were not restored. |
| PASS | Structural extraction and rebuilt model | All 13 extracted PDB fragments, seed FASTA and Stockholm alignment were byte-identical to the archive. All rebuilt HMM parameters matched after excluding build-date metadata. The model file hash changes because the native builder writes a new date. See `reproduction_verification.json`. |
| PASS | Original evaluation | All sampled validation/test accessions, labels, similarity groups and decision scores matched. Original selected cutoff remained 8.7 bits; test TP174/FN0/FP0/TN1740. |
| PASS | Complete native search | All full-database protein classifications and domain membership at 8.7 and 17.2 bits matched. Protein counts426/398 and domain counts514/485 were reproduced. Both searches reported272598 sequence rows and314336 domain rows. A small number of printed fields differed between native CPU builds:20 sequence rows and27 domain rows. One weak best-domain score was -2.0 instead of -2.1; no classification changed. Raw outputs are therefore not claimed byte-identical. |
| PASS | Final-report expanded evaluation | The separately executed reanalysis retains threshold17.2 and test TP172/FN2/FP2/TN283816; sensitivity/precision98.85%, MCC0.9885. Earlier calculation-cell replay reproduced the saved JSON results exactly. Fresh native classification at17.2 agrees with the archive. |
| PASS | Full-scan report summaries | Archived-score recalculation:398 proteins,485 domains,52 multidomain;395 PF00014 plus3 additional;3 PF00014 missed. Taxonomic counts223/97/57/17/4 sum398. Fresh native positive and above-threshold domain memberships match. |
| PASS | Separation checks | The previous direct-overlap audit over566420 evaluation proteins found no training accession or exact complete seed-fragment occurrence; no shared validation/test accession, exact full-sequence hash or retained group. The reproduction retained the original sampled split. |
| PASS | Pfam comparison | Both fixed models were independently searched in the comparison. All283992 structural test decisions matched. PF00014.31 calibration20.3: test TP174/FN0/FP1/TN283817. Published GA21 was applied to native scores; validation Q8WPG5 native20.956762 was correctly below21 despite its rounded21.0 display. |
| PASS | Structural review | All78 pairwise Cα fits reproduced archived RMSDs; mean0.920Å. Native ChimeraX1.11.1 Matchmaker and independent coordinate checks support Figures7–8. No fragment was removed on the basis of this review. |
| PASS | Report contents | Current report: five pages, eight figures, two tables, structured abstract under150 words, introduction, methods, results, discussion, conclusion and references. Lecturer brief pp15–23 covers these core tasks and requires appropriate citations. |
| WARNING | Scientific scope | PF00014 labels are an annotation reference. Exact overlap exclusion does not remove all homologous relationships. Expanded evaluation is retrospective; reference-model seed overlap was not excluded. These disclosed limits were not altered by the execution repair. |
| INFO | Tested environment | Local macOS execution used MUSTANG3.2.4, HMMER3.4, MMseqs2 15.6f452 and Biopython1.88. The Ubuntu archive records MMseqs2 15-6f452+ds-2. Native package/build differences explain why all raw printed fields are not required to match. Fresh Colab installation/authorization was not tested. |

## Reproduction files

- `../REPRODUCING.md` gives the archive checksum, download link, local command and Colab setup.
- `../Code/environment.yml` lists the software environment.
- `../Code/run_original_notebook.py` executes every code cell and records failures.
- `../Code/verify_original_reproduction.py` compares a completed clean run against the original archive.
- `original_workflow_execution.json` records all71 successful cells.
- `reproduction_verification.json` records the model, input, membership and classification comparisons.

No missing private historical snapshot remains in the documented local reproduction path. The required source archive is included in the complete submission and linked from the public repository. Submission portal size/extension rules remain a submission-time administrative check. The audit establishes reproducibility within its stated scope; it does not establish an external independent biological test.
