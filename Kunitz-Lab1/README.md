# Detecting Kunitz domains with a profile HMM

**Elham Hadadi · Laboratory of Bioinformatics I · University of Bologna**

A structure-based approach to identifying Kunitz domains in protein sequences. This project builds a profile hidden Markov model (HMM) from selected PDB structures, evaluates its predictions against Swiss-Prot annotations, and examines the distribution of predicted domains across the database.

[**Read the report →**](Final_Report/Elham_Hadadi_Lab1_Visual_Report.pdf) · [Original notebook](Notebooks/Lab1.ElhamHadadi.snapshot-fixed.ipynb) · [Reanalysis notebook](Notebooks/Lab1_ElhamHadadi_Reanalysis_Executed.ipynb)

## Approach

1. Select structurally resolved Kunitz fragments, checking coordinate completeness, experimental quality and sequence redundancy.
2. Align the selected fragments with **MUSTANG** and build a profile HMM with **HMMER**.
3. Evaluate the model using reviewed, nonfragment Swiss-Prot proteins, with **PF00014 annotation** defining the reference labels.
4. Select a score threshold using validation MCC, apply it to the test set, and inspect false positives and false negatives.
5. Search complete Swiss-Prot and summarize the predicted proteins and their taxonomic distribution.

The final structural seed contains **13 fragments from 10 UniProt accessions**. Its alignment and HMM both have **53 positions**.

![Sequence alignment derived from the structural superposition of the 13 seed fragments](Supplementary/model_evidence_figures/structural_alignment.png)
*The MUSTANG-derived seed alignment. Cysteines are highlighted in yellow.*

## Results

The initial benchmark used a random sample of reference negatives. A retrospective reanalysis retained the same HMM and positive split while expanding evaluation to all **566,072 nonfragment reference negatives**. Validation MCC selected a best-domain score threshold of **17.2 bits**.

| Test-set measure | Result |
| :--- | ---: |
| Sensitivity | 98.85% |
| Precision | 98.85% |
| Specificity | 99.9993% |
| Matthews correlation coefficient (MCC) | 0.9885 |
| True positives / false negatives | 172 / 2 |
| True negatives / false positives | 283,816 / 2 |

Applying this threshold to **575,748 Swiss-Prot proteins** identified **398 proteins** with **485 above-threshold domains**. Of the predicted proteins, 395 had PF00014 annotation and three had other archived Kunitz evidence. Three PF00014-annotated proteins were missed.

![Taxonomic distribution of the 398 predicted proteins](Supplementary/model_evidence_figures/taxonomic_distribution.png)
*Counts describe this Swiss-Prot search; they are not normalized estimates of biological prevalence.*

## Files

| Folder | Contents |
| :--- | :--- |
| [Final_Report](Final_Report/) | Four-page report in PDF and Word format |
| [Notebooks](Notebooks/) | Executed original workflow and retrospective reanalysis |
| [Code](Code/) | Calculation, figure-generation and consistency-check scripts |
| [Supplementary](Supplementary/) | HMM, alignment, selected raw search outputs, result tables and figures |
| [Audit](Audit/) | Verification checklist and execution log |

## Reproducing the reanalysis

The reanalysis notebook recalculates results from saved HMMER outputs; it does not rebuild the model or rerun the searches.

1. Download the [original project backup](https://drive.google.com/file/d/1WBFyxo18TtFOom_M_cH_BvdCJNNUbiA9/view), named `kunitz_project_backup_20261107T213430_789701Z.zip`.
2. Place it in `~/Downloads` for a local Python 3 notebook session, or in `MyDrive/kunitz_project_backups` for Google Colab.
3. Open the [reanalysis notebook](Notebooks/Lab1_ElhamHadadi_Reanalysis_Executed.ipynb) in a fresh working directory and run its cells in order. It verifies the archived inputs and recalculates the threshold, evaluation metrics and complete-database results.

The backup is stored separately because of its size; Drive access may be required. Large tables in this repository are gzip-compressed. The supporting scripts use the directory layout documented in the notebooks.

The original workflow used **Biopython, MUSTANG 3.2.4, HMMER 3.4 and MMseqs2**. A fresh end-to-end execution of the original notebook also requires an older November 3 snapshot that is not included here. See the [reproducibility checklist](Audit/Submission_Readiness_Audit.md) for the exact dependency and verification scope.

## Interpretation and limitations

PF00014 is an annotation-based reference rather than experimental ground truth. The two test false positives have other Kunitz evidence, so disagreement with these labels does not necessarily indicate a biological false match.

Direct training-accession and exact seed-fragment overlap checks passed. Related sequences can nevertheless remain across the evaluation sets. The retrospective results therefore support performance on the defined dataset, rather than establishing performance on a new external dataset or on remote homologs. The complete-database search also includes training and evaluation proteins.
