# Pfam PF00014 reference comparison

Elham Hadadi — Laboratory of Bioinformatics I

## Purpose and design
The existing structural model is compared with the official InterPro/Pfam PF00014.31 profile, downloaded on 18 November 2026 (InterPro 111.0). The reference has 53 match states and 98 seed sequences. The structural model is unchanged: 53 match states, 13 structural seed fragments and fixed 17.2-bit best-domain threshold.

Both models were searched afresh against the same 566,420 existing evaluation proteins, preserving the Swiss-Prot 2026_03 labels and validation/test membership. PyHMMER 0.12.3 was used with bias filtering disabled, F1=F2=F3=1, E=domE=1000000, seed42 and two configured workers. Search ran from disk to limit memory usage. All 283,992 structural test decisions exactly reproduced the archived classification.

The primary Pfam threshold was selected solely from expanded validation MCC. The widest optimal interval was (19.6, 21.0], with midpoint 20.3 bits; validation had 174 TP, 0 FN, 0 FP and 282,254 TN. This cutoff was saved before computing test metrics. No structural threshold was retuned.

## Test results
Test membership: 174 positives and 283,818 negatives.

| Model | Domain cutoff | TP | FN | FP | Sensitivity | Precision | MCC |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Structural HMM | 17.2 | 172 | 2 | 2 | 98.85% | 98.85% | 0.9885 |
| Pfam validation-calibrated | 20.3 | 174 | 0 | 1 | 100.00% | 99.43% | 0.9971 |
| Pfam curated GA | 21.0 | 174 | 0 | 0 | 100.00% | 100.00% | 1.0000 |

The primary calibrated comparison favors Pfam on this reference: it recovers both structural-model misses (A0A1Q1NL17 andO62247), rejects Q09JW3 and retains P56409 as an annotation-defined FP. Full scores and differences are included in computed/.

## Published gathering cutoffs
Pfam supplies both sequence and domain GA cutoffs of 21 bits. This secondary rule was applied without tuning. It gave 173 TP / 1 FN / 0 FP / 282,254 TN on validation and 174 TP / 0 FN / 0 FP / 283,818 TN on test. A boundary check was necessary: Q8WPG5 is printed as 21.0 bits but has native sequence/domain score 20.956762313842773, below 21. This validation example is a GA miss. Native boundary scores are included; the primary one-decimal midpoint comparison is unchanged.

## Interpretation
This is a comparison against an established reference model, not an external or newly independent test. PF00014 annotations define the labels, so close agreement by the Pfam reference is expected and is not experimental biological validation. Possible overlap between the Pfam seed and benchmark was not excluded. The main benchmark's homology and retrospective limitations remain. The result does not justify claims that the models are equivalent or that fewer training sequences are intrinsically superior.

## Recalculate
Run `python recompute_comparison.py` using Python 3, or open Pfam_Comparison_Executed.ipynb with this folder as its working directory and Run all. The calculation needs only the standard library. Inputs include compressed raw tblout for both models, benchmark membership and original decisions; outputs are written under computed/. The portable calculation was executed and exactly reproduced all saved comparison results.

## Repeat the searches
The included search_models.py and check_ga_boundary.py are original run sources and retain the task's work/ directory layout. Full search reproduction additionally needs the original Swiss-Prot 2026_03 FASTA from the project backup, existing assignments.csv, the two HMM files and PyHMMER 0.12.3. See search_record.json for hashes and settings. No executable from the project backup was run. The complete source database is linked from the project's main README rather than duplicated here.

## Sources
- Model: https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00014/?annotation=hmm
- Metadata: https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00014/
- Search settings: https://pyhmmer.readthedocs.io/en/stable/api/plan7/pli.html
