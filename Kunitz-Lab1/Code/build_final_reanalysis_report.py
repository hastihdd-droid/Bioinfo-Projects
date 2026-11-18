from pathlib import Path
import json,csv,sys,os
sys.path.insert(0,str(Path('work/report_packages').resolve()))
os.environ['MPLCONFIGDIR']=str(Path('work/mplconfig').resolve())
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from docx import Document
from docx.shared import Cm,Pt
r=json.loads(Path('work/reanalysis/results.json').read_text());f=json.loads(Path('work/reanalysis/full_scan_reanalysis.json').read_text())
figs=Path('work/final_reanalysis_figures');figs.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
# Actual validation threshold sweep.
with open('work/reanalysis/threshold_sweep.csv') as h:sw=list(csv.DictReader(h))
x=[(float(a['lower'])+float(a['upper']))/2 for a in sw];y=[float(a['MCC']) for a in sw]
fig,ax=plt.subplots(figsize=(4,2.6));ax.step([float(a['lower']) for a in sw]+[float(sw[-1]['upper'])],y+[y[-1]],where='post',color='#376d91',lw=1.7);ax.axvline(17.2,color='#c38235',ls='--',label='Selected: 17.2 bits');ax.set(xlim=(0,40),ylim=(.94,1.002),xlabel='Best-domain score threshold (bits)',ylabel='Validation MCC');ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=9);fig.tight_layout();fig.savefig(figs/'threshold.png',dpi=250);plt.close(fig)
fig,ax=plt.subplots(figsize=(4,2.8));a=np.array([[283816,2],[2,172]]);ax.imshow([[.4,.05],[.05,.4]],cmap='Blues',vmin=0,vmax=1)
for (i,j),v in np.ndenumerate(a):ax.text(j,i,f'{v:,}',ha='center',va='center',weight='bold',fontsize=14)
ax.set(xticks=[0,1],xticklabels=['Negative','Positive'],yticks=[0,1],yticklabels=['Negative','Positive'],xlabel='Predicted class',ylabel='Reference class');fig.tight_layout();fig.savefig(figs/'test.png',dpi=250);plt.close(fig)
fig,ax=plt.subplots(figsize=(4,2.6));ax.bar([0,1,2],[395,3,3],color=['#376d91','#aeb5bb','#c38235']);ax.set_xticks([0,1,2],['PF00014\ndetected','PF00014\nmissed','Detected without\nPF00014']);ax.set_ylabel('Proteins');ax.set_ylim(0,440)
for i,v in enumerate([395,3,3]):ax.text(i,v+6,str(v),ha='center')
fig.tight_layout();fig.savefig(figs/'full_scan.png',dpi=250);plt.close(fig)
d=Document('outputs/Elham_Hadadi_Lab1_Report_Revised.docx');p=d.paragraphs
p[3].text=''
parts=[('Motivation: ','Can a structure-derived profile HMM detect Kunitz-related domains in Swiss-Prot? '),('Results: ','A 53-state HMM was built from 13 structural fragments. An initial sampled-negative benchmark gave error-free classification at 8.7 bits. A retrospective reanalysis expanded evaluation to all 566,072 nonfragment reference negatives while retaining the original positive split. Validation MCC selected a threshold of 17.2 bits. The test contained 172 true positives, two false negatives, two false positives and 283,816 true negatives, giving sensitivity and precision of 98.85% and MCC of 0.9885. Applying this threshold to complete Swiss-Prot identified 398 proteins: 395 with PF00014 annotation and three additional candidates. '),('Availability and Implementation: ','The original Colab notebook and a separate reproducible reanalysis notebook document the workflows. The supplementary archive contains the verified inputs, calculation and outputs.')]
for k,v in parts:p[3].add_run(k).bold=True;p[3].add_run(v)
assert len(p[3].text.split())<=150
p[15].text='There were initially 368 positive candidates. Training-accession overlap and exact training-fragment matches excluded 15 proteins; a further five met a local similarity criterion of at least 90% identity and 90% seed coverage. This left 348 positives. The initial benchmark sampled 3,480 of 566,072 reference negatives with seed 42. To examine dependence on this sampling, the reanalysis used the complete nonfragment negative pool.'
p[16].text='The original 174 validation and 174 test positives were retained. Original groups based on exact duplicates and MMseqs2 matches at 80% identity over at least 40 paired standard residues remained intact [3]. For added negatives, exact complete-sequence duplicates were grouped globally; original negative groups kept their previous allocation. Each new group was assigned using the parity of SHA256 of seed 42 and its smallest accession. This gave 282,254 validation and 283,818 test negatives. No conflicting-label exact-sequence groups were found. Added negatives were not subjected to a new homology clustering step; lower-identity relationships can remain across the split.'
p[18].text='The decision score was the highest reported domain bit score; absent hits were negative. The archived hmmsearch runs used --max --noali -E 1000000 --domE 1000000, --seed 42 and --cpu 2. The broad reporting limits were separate from the classification threshold. Raw tblout scores were reparsed, with UniProt identifiers normalized to accessions. All scores for the original evaluation proteins agreed between their original searches and the complete scan.'
p[19].text='The reanalysis protocol was recorded before assigning the additional negatives or calculating its metrics. The model remained unchanged. Validation MCC selected the midpoint of the widest optimal interval, (13.6, 20.8] bits, giving 17.2 bits. This threshold was recorded before applying it to the reanalysis test set. Because the original results had already been inspected, this is a retrospective reanalysis, not a new external test. No alternative split seeds were tried to obtain a preferred result.'
p[21].text='3.1 Validation and test results'
p[23].text='Figure 1. Test confusion matrix at 17.2 bits. Two reference positives were missed and two reference negatives were predicted positive.'
p[24].text='Validation gave 173 true positives, one false negative, no false positives and 282,254 true negatives (MCC 0.9971). The fixed 17.2-bit threshold on test gave sensitivity 172/174 = 98.85%, precision 172/174 = 98.85%, specificity 283,816/283,818 = 99.9993% and MCC 0.9885.'
p[25].text='The initial sampled-negative test had no errors at 8.7 bits. Reusing 8.7 bits on the expanded test retained all positives but produced nine reference false positives (precision 95.08%, MCC 0.9751). Thus, the change in sensitivity reflects recalibration on a larger validation background, not a correction to the original arithmetic.'
p[26].text='3.2 Complete Swiss-Prot search'
p[27].text='Applying the recalibrated 17.2-bit rule to the recorded search of all 575,748 Swiss-Prot proteins identified 398 positive proteins with 485 domains above threshold; 52 proteins had multiple domains. The model and raw search output were unchanged. This full scan includes training and evaluation proteins and is not an independent test.'
p[29].text='Figure 2. Complete-database comparison at 17.2 bits: 395 of 398 PF00014-annotated proteins detected, three missed and three additional predictions.'
p[30].text='The 398-protein PF00014 reference included 368 nonfragments and 30 fragments. Three annotated nonfragment proteins were missed: A0A1Q1NL17 (11.2 bits), D3GGZ8 (11.3) and O62247 (13.9). Among the 398 predictions, 395 had PF00014 annotation (99.25% annotation agreement). This is not a definitive estimate of biological precision.'
p[31].text='3.3 Threshold selection and error review'
p[32].text='Figure 3 shows validation MCC across score thresholds. The higher cutoff reduced reference false positives at the cost of missing low-scoring annotated positives. The two test false negatives were A0A1Q1NL17 and O62247; D3GGZ8 was the validation false negative. Their scores are below the fixed cutoff, so the observed misses are consistent with the stated rule.'
p[34].text='Figure 3. Validation MCC in the expanded-negative reanalysis. The dashed line marks the midpoint selected from the optimal interval. Test results were not used to select this cutoff.'
p[35].text='The two test false positives relative to PF00014 were P56409 (ornithodorin, 18.0 bits) and Q09JW3 (monogrin 2, 17.2). Both had BPTI/Kunitz domain evidence in the archived UniProt records despite lacking PF00014. They remain false positives under the fixed reference definition; this does not establish that they are biological false matches. The third additional full-scan candidate, P84555, was an explicitly annotated Kunitz fragment. These records provide existing annotation support, not new experimental confirmation.'
p[37].text='At 17.2 bits, predictions comprised 223 vertebrate, 97 arthropod, 57 cnidarian, 17 other metazoan and four nematode proteins. These mutually exclusive bins are descriptive and are not equivalent taxonomic ranks. All predictions were in metazoan bins at this cutoff; the weaker fungal and bacterial predictions from the initial 8.7-bit scan did not remain.'
p[38].text='There were 31 explicitly flagged fragments among the predictions: 30 with PF00014 and P84555 without it. Fragment status was excluded from evaluation datasets but retained in the complete-database application.'
p[40].text='Using the complete negative pool exposed the dependence of threshold calibration on background sampling. The original perfect benchmark was numerically correct, but random negative sampling missed rare reference-negative matches. The larger validation background selected a more conservative threshold, improving test precision while reducing sensitivity.'
p[41].text='The main limitation remains the annotation-based reference. PF00014 absence does not prove biological absence, as illustrated by the two test false positives with other Kunitz evidence. These labels were not changed after reviewing errors. The model detects family-related sequence similarity; it does not establish inhibitory activity or complete domain boundaries.'
p[42].text='The reanalysis is retrospective and reuses the positive split. The original similarity rule does not remove all homologous relationships: only nine test positives had less than 40% identity to their closest coverage-qualified training match. Seven of these nine were detected at 17.2 bits, compared with all nine at 8.7 bits. This small subset illustrates the sensitivity trade-off but cannot establish broad remote-homology performance. Newly added negatives were grouped for exact duplication but not exhaustively for homology. A future external dataset would be needed for a fresh independent assessment.'
p[44].text='Expanding negative evaluation data changed the validation-selected threshold from 8.7 to 17.2 bits. At the recalibrated threshold, test sensitivity and precision were 98.85% with MCC 0.9885. The complete Swiss-Prot application identified 398 candidates, including 395 PF00014-annotated proteins and three additional predictions. The model remains useful for candidate detection, with explicit trade-offs between missed positives and annotation-defined false positives.'
p[46].text='The original analysis used Biopython, MUSTANG 3.2.4, HMMER 3.4 and MMseqs2 15-6f452+ds-2 in Colab, using the verified Swiss-Prot 2026_03 archive. Its 70 code cells completed in the recorded run. The subsequent reanalysis reparsed the archived HMMER output and verified extracted files against the backup manifest; HMMER and the structural model were not rerun or retrained. The supplementary notebook records the fixed protocol, split, recalculation and updated full-scan results. Original results are retained for comparison.'
# Resize the existing dataset table, maintaining all counts.
t=d.tables[0]
t.rows[1].cells[2].text='282,254';t.rows[2].cells[2].text='283,818'
for row in t.rows:
 for c in row.cells:
  for pp in c.paragraphs:
   for run in pp.runs:run.font.name='Times New Roman';run.font.size=Pt(9)
# Preserve three figures, replace with verified updated results.
from PIL import Image
for shape,name in zip(d.inline_shapes,['test.png','full_scan.png','threshold.png']):
 b=shape._inline.graphic.graphicData.pic.blipFill.blip
 d.part.related_parts[b.embed]._blob=(figs/name).read_bytes()
 with Image.open(figs/name) as im:ratio=im.height/im.width
 shape.width=Cm(8.25);shape.height=int(shape.width*ratio)
for para in d.paragraphs:
 if para.style.name=='Normal' and para.text and not para.text.startswith('Elham Hadadi'):
  para.paragraph_format.space_after=Pt(4)
# Current hyperlink to original notebook remains, labelled as original workflow.
for h in d._element.body.iter():
 if h.tag.endswith('}t') and h.text=='Open the executed Colab notebook':h.text='Open the original executed Colab notebook'
d.save('outputs/Elham_Hadadi_Lab1_Final_Reanalysis.docx')
print('Abstract words:',len(p[3].text.split()))
