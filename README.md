# LM_Bioinformatics_HMM_KunitzDomain
This is a computational pipeline for building a profile Hidden Markov Model (HMM) for the Kunitz-type protease inhibitor domain (Pfam ID: PF00014). Kunitz domains are known to inhibit the proteolytic activity of a protein family called proteases. Some examples of Kunitz-type protease inhibitors include aprotinin (bovine pancreatic trypsin inhibitor, BPTI), Alzheimer's amyloid precursor protein (APP), and tissue factor pathway inhibitor (TFPI).
An HMM model was constructed and calibrated using structure-based multiple alignemnt (PDBeFOLD) of representative Kunitz-type protease inhibitors.To evaluate its performance, the model was tested against both positive sequences (known Kunitz domains) and negative controls (non-Kunitz proteins).

This project was carried out as part of the Laboratory of Bioinformatics 1 during my MSc in Bioinformatics at the University of Bologna (Alma Mater Studiorum), aiming  to integrate approaches from structural bioinformatics, sequence analysis, and statistical validation of predictive models. 

## Needed Packages
We recommend using [Conda](https://docs.conda.io/en/latest/) to manage the project environment and dependencies. The following packages were installed:

**CD-HIT** – for clustering and comparing large sets of sequences
```bash
conda install -c bioconda cd-hit
```
**HMMER** – for searching sequence databases using profile hidden Markov models
```bash
conda install -c bioconda hmmer
```

**BLAST+** – for performing sequence similarity searches
```bash
conda install -c bioconda blast
```

**Biopython** – for working with biological data in Python
```bash
conda install -c conda-forge biopython
```

This setup ensures compatibility and simplifies dependency management throughout the project.

## Repository contents (da cambiare il titolo)
file
file
file (con descrizione)
.
.
.
.
.

## Workflow
1
2
3
4
.
.
.

## Outputs

file
file
file
file (descrizione)
.
.
.
.
.
.

