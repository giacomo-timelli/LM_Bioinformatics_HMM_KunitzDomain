# LM_Bioinformatics_HMM_KunitzDomain
This is a computational pipeline for building a profile Hidden Markov Model (HMM) for the Kunitz-type protease inhibitor domain (Pfam ID: PF00014). Kunitz domains are known to inhibit the proteolytic activity of a protein family called proteases. Some examples of Kunitz-type protease inhibitors include aprotinin (bovine pancreatic trypsin inhibitor, BPTI), Alzheimer's amyloid precursor protein (APP), and tissue factor pathway inhibitor (TFPI).
An HMM model was constructed and calibrated using structure-based multiple alignemnt ([PDBeFOLD](https://www.ebi.ac.uk/msd-srv/ssm/)) of representative Kunitz-type protease inhibitors.To evaluate its performance, the model was tested against both positive sequences (known Kunitz domains) and negative controls (non-Kunitz proteins).

This project was carried out as part of the Laboratory of Bioinformatics 1 during my **MSc in Bioinformatics at the University of Bologna** (Alma Mater Studiorum), aiming  to integrate approaches from structural bioinformatics, sequence analysis, and statistical validation of predictive models. 

## Needed Packages
We recommend using [conda](https://docs.conda.io/en/latest/) to manage the project environment and dependencies. The following packages were installed:

**CD-HIT** – for clustering and reducing redundancy of protein sequences.
```bash
conda install -c bioconda cd-hit
```
**HMMER** – to build and search Hidden Markov Models (HMMs) for detecting protein domains.
```bash
conda install -c bioconda hmmer
```

**BLAST+** – for performing sequence similarity using blastp.
```bash
conda install -c bioconda blast
```

**Biopython** – necessary for executing get_seq.py, which performs sequence extraction.
```bash
conda install -c conda-forge biopython
```

This setup ensures compatibility and simplifies dependency management throughout the project.

## Workflow
### 1.Set up the conda environment
Make sure you have `conda` installed, create a new environment and activate it. Then, install the required dependencies using `conda` as specified in the [Required Packages](https://github.com/giacomo-timelli/LM_Bioinformatics_HMM_KunitzDomain/blob/main/README.md#needed-packages) section.
```bash
conda create --name my_env
conda activate my_env
```

### 2.Download RCSB PDB and Swiss-Prot datasets
Access the [RCSB PDB](https://www.rcsb.org/) database and type this query:
```
Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND Polymer Entity Sequence Length >= 45
```
> **Note:** After running the query, go to the **"Create Custom Report"** section and make sure to select the following attributes:  
> **PDB ID**, **Entity ID**, **Sequence**, **Auth Asym ID**, and **Annotation Identifier**.
>
>  Eventually download the *csv file*.


Access the [UniProt](https://www.uniprot.org/) database, filter for Swiss-Prot reviewed entries and download the complete set of protein sequences in FASTA format.
The downloaded file, named `uniprot_sprot.fasta`, should be placed in your working directory. It will be used to extract protein sequences that do not contain Kunitz domains.

### 3.Retrieve Representative Structural IDs
Run the script:
```bash
bash script_representative_kunitz.sh
bash filtered_fasta.sh
```
This will create the file  `tmp_pdb_kunitz_clstr_rep_clean_ids.txt`
> **Note:** Before submitting the list to PDBeFold, sequences are automatically filtered using the `filtered_fasta.sh` script, which removes entries that are too long, too short, or contain unstructured terminal regions. This ensures compatibility with PDBeFold and preserves the quality of the alignment and HMM.





















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


## Repository overview
file
file
file (con descrizione)
.
.
.
.
.


