# Interferon promoter motif conservation and comparative regulatory analysis

This repository contains the comparative promoter-analysis pipeline used to examine predicted transcription-factor motif conservation across mammalian *IFNE*, *IFNA1* and *IFNB* regulatory regions.

The workflow uses **hg38-anchored UCSC Multiz alignments**, applies species-level alignment QC, extracts orthologous promoter sequence for retained taxa, performs motif scanning with FIMO, merges motif hits across species and transcription factors, and generates comparative prevalence, phylogenetic presence/absence and positional-conservation plots.

> Requirements include **MEME Suite** (including FIMO), Python, Biopython, pandas, NumPy, matplotlib and motif PWMs from sources such as JASPAR and HOCOMOCO.

## Overview

The analysis proceeds as:

```text
UCSC Multiz MAF
    ↓
MAF QC
    ↓
QC-passing species list
    ↓
Gap-free FASTA extraction
    ↓
FIMO motif scanning
    ↓
Merge motif hits across species and TFs
    ↓
Comparative IFNE / IFNA1 / IFNB plotting
```

The promoter sequences represent genomic regions **orthologous to the human promoter interval**. They are not defined independently relative to species-specific transcription start sites.

## Input promoter alignments

Promoter regions are obtained from hg38-anchored UCSC Multiz MAF alignments for:

- *IFNE*
- *IFNA1*
- *IFNB*

The analysed interval is defined relative to the annotated human transcription start site and projected across species through the Multiz alignment.

Example files:

```text
data/maf/IFNE_100way.maf
data/maf/IFNA1_100way.maf
data/maf/IFNB_100way.maf
```

## 1. Quality-control the MAF alignments

Before extracting FASTA sequences, run the MAF QC script independently for each promoter.

The QC evaluates each species across alignment columns in which the human reference contains sequence.

Current exclusion criteria are:

- coverage relative to the hg38-aligned interval < 0.50;
- N fraction > 0.10;
- callable-window fraction < 0.30.

Callable windows are 100 bp windows containing:

- at least 70% ungapped sequence;
- no more than 10% Ns.

Long contiguous gaps are recorded as a warning but do not, by themselves, cause species exclusion.

Example:

```bash
python maf_qc.py \
  --maf ../data/maf/IFNE_100way.maf \
  --out-prefix ../analysis/qc/ifne_promoter_qc \
  --write-aligned-fasta
```

Run the same procedure separately for *IFNA1* and *IFNB*.

The QC script produces:

```text
*.metrics.tsv
*.excluded.tsv
*.kept.txt
*.aligned.fa
```

The `*.kept.txt` file is used as the whitelist for downstream FASTA extraction.

## 2. Parse MAF files into species FASTAs

Use `parse_maf.py` to extract gap-free orthologous sequence for species that passed QC.

Example for *IFNE*:

```bash
python parse_maf.py \
  -i ../data/maf/IFNE_100way.maf \
  -o ../analysis/fasta/IFNE \
  --keep-species ../analysis/qc/ifne_promoter_qc.kept.txt
```

Equivalent commands should be run for *IFNA1* and *IFNB*.

`parse_maf.py` removes and recreates the specified output directory before writing new FASTA files, ensuring that FASTAs from previously excluded taxa do not remain in the directory.

Each output FASTA contains the gap-free sequence recovered from the relevant MAF blocks for that species.

## 3. Motif files

PWM motif files should be placed in:

```text
motifs/
```

Motifs may be obtained from sources such as:

- JASPAR
- HOCOMOCO

Where equivalent motifs are present in multiple databases, downstream merging should retain a clear TF-family label and avoid interpreting database-specific motif duplication as independent biological evidence.

## 4. FIMO scanning

Run motif scans for each retained species and promoter.

Example:

```bash
bash scripts/scan_all.sh
```

FIMO output should retain, where available:

- species
- motif identifier / TF label
- sequence name
- start
- end
- strand
- score
- p-value
- q-value

The same motif panel should be applied to *IFNE*, *IFNA1* and *IFNB* wherever possible.

## 5. Merge FIMO results

FIMO hits are first combined across species and transcription factors within each promoter dataset.

The final structure is:

```text
IFNE_fimo_combined.tsv
IFNA1_fimo_combined.tsv
IFNB_fimo_combined.tsv
```

Each combined table contains all retained species and all scanned transcription-factor motifs for that promoter.

The plotting workflow expects, at minimum:

```text
species
motif / TF identifier
p-value
```

and, for positional analysis:

```text
start
end
```

Optional but recommended fields include:

```text
q-value
sequence_name
strand
```

## 6. Comparative promoter figure

Generate the comparative promoter figure with:

```bash
python make_promoter_figure.py \
  --config ../config/config.yaml \
  --ifne-maf ../data/maf/IFNE_100way.maf \
  --ifne-fimo ../analysis/tables/IFNE_fimo_combined.tsv \
  --ifna-maf ../data/maf/IFNA1_100way.maf \
  --ifna-fimo ../analysis/tables/IFNA1_fimo_combined.tsv \
  --ifnb-maf ../data/maf/IFNB_100way.maf \
  --ifnb-fimo ../analysis/tables/IFNB_fimo_combined.tsv \
  --tree ../data/trees/hg38_100way.nh \
  --ucsc-json ../data/trees/ucsc_genomes_official.json \
  --out-prefix ../analysis/figures/IFNE_promoter_figure
```

The figure contains:

- **Panel A:** motif prevalence across *IFNE*, *IFNA1* and *IFNB*;
- **Panel B:** phylogenetic species-by-motif presence/absence matrices;
- **Panel C:** positional motif maps within the *IFNE* orthologous promoter interval.

For direct comparison between promoters, Panels A and B should use the shared set of QC-passing species present across all three promoter datasets.

## Motif groups

The current analysis spans motifs representing several regulatory programmes.

### ETS / epithelial regulation

- ELF1
- ELF3
- ELF5
- ETS1

### Steroid-hormone regulation

- PGR
- AR
- ESR1
- ESR2

### Antiviral / inflammatory regulation

- IRF1
- IRF7
- IRF9
- NF-κB / RELA-family motifs

### Cancer / proliferation-associated candidates

Exploratory screens included:

- TP53-family motifs
- E2F1 / E2F4
- TEAD1–4
- AP-1 / JUN / FOS
- MAFF / small-MAF / NRF2-associated motifs
- FOXO1 / FOXO3

### Trophoblast-associated candidates

Exploratory screens included:

- TFAP2C / AP-2-family motifs
- GATA-family motifs
- TEAD4
- ELF5
- CDX2
- EOMES
- ESRRB
- MAFK
- ASCL2
- SMAD-family motifs
- GCM1
- KLF5
- MSX2
- EPAS1 / HIF2A
- PPARG

Motif matches are interpreted as **predicted cis-regulatory elements**, not as evidence of transcription-factor binding or occupancy.

## Relevant QC parameters

Current parameters are:

```text
minimum coverage vs hg38:        0.50
maximum N fraction:              0.10
callable window size:            100 bp
minimum ungapped fraction:       0.70
maximum N fraction per window:   0.10
minimum callable-window fraction: 0.30
```

Long gap runs are recorded for QC purposes but are not used as an automatic exclusion criterion.

FIMO significance thresholds should be defined consistently in `config/config.yaml` and reported explicitly in the Methods.

_Last updated: 2026-11-17_
