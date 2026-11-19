# IFNE promoter motif conservation analysis

Comparative analysis of predicted transcription-factor motifs in orthologous *IFNE*, *IFNA1* and *IFNB* promoter regions across eutherian mammals.

The workflow uses hg38-anchored UCSC Multiz MAF alignments, removes poorly aligned or structurally problematic taxa, scans retained promoter sequences with FIMO, and generates comparative motif prevalence and positional-conservation plots.

## Pipeline

Run the steps below separately for *IFNE*, *IFNA1* and *IFNB*.

### 1. MAF structural sanity check

```bash
python maf_check.py ../data/maf/IFNE_100way_vertebrates_hg38.maf \
  --tsv-out ../analysis/qc/IFNE_maf_sanity.tsv
```

Repeat for IFNA1 and IFNB.

`FAIL` taxa are removed downstream. `WARN` taxa are retained.

---

### 2. Promoter alignment QC

```bash
python maf_qc.py \
  --maf ../data/maf/IFNE_100way_vertebrates_hg38.maf \
  --out-prefix ../analysis/qc/ifne_promoter_qc \
  --write-aligned-fasta
```

Current exclusion criteria:

- coverage versus hg38 < 0.50
- N fraction > 0.10
- callable-window fraction < 0.30

Long gaps are recorded but do not automatically cause exclusion.

Repeat for IFNA1 and IFNB.

---

### 3. Build the final species list

The final list contains:

```text
QC-passing species
∩ eutherian mammals
− MAF sanity-check FAIL taxa
```

For IFNE:

```bash
python build_final_keep.py \
  --qc-kept ../analysis/qc/ifne_promoter_qc.kept.txt \
  --keep-list ../config/eutherian_species_hg38_100way.txt \
  --sanity-tsv ../analysis/qc/IFNE_maf_sanity.tsv \
  --out ../analysis/qc/ifne_final.kept.txt
```

Repeat for IFNA1 and IFNB using their corresponding QC and sanity-check files.

---

### 4. Extract promoter FASTAs

```bash
python parse_maf.py \
  -i ../data/maf/IFNE_100way_vertebrates_hg38.maf \
  -o ../analysis/fasta/IFNE \
  --keep-species ../analysis/qc/ifne_final.kept.txt
```

Repeat for IFNA1 and IFNB.

`parse_maf.py` deletes and recreates the output directory before writing FASTAs, so old excluded species cannot remain from previous runs.

---

### 5. Run FIMO

Motif files are stored under:

```text
motifs/jaspar/
motifs/hocomoco/
```

Run:

```bash
bash scan_all_filtered.sh \
  ../analysis/fasta/IFNE \
  ../analysis/fimo_out/IFNE
```

Repeat for IFNA1 and IFNB.

Use the same motif set and FIMO settings for all three promoters.

---

### 6. Combine significant FIMO hits

```bash
python combine_fimo_all_hits.py \
  --fimo-dir ../analysis/fimo_out/IFNE \
  --out ../analysis/tables/IFNE_fimo_combined.tsv \
  --config ../config/config.yaml
```

Repeat for IFNA1 and IFNB.

This script:

- uses FIMO's own p- and q-values;
- applies the thresholds defined in `config/config.yaml`;
- retains all significant positional hits;
- removes only exact duplicate calls.

Current main significance threshold:

```text
q <= 0.05
```

---

### 7. Generate the comparative figure

```bash
python make_promoter_figure_final.py \
  --config ../config/config.yaml \
  --ifne-maf ../data/maf/IFNE_100way_vertebrates_hg38.maf \
  --ifne-fimo ../analysis/tables/IFNE_fimo_combined.tsv \
  --ifne-keep-species ../analysis/qc/ifne_final.kept.txt \
  --ifna-maf ../data/maf/IFNA1_100way_vertebrates_hg38.maf \
  --ifna-fimo ../analysis/tables/IFNA1_fimo_combined.tsv \
  --ifna-keep-species ../analysis/qc/ifna_final.kept.txt \
  --ifnb-maf ../data/maf/IFNB_100way_vertebrates_hg38.maf \
  --ifnb-fimo ../analysis/tables/IFNB_fimo_combined.tsv \
  --ifnb-keep-species ../analysis/qc/ifnb_final.kept.txt \
  --tree ../data/trees/hg38_100way.nh \
  --ucsc-json ../data/trees/ucsc_genomes_official.json \
  --out-prefix ../analysis/figures/IFNE_promoter_figure
```

The figure contains:

- **Panel A:** motif prevalence across *IFNE*, *IFNA1* and *IFNB*
- **Panel B:** phylogenetic motif presence/absence matrices
- **Panel C:** positional motif conservation within the *IFNE* orthologous promoter interval

Panels A and B use the species shared across all three final keep-lists.

Panel C uses all species in the IFNE final keep-list.

A retained species with no significant motif hit remains in the analysis and is scored as motif absence.

## Notes

Promoter coordinates are defined relative to the human hg38 interval and compared across homologous Multiz alignment positions. Motif hits represent **predicted cis-regulatory motifs**, not demonstrated transcription-factor binding.

*Last updated: 2026-11-19*
