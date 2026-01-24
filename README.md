# IFNE promoter motif turnover & hormone-modulation analysis

This repository is a scaffold to quantify **motif turnover** across mammals in the *IFNE* promoter and relate it to **hormone modulation**.
It standardizes FIMO scans (PGR/AR/ER, ETS, Sp1, NF-κB, GATA, IRFs), computes **composite modules**, and derives **gain/loss** across a phylogeny, with figures (heatmaps & lollipop maps).

> You will need the **MEME Suite** (FIMO) installed locally and motif PWMs (JASPAR/HOCOMOCO). The code does **not** download or run FIMO here.

## Quick start

1. Place species promoter FASTAs in `data/fasta/`, one FASTA per species (e.g., `human.fa`, `mouse.fa`). Each FASTA should contain the **−1000..+100 bp** window relative to the species' TSS.
2. Put PWM motif files (JASPAR/HOCOMOCO MEME-format) into `motifs/`. 
3. (Optional) Build species-specific 1st-order **Markov background** models (order=1) from each promoter set and save to `data/background/<species>.bg`. Otherwise, FIMO will estimate background.
4. Edit `config/config.yaml` with your species list, paths, thresholds, and composite spacing rules.
5. Run scans:
   ```bash
   bash scripts/scan_all.sh
   ```
   The
   ```python
   python merge_tf.py
   ```
   And finally for plotting:
   ```python
   python plot_promoter_alignment.py
   ```
7. Open the notebook:
   ```bash
   jupyter lab notebooks/IFNE_motif_turnover.ipynb
   ```
8. The notebook parses FIMO TSVs from `data/fimo_out/`, builds a composite module table and produces publication-ready plots.

## Included motifs (suggested JASPAR IDs)
- **Steroid receptors**: PGR/NR3C3, AR/NR3C4, ESR1/ESR2 (ERE), GR
- **Partners**: SP1, JUN/FOS (AP-1), NFKB1/RELA, GATA2/GATA3
- **Epithelial program**: ETS1/ELK1/ELF3 (ETS family), GRHL/KLF, TEAD (YAP/TAZ)
- Confirm absence of IRF/ISRE

See `motifs/MOTIF_CATALOG.md` for a starter list of IDs.

## Relevant parameters
- Primary threshold: `p <= 1e-3` and BH `q <= 0.05`; sensitivity at `q <= 0.2`.
- Composite spacing bins: `0–10`, `10–30`, `30–60` bp (configurable).

_Last updated: 2025-11-12_
