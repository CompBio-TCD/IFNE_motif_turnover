#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_species_from_maf(maf_path):
    """
    Parse species IDs from MAF 's' lines.
    Example: hg38.chr9 -> hg38
    """
    species = set()
    with open(maf_path) as f:
        for line in f:
            if not line.startswith("s "):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            src = fields[1]
            sp = src.split(".")[0]
            species.add(sp)
    return sorted(species)


def load_alignment_from_maf(maf_path, species_prefixes=None):
    """
    Build one concatenated aligned sequence per species from a MAF.

    - Concatenates all alignment blocks
    - Fills species absent from a block with gaps for that block length

    Returns:
        aln_seqs: dict species -> concatenated aligned sequence
    """
    if species_prefixes is None:
        species_prefixes = get_species_from_maf(maf_path)

    aln_seqs = {sp: [] for sp in species_prefixes}

    for block in AlignIO.parse(maf_path, "maf"):
        present = set()

        # add aligned sequence for species present in this block
        for rec in block:
            sp = rec.id.split(".")[0]
            if sp in aln_seqs:
                aln_seqs[sp].append(str(rec.seq).upper())
                present.add(sp)

        # fill absent species with gaps for this block
        block_len = block.get_alignment_length()
        for sp in aln_seqs:
            if sp not in present:
                aln_seqs[sp].append("-" * block_len)

    # concatenate blocks
    for sp in aln_seqs:
        aln_seqs[sp] = "".join(aln_seqs[sp])

    return aln_seqs


def longest_run(seq, target_chars):
    best = 0
    cur = 0
    for ch in seq:
        if ch in target_chars:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 0
    return best


def ungapped_length(seq):
    return sum(1 for ch in seq if ch != "-")


def ungapped_fraction(seq):
    if len(seq) == 0:
        return 0.0
    return ungapped_length(seq) / len(seq)


def n_fraction(seq):
    """
    Fraction of Ns among non-gap bases.
    """
    ungapped = [ch for ch in seq if ch != "-"]
    if len(ungapped) == 0:
        return 1.0
    return sum(1 for ch in ungapped if ch.upper() == "N") / len(ungapped)


def callable_window_fraction(seq, window_size=100, min_ungapped_frac=0.7, max_n_frac=0.1):
    """
    Sliding-window QC on a sequence string.

    A window is callable if:
      - enough ungapped sequence
      - not too many Ns among non-gap bases

    Returns:
      fraction of windows callable
    """
    if len(seq) < window_size:
        return 0.0

    n_windows = 0
    n_callable = 0

    for i in range(0, len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        ug_frac = ungapped_fraction(window)
        n_frac = n_fraction(window)
        n_windows += 1
        if ug_frac >= min_ungapped_frac and n_frac <= max_n_frac:
            n_callable += 1

    return n_callable / n_windows if n_windows > 0 else 0.0


def evaluate_alignment_sequences_vs_reference(
    aln_seqs,
    ref_species="hg38",
    min_coverage_vs_ref=0.50,
    max_gap_run=150,
    max_n_frac=0.10,
    use_callable_windows=True,
    callable_window_size=100,
    callable_min_ungapped_frac=0.70,
    callable_max_n_frac=0.10,
    min_callable_window_frac=0.30,
):
    """
    Evaluate species relative to the reference-aligned span only.

    This is the key fix for UCSC multiz MAF output:
    we only assess sequence across columns where the reference species
    (typically hg38) has aligned sequence.

    Returns:
      list of per-species QC result dicts
    """
    if ref_species not in aln_seqs:
        raise ValueError(f"Reference species {ref_species} not found in alignment")

    ref_seq = aln_seqs[ref_species]
    ref_mask = [b != "-" for b in ref_seq]
    ref_aligned_len = sum(ref_mask)

    if ref_aligned_len == 0:
        raise ValueError(f"Reference species {ref_species} has zero aligned bases")

    results = []

    for species, seq in aln_seqs.items():
        # restrict to columns where reference has sequence
        seq_on_ref = "".join(ch for ch, keep in zip(seq, ref_mask) if keep)

        coverage_vs_ref = ungapped_length(seq_on_ref) / ref_aligned_len
        longest_gap = longest_run(seq_on_ref, {"-"})
        n_frac = n_fraction(seq_on_ref)

        if use_callable_windows:
            callable_frac = callable_window_fraction(
                seq_on_ref,
                window_size=callable_window_size,
                min_ungapped_frac=callable_min_ungapped_frac,
                max_n_frac=callable_max_n_frac,
            )
        else:
            callable_frac = ""

        reasons = []
        warnings = []

        # never exclude the reference species
        if species != ref_species:
            if coverage_vs_ref < min_coverage_vs_ref:
                reasons.append(f"low_coverage_vs_{ref_species}:{coverage_vs_ref:.3f}<{min_coverage_vs_ref}")
            if longest_gap > max_gap_run:
                warnings.append(f"long_gap_run:{longest_gap}>{max_gap_run}")
                #reasons.append(f"long_gap_run:{longest_gap}>{max_gap_run}")
            if n_frac > max_n_frac:
                reasons.append(f"high_N_fraction:{n_frac:.3f}>{max_n_frac}")
            if use_callable_windows and callable_frac < min_callable_window_frac:
                reasons.append(f"low_callable_window_fraction:{callable_frac:.3f}<{min_callable_window_frac}")

        results.append({
            "species": species,
            "ref_species": ref_species,
            "ref_aligned_length": ref_aligned_len,
            "coverage_vs_ref": round(coverage_vs_ref, 4),
            "longest_gap_run": longest_gap,
            "N_fraction": round(n_frac, 4),
            "callable_window_fraction": round(callable_frac, 4) if callable_frac != "" else "",
            "exclude": len(reasons) > 0,
            "reasons": ";".join(reasons) if reasons else ""
            "warnings" ";".join(warnings) if  warnings else ""
        })

    return results


def write_outputs(results, out_prefix, aln_seqs=None, write_aligned_fasta=False):
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    metrics_path = out_prefix.with_suffix(".metrics.tsv")
    excluded_path = out_prefix.with_suffix(".excluded.tsv")
    kept_path = out_prefix.with_suffix(".kept.txt")

    fieldnames = [
        "species",
        "ref_species",
        "ref_aligned_length",
        "coverage_vs_ref",
        "longest_gap_run",
        "N_fraction",
        "callable_window_fraction",
        "exclude",
        "reasons",
        "warnings"
    ]

    with open(metrics_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in sorted(results, key=lambda x: x["species"]):
            w.writerow(r)

    with open(excluded_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["species", "reasons"], delimiter="\t")
        w.writeheader()
        for r in sorted(results, key=lambda x: x["species"]):
            if r["exclude"]:
                w.writerow({
                    "species": r["species"],
                    "reasons": r["reasons"]
                })

    with open(kept_path, "w") as f:
        for r in sorted(results, key=lambda x: x["species"]):
            if not r["exclude"]:
                f.write(r["species"] + "\n")

    print(f"Wrote metrics:  {metrics_path}")
    print(f"Wrote excluded: {excluded_path}")
    print(f"Wrote kept:     {kept_path}")

    if write_aligned_fasta and aln_seqs is not None:
        fasta_path = out_prefix.with_suffix(".aligned.fa")
        records = [
            SeqRecord(Seq(seq), id=sp, description="")
            for sp, seq in sorted(aln_seqs.items())
        ]
        SeqIO.write(records, fasta_path, "fasta")
        print(f"Wrote aligned FASTA: {fasta_path}")


def main():
    ap = argparse.ArgumentParser(
        description="QC species in a UCSC MAF alignment relative to the hg38-aligned span."
    )
    ap.add_argument("--maf", required=True, help="Input MAF file")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--ref-species", default="hg38", help="Reference species for aligned-span QC (default: hg38)")
    ap.add_argument("--min-coverage-vs-ref", type=float, default=0.50)
    ap.add_argument("--max-gap-run", type=int, default=150)
    ap.add_argument("--max-n-frac", type=float, default=0.10)
    ap.add_argument("--no-callable-windows", action="store_true", help="Disable sliding-window QC")
    ap.add_argument("--callable-window-size", type=int, default=100)
    ap.add_argument("--callable-min-ungapped-frac", type=float, default=0.70)
    ap.add_argument("--callable-max-n-frac", type=float, default=0.10)
    ap.add_argument("--min-callable-window-frac", type=float, default=0.30)
    ap.add_argument("--write-aligned-fasta", action="store_true", help="Also write concatenated aligned FASTA")

    args = ap.parse_args()

    species_list = get_species_from_maf(args.maf)
    print(f"Found {len(species_list)} species in MAF")

    aln_seqs = load_alignment_from_maf(args.maf, species_list)

    aln_lengths = {len(v) for v in aln_seqs.values()}
    if len(aln_lengths) != 1:
        raise RuntimeError(f"Not all aligned sequences have the same length: {sorted(aln_lengths)}")
    print(f"Concatenated alignment length: {list(aln_lengths)[0]}")

    results = evaluate_alignment_sequences_vs_reference(
        aln_seqs=aln_seqs,
        ref_species=args.ref_species,
        min_coverage_vs_ref=args.min_coverage_vs_ref,
        max_gap_run=args.max_gap_run,
        max_n_frac=args.max_n_frac,
        use_callable_windows=not args.no_callable_windows,
        callable_window_size=args.callable_window_size,
        callable_min_ungapped_frac=args.callable_min_ungapped_frac,
        callable_max_n_frac=args.callable_max_n_frac,
        min_callable_window_frac=args.min_callable_window_frac,
    )

    kept = sum(1 for r in results if not r["exclude"])
    excl = sum(1 for r in results if r["exclude"])

    print(f"Retained: {kept}")
    print(f"Excluded: {excl}")

    if excl:
        print("\nExcluded species and reasons:")
        for r in sorted(results, key=lambda x: x["species"]):
            if r["exclude"]:
                print(f"  {r['species']}: {r['reasons']}")

    write_outputs(
        results=results,
        out_prefix=args.out_prefix,
        aln_seqs=aln_seqs,
        write_aligned_fasta=args.write_aligned_fasta
    )


if __name__ == "__main__":
    main()
