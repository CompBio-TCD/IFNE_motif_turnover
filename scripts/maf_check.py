"""
maf_sanity_check.py

Sanity-check a targeted UCSC/Multiz MAF before concatenating per-species
sequence for motif scanning.

Core checks
-----------
For each species:
  1. Exactly one MAF source chromosome/scaffold.
  2. Consistent strand across all blocks.
  3. Species does not occur twice within a single alignment block.
  4. Source intervals do not overlap.
  5. Source intervals are monotonic with MAF block order.
  6. Report total aligned bp and genomic span.
  7. Flag suspiciously large genomic spans.

Optional exact-reference validation
-----------------------------------
If --reference-map is supplied, the script reconstructs the species sequence
from the MAF (gaps removed, concatenated in MAF block order), extracts the
corresponding genomic segments from the SAME reference assembly using
`samtools faidx`, and checks exact sequence identity.

reference-map format (tab-separated, no header required):
    hg38    /path/to/hg38.fa
    mm10    /path/to/mm10.fa
    bosTau9 /path/to/bosTau9.fa

This exact-reference check is stronger than BLAST for confirming that the
sequence reconstructed from the MAF is exactly what the source assembly
contains. BLAST can still be used independently afterward to assess uniqueness.

Exit codes
----------
0 : all species PASS (warnings allowed only for span threshold)
1 : one or more species CHECK/FAIL
2 : usage / input / external-tool error

Example
-------
python maf_sanity_check.py IFNE_100way.maf

python maf_sanity_check.py IFNE_100way.maf \
    --max-span 5000 \
    --reference-map references.tsv \
    --tsv-out IFNE_MAF_QC.tsv

Notes
-----
- Multiple MAF blocks are normal.
- A minus-strand species is not a problem by itself.
- The dangerous case is concatenating sequence from multiple loci/scaffolds,
  mixed strands, duplicate species entries, or incoherent genomic ordering.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class MafRecord:
    block: int
    species: str
    src: str
    maf_start: int
    size: int
    strand: str
    src_size: int
    text: str
    fwd_start: int
    fwd_end: int


def parse_args():
    p = argparse.ArgumentParser(
        description="Sanity-check a targeted MAF before per-species concatenation."
    )
    p.add_argument("maf", help="Input UCSC/Multiz MAF")
    p.add_argument(
        "--max-span",
        type=int,
        default=5000,
        help="Warn if one species spans more than this many genomic bp "
             "(default: 5000). Use 0 to disable."
    )
    p.add_argument(
        "--reference-map",
        default=None,
        help="Optional TSV mapping species IDs to SAME-ASSEMBLY reference FASTAs."
    )
    p.add_argument(
        "--tsv-out",
        default=None,
        help="Optional output TSV summary."
    )
    p.add_argument(
        "--strict-span",
        action="store_true",
        help="Treat span > --max-span as FAIL instead of WARN."
    )
    p.add_argument(
        "--show-intervals",
        action="store_true",
        help="Print all source intervals for every species, not only problematic ones."
    )
    return p.parse_args()


def revcomp(seq: str) -> str:
    table = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return seq.translate(table)[::-1]


def parse_maf(maf_path: Path) -> Tuple[Dict[str, List[MafRecord]], Dict[str, int]]:
    records: Dict[str, List[MafRecord]] = defaultdict(list)
    dup_in_block: Dict[str, int] = defaultdict(int)

    block_id = -1
    seen_in_block = set()
    in_block = False

    with maf_path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if not line.strip():
                if in_block:
                    in_block = False
                    seen_in_block = set()
                continue

            if line.startswith("a"):
                block_id += 1
                in_block = True
                seen_in_block = set()
                continue

            if not line.startswith("s "):
                continue

            # Some MAFs may omit 'a' before first s-line.
            if block_id < 0:
                block_id = 0
                in_block = True

            fields = line.split()
            if len(fields) < 7:
                raise ValueError(f"Malformed MAF s-line:\n{line}")

            src = fields[1]
            start = int(fields[2])
            size = int(fields[3])
            strand = fields[4]
            src_size = int(fields[5])
            text = fields[6]

            species = src.split(".", 1)[0]

            if species in seen_in_block:
                dup_in_block[species] += 1
            seen_in_block.add(species)

            if strand == "+":
                fwd_start = start
                fwd_end = start + size
            elif strand == "-":
                # MAF start is relative to reverse-complement coordinates.
                fwd_start = src_size - start - size
                fwd_end = src_size - start
            else:
                raise ValueError(f"Unexpected strand '{strand}' in:\n{line}")

            records[species].append(
                MafRecord(
                    block=block_id,
                    species=species,
                    src=src,
                    maf_start=start,
                    size=size,
                    strand=strand,
                    src_size=src_size,
                    text=text,
                    fwd_start=fwd_start,
                    fwd_end=fwd_end,
                )
            )

    if not records:
        raise ValueError("No MAF sequence ('s') lines found.")

    return records, dup_in_block


def load_reference_map(path: Path) -> Dict[str, Path]:
    refs = {}
    with path.open() as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                raise ValueError(
                    f"{path}:{lineno}: expected 'species<TAB>fasta'"
                )
            species, fasta = parts[0], parts[1]
            refs[species] = Path(fasta).expanduser().resolve()
    return refs


def samtools_fetch(fasta: Path, region: str) -> str:
    cmd = ["samtools", "faidx", str(fasta), region]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"samtools faidx failed for {region} in {fasta}:\n{proc.stderr.strip()}"
        )
    lines = proc.stdout.splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def reconstruct_maf_sequence(recs: List[MafRecord]) -> str:
    # Same basic assumption as the existing pipeline:
    # concatenate gap-stripped sequence in MAF block order.
    return "".join(
        r.text.replace("-", "")
        for r in sorted(recs, key=lambda x: x.block)
    ).upper()


def reconstruct_reference_sequence(recs: List[MafRecord], fasta: Path) -> str:
    """
    Reconstruct the sequence represented by the MAF using exact source coordinates
    from the same assembly.

    For + strand: fetch forward sequence.
    For - strand: fetch forward genomic sequence, then reverse-complement each
    MAF block so it matches the orientation stored in the MAF.
    """
    pieces = []
    for r in sorted(recs, key=lambda x: x.block):
        # samtools uses 1-based inclusive coordinates.
        region = f"{r.src}:{r.fwd_start + 1}-{r.fwd_end}"
        seq = samtools_fetch(fasta, region).upper()
        if r.strand == "-":
            seq = revcomp(seq)
        pieces.append(seq)
    return "".join(pieces)


def interval_overlap_count(recs: List[MafRecord]) -> Tuple[int, List[str]]:
    msgs = []
    count = 0
    ordered = sorted(recs, key=lambda r: (r.src, r.fwd_start, r.fwd_end))
    for prev, curr in zip(ordered[:-1], ordered[1:]):
        if curr.src != prev.src:
            continue
        if curr.fwd_start < prev.fwd_end:
            count += 1
            msgs.append(
                f"{prev.src}:{prev.fwd_start}-{prev.fwd_end} overlaps "
                f"{curr.src}:{curr.fwd_start}-{curr.fwd_end}"
            )
    return count, msgs


def monotonicity_problem(recs: List[MafRecord]) -> str | None:
    if len(recs) < 2:
        return None

    strands = {r.strand for r in recs}
    sources = {r.src for r in recs}

    if len(strands) != 1 or len(sources) != 1:
        return None

    ordered = sorted(recs, key=lambda r: r.block)
    starts = [r.fwd_start for r in ordered]
    strand = ordered[0].strand

    if strand == "+":
        bad = any(b < a for a, b in zip(starts[:-1], starts[1:]))
        if bad:
            return "forward genomic coordinates decrease with MAF block order"
    else:
        bad = any(b > a for a, b in zip(starts[:-1], starts[1:]))
        if bad:
            return "forward genomic coordinates increase with MAF block order on '-' strand"

    return None


def exact_reference_check(
    species: str,
    recs: List[MafRecord],
    fasta: Path,
) -> Tuple[str, str]:
    if not fasta.exists():
        return "FAIL", f"reference FASTA not found: {fasta}"

    maf_seq = reconstruct_maf_sequence(recs)

    try:
        ref_seq = reconstruct_reference_sequence(recs, fasta)
    except Exception as e:
        return "FAIL", str(e)

    if len(maf_seq) != len(ref_seq):
        return (
            "FAIL",
            f"length mismatch: MAF reconstructed={len(maf_seq)}, "
            f"reference reconstructed={len(ref_seq)}",
        )

    mismatches = sum(a != b for a, b in zip(maf_seq, ref_seq))
    if mismatches:
        pct = 100.0 * (len(maf_seq) - mismatches) / len(maf_seq)
        return (
            "FAIL",
            f"{mismatches} mismatch(es); exact identity={pct:.6f}%",
        )

    return "PASS", f"exact sequence match ({len(maf_seq)} bp; 100% identity)"


def main():
    args = parse_args()
    maf_path = Path(args.maf).expanduser().resolve()

    if not maf_path.exists():
        print(f"ERROR: MAF not found: {maf_path}", file=sys.stderr)
        return 2

    try:
        records, dup_in_block = parse_maf(maf_path)
    except Exception as e:
        print(f"ERROR parsing MAF: {e}", file=sys.stderr)
        return 2

    refs = {}
    if args.reference_map:
        refmap_path = Path(args.reference_map).expanduser().resolve()
        try:
            refs = load_reference_map(refmap_path)
        except Exception as e:
            print(f"ERROR reading reference map: {e}", file=sys.stderr)
            return 2

        if shutil.which("samtools") is None:
            print(
                "ERROR: --reference-map was supplied but 'samtools' is not in PATH.",
                file=sys.stderr,
            )
            return 2

    rows = []
    detailed = {}
    any_fail = False

    header = (
        f"{'species':24s} {'sources':>7s} {'blocks':>7s} {'strand':>8s} "
        f"{'aligned_bp':>10s} {'span_bp':>10s} {'ref':>6s}  status"
    )
    print(header)
    print("-" * len(header))

    for species in sorted(records):
        recs = records[species]
        sources = sorted({r.src for r in recs})
        strands = sorted({r.strand for r in recs})
        blocks = sorted({r.block for r in recs})
        aligned_bp = sum(r.size for r in recs)

        if len(sources) == 1:
            min_start = min(r.fwd_start for r in recs)
            max_end = max(r.fwd_end for r in recs)
            span = max_end - min_start
        else:
            span = None

        fatal = []
        warnings_ = []

        if len(sources) != 1:
            fatal.append(f"multiple MAF sources: {', '.join(sources)}")

        if len(strands) != 1:
            fatal.append(f"mixed strands: {', '.join(strands)}")

        if dup_in_block.get(species, 0):
            fatal.append(
                f"species appears more than once within "
                f"{dup_in_block[species]} MAF block(s)"
            )

        overlap_n, overlap_msgs = interval_overlap_count(recs)
        if overlap_n:
            fatal.append(f"{overlap_n} overlapping genomic interval pair(s)")

        mono = monotonicity_problem(recs)
        if mono:
            fatal.append(mono)

        if (
            args.max_span
            and span is not None
            and span > args.max_span
        ):
            msg = (
                f"genomic span {span} bp exceeds --max-span {args.max_span} bp"
            )
            if args.strict_span:
                fatal.append(msg)
            else:
                warnings_.append(msg)

        ref_status = "NA"
        ref_message = ""

        if refs:
            if species not in refs:
                ref_status = "MISS"
                warnings_.append("no reference FASTA supplied for species")
            else:
                ref_status, ref_message = exact_reference_check(
                    species, recs, refs[species]
                )
                if ref_status != "PASS":
                    fatal.append(f"reference check: {ref_message}")

        if fatal:
            status = "FAIL"
            any_fail = True
        elif warnings_:
            status = "WARN"
        else:
            status = "PASS"

        strand_text = "/".join(strands)
        span_text = "NA" if span is None else str(span)

        print(
            f"{species:24s} {len(sources):7d} {len(blocks):7d} "
            f"{strand_text:>8s} {aligned_bp:10d} {span_text:>10s} "
            f"{ref_status:>6s}  {status}"
        )

        details = []
        details.extend(f"FAIL: {x}" for x in fatal)
        details.extend(f"WARN: {x}" for x in warnings_)
        details.extend(f"OVERLAP: {x}" for x in overlap_msgs)
        if ref_message and ref_status == "PASS":
            details.append(f"REF: {ref_message}")

        detailed[species] = details

        rows.append({
            "species": species,
            "n_sources": len(sources),
            "sources": ",".join(sources),
            "n_blocks": len(blocks),
            "strands": ",".join(strands),
            "aligned_bp": aligned_bp,
            "genomic_span_bp": "" if span is None else span,
            "duplicate_in_block_count": dup_in_block.get(species, 0),
            "overlap_pair_count": overlap_n,
            "monotonic_order": "FAIL" if mono else "PASS",
            "reference_check": ref_status,
            "reference_message": ref_message,
            "status": status,
            "messages": " | ".join(details),
        })

    print("\nDetailed findings")
    print("-----------------")
    found_details = False
    for species in sorted(detailed):
        if not detailed[species] and not args.show_intervals:
            continue
        found_details = True
        print(f"\n{species}")
        for msg in detailed[species]:
            print(f"  - {msg}")

        if args.show_intervals or any(
            x.startswith(("FAIL:", "WARN:", "OVERLAP:"))
            for x in detailed[species]
        ):
            for r in sorted(records[species], key=lambda x: x.block):
                print(
                    f"    block={r.block:4d}  src={r.src:24s} "
                    f"strand={r.strand}  "
                    f"fwd={r.fwd_start}-{r.fwd_end}  size={r.size}"
                )

    if not found_details:
        print("No problems detected.")

    if args.tsv_out:
        import csv
        out_path = Path(args.tsv_out).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nWrote TSV: {out_path}")

    n_pass = sum(r["status"] == "PASS" for r in rows)
    n_warn = sum(r["status"] == "WARN" for r in rows)
    n_fail = sum(r["status"] == "FAIL" for r in rows)

    print("\nSummary")
    print("-------")
    print(f"Species: {len(rows)}")
    print(f"PASS:    {n_pass}")
    print(f"WARN:    {n_warn}")
    print(f"FAIL:    {n_fail}")

    if refs:
        print(
            "\nReference validation used exact source-coordinate extraction from "
            "the supplied same-assembly FASTAs."
        )

    if any_fail:
        print(
            "\nRESULT: FAIL — do not blindly concatenate all species until flagged "
            "cases are resolved."
        )
        return 1

    print(
        "\nRESULT: PASS/WARN — no fatal multi-locus/strand/order problems detected."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
