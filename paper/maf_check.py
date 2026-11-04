#!/usr/bin/env python3

import argparse
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(
        description="Check whether each species maps to one coherent genomic interval in a MAF."
    )
    p.add_argument("maf", help="Input MAF file")
    return p.parse_args()


def main():
    args = parse_args()

    # species -> list of records
    records = defaultdict(list)

    # duplicate species within individual MAF blocks
    duplicate_in_block = defaultdict(int)

    block_id = 0
    current_species = set()

    with open(args.maf) as fh:
        for line in fh:
            line = line.rstrip()

            # blank line marks end of alignment block
            if not line:
                block_id += 1
                current_species = set()
                continue

            if not line.startswith("s "):
                continue

            fields = line.split()

            # MAF s-line format:
            # s src start size strand srcSize text
            src = fields[1]
            start = int(fields[2])
            size = int(fields[3])
            strand = fields[4]
            src_size = int(fields[5])

            species = src.split(".", 1)[0]

            if species in current_species:
                duplicate_in_block[species] += 1

            current_species.add(species)

            # Convert to forward-strand genomic coordinates.
            #
            # For +:
            #   [start, start + size)
            #
            # For -:
            #   MAF start is relative to reverse complement,
            #   so convert back to forward genomic coordinates.
            if strand == "+":
                fwd_start = start
                fwd_end = start + size
            else:
                fwd_start = src_size - start - size
                fwd_end = src_size - start

            records[species].append({
                "block": block_id,
                "src": src,
                "strand": strand,
                "maf_start": start,
                "size": size,
                "src_size": src_size,
                "fwd_start": fwd_start,
                "fwd_end": fwd_end,
            })

    problems = defaultdict(list)

    print(
        f"{'species':24s} "
        f"{'sources':>7s} "
        f"{'blocks':>7s} "
        f"{'strands':>8s} "
        f"{'span_bp':>10s} "
        f"status"
    )
    print("-" * 90)

    for species in sorted(records):
        recs = records[species]

        sources = sorted(set(r["src"] for r in recs))
        strands = sorted(set(r["strand"] for r in recs))
        blocks = sorted(set(r["block"] for r in recs))

        # 1. More than one genomic source/scaffold
        if len(sources) > 1:
            problems[species].append(
                f"multiple sources: {', '.join(sources)}"
            )

        # 2. Strand inconsistency
        if len(strands) > 1:
            problems[species].append(
                f"mixed strands: {', '.join(strands)}"
            )

        # 3. Duplicate species occurrence inside one MAF block
        if duplicate_in_block.get(species, 0) > 0:
            problems[species].append(
                f"appears multiple times within "
                f"{duplicate_in_block[species]} block(s)"
            )

        # Only do interval continuity checks sensibly if one source
        span = None

        if len(sources) == 1:
            # sort intervals by forward genomic coordinate
            sorted_recs = sorted(
                recs,
                key=lambda r: (r["fwd_start"], r["fwd_end"])
            )

            min_start = min(r["fwd_start"] for r in sorted_recs)
            max_end = max(r["fwd_end"] for r in sorted_recs)

            span = max_end - min_start

            # 4. Check for overlapping source intervals
            for prev, curr in zip(sorted_recs[:-1], sorted_recs[1:]):
                if curr["fwd_start"] < prev["fwd_end"]:
                    problems[species].append(
                        "overlapping genomic intervals: "
                        f"{prev['fwd_start']}-{prev['fwd_end']} and "
                        f"{curr['fwd_start']}-{curr['fwd_end']}"
                    )

            # 5. Check whether block order follows genomic order
            #
            # For + strand, genomic coordinates should generally increase
            # with block order.
            #
            # For - strand, they should generally decrease in the MAF's
            # orientation.
            by_block = sorted(recs, key=lambda r: r["block"])

            if len(by_block) > 1 and len(strands) == 1:
                strand = strands[0]
                starts = [r["fwd_start"] for r in by_block]

                if strand == "+":
                    if any(b < a for a, b in zip(starts[:-1], starts[1:])):
                        problems[species].append(
                            "genomic coordinates are not monotonic "
                            "with MAF block order on + strand"
                        )

                elif strand == "-":
                    if any(b > a for a, b in zip(starts[:-1], starts[1:])):
                        problems[species].append(
                            "genomic coordinates are not monotonic "
                            "with MAF block order on - strand"
                        )

        status = "PASS" if not problems[species] else "CHECK"

        span_str = str(span) if span is not None else "NA"

        print(
            f"{species:24s} "
            f"{len(sources):7d} "
            f"{len(blocks):7d} "
            f"{'/'.join(strands):>8s} "
            f"{span_str:>10s} "
            f"{status}"
        )

    print("\nDetailed checks")
    print("----------------")

    any_problem = False

    for species in sorted(problems):
        if not problems[species]:
            continue

        any_problem = True

        print(f"\n{species}")
        for p in problems[species]:
            print(f"  - {p}")

        print("  intervals:")
        for r in sorted(
            records[species],
            key=lambda x: (x["src"], x["fwd_start"])
        ):
            print(
                f"    block={r['block']:4d} "
                f"src={r['src']:25s} "
                f"strand={r['strand']} "
                f"fwd={r['fwd_start']}-{r['fwd_end']} "
                f"size={r['size']}"
            )

    if not any_problem:
        print("PASS: no obvious multi-locus, strand, overlap, or ordering problems detected.")

    print("\nInterpretation")
    print("--------------")
    print("Multiple MAF blocks are normal.")
    print("The main red flags are:")
    print("  - >1 source chromosome/scaffold per species")
    print("  - mixed strand orientations")
    print("  - overlapping source intervals")
    print("  - non-monotonic genomic ordering")
    print("  - duplicate copies of one species within a single block")


if __name__ == "__main__":
    main()
