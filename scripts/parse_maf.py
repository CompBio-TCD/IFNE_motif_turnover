#!/usr/bin/env python3

import argparse
import os
import shutil
from collections import defaultdict

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a MAF file to one gap-free FASTA per species "
            "(one record per chromosome), optionally retaining only "
            "species listed in a QC keep file."
        )
    )

    parser.add_argument(
        "-i", "--maf",
        required=True,
        help="Input MAF file"
    )

    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help=(
            "Output directory for species FASTA files. "
            "WARNING: if this directory already exists, it will be deleted "
            "and recreated before writing output."
        )
    )

    parser.add_argument(
        "--keep-species",
        default=None,
        help=(
            "Optional text file containing one species ID per line to retain, "
            "e.g. the *.kept.txt file produced by the MAF QC script."
        )
    )

    return parser.parse_args()


def load_keep_species(path):
    """
    Load a species whitelist.

    Expected format:
        hg38
        panTro4
        mm10
        ...

    Blank lines and lines beginning with '#' are ignored.
    """
    if path is None:
        return None

    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Keep-species file not found: {path}"
        )

    with open(path) as fh:
        keep = {
            line.strip()
            for line in fh
            if line.strip() and not line.lstrip().startswith("#")
        }

    if not keep:
        raise ValueError(
            f"Keep-species file contains no species IDs: {path}"
        )

    return keep


def main():
    args = parse_args()

    maf_file = os.path.abspath(args.maf)
    outdir = os.path.abspath(args.outdir)

    if not os.path.isfile(maf_file):
        raise FileNotFoundError(
            f"MAF file not found: {maf_file}"
        )

    keep_species = load_keep_species(args.keep_species)

    print("MAF file:", maf_file)
    print("Output directory:", outdir)

    if keep_species is not None:
        print(
            f"QC filtering enabled: keeping "
            f"{len(keep_species)} species from whitelist"
        )
    else:
        print("QC filtering disabled: retaining all species found in MAF")

    # ------------------------------------------------------------
    # Start with a completely clean output directory
    # ------------------------------------------------------------

    if os.path.exists(outdir):
        print(f"Removing existing output directory: {outdir}")
        shutil.rmtree(outdir)

    os.makedirs(outdir, exist_ok=True)

    # ------------------------------------------------------------
    # Parse MAF
    #
    # Structure:
    #     species -> chromosome -> list of ungapped sequence chunks
    # ------------------------------------------------------------

    species_chrom_seqs = defaultdict(
        lambda: defaultdict(list)
    )

    species_seen = set()
    species_skipped_qc = set()

    for block in AlignIO.parse(maf_file, "maf"):

        for rec in block:

            parts = rec.id.split(".", 1)

            if len(parts) != 2:
                continue

            species, chrom = parts

            species_seen.add(species)

            # Skip taxa that failed QC
            if (
                keep_species is not None
                and species not in keep_species
            ):
                species_skipped_qc.add(species)
                continue

            # Remove alignment gaps before writing FASTA
            seq_str = str(rec.seq).replace("-", "").upper()

            if not seq_str:
                continue

            species_chrom_seqs[species][chrom].append(
                seq_str
            )

    # ------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------

    print(f"\nSpecies found in MAF: {len(species_seen)}")

    if keep_species is not None:
        retained_from_maf = (
            species_seen & keep_species
        )

        missing_from_maf = (
            keep_species - species_seen
        )

        print(
            f"Species retained after QC: "
            f"{len(retained_from_maf)}"
        )

        print(
            f"Species skipped because they failed QC: "
            f"{len(species_skipped_qc)}"
        )

        if missing_from_maf:
            print(
                f"WARNING: {len(missing_from_maf)} species "
                f"were listed in the keep file but were not "
                f"found in the MAF:"
            )

            for species in sorted(missing_from_maf):
                print(f"  {species}")

    # ------------------------------------------------------------
    # Write one FASTA per retained species
    # ------------------------------------------------------------

    n_written = 0

    for species in sorted(species_chrom_seqs):

        chrom_dict = species_chrom_seqs[species]

        out_fa = os.path.join(
            outdir,
            f"{species}.fa"
        )

        print("Writing:", out_fa)

        n_records = 0

        with open(out_fa, "w") as out_handle:

            for chrom in sorted(chrom_dict):

                chunks = chrom_dict[chrom]

                full_seq = "".join(chunks)

                if not full_seq:
                    continue

                rec_id = f"{species}_{chrom}"

                rec = SeqRecord(
                    Seq(full_seq),
                    id=rec_id,
                    description=""
                )

                SeqIO.write(
                    rec,
                    out_handle,
                    "fasta"
                )

                n_records += 1

        print(
            f"  Wrote {n_records} chromosome record(s)"
        )

        n_written += 1

    # ------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------

    print("\nDone.")

    print(
        f"FASTA files written: {n_written}"
    )

    print(
        f"Output directory: {outdir}"
    )


if __name__ == "__main__":
    main()
