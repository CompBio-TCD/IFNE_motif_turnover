import argparse
import os
from collections import defaultdict
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert MAF to one gap-free FASTA per species (one record per chromosome)."
    )
    parser.add_argument(
        "-i", "--maf",
        required=True,
        help="Input MAF file"
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory for species FASTA files"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    maf_file = os.path.abspath(args.maf)
    outdir = os.path.abspath(args.outdir)

    print("MAF file:", maf_file)
    print("Output directory:", outdir)

    # Explicitly create output directory
    os.makedirs(outdir, exist_ok=True)

    # species -> chrom -> list of sequence chunks
    species_chrom_seqs = defaultdict(lambda: defaultdict(list))

    for block in AlignIO.parse(maf_file, "maf"):
        for rec in block:
            parts = rec.id.split(".", 1)
            if len(parts) != 2:
                continue

            species, chrom = parts
            seq_str = str(rec.seq).replace("-", "")
            if not seq_str:
                continue

            species_chrom_seqs[species][chrom].append(seq_str)

    # Write one FASTA per species
    for species, chrom_dict in species_chrom_seqs.items():
        out_fa = os.path.join(outdir, f"{species}.fa")

        print("Writing:", out_fa)  # <-- hard diagnostic

        with open(out_fa, "w") as out_handle:
            for chrom, chunks in chrom_dict.items():
                full_seq = "".join(chunks)
                if not full_seq:
                    continue

                rec_id = f"{species}_{chrom}"
                rec = SeqRecord(Seq(full_seq), id=rec_id, description="")
                SeqIO.write(rec, out_handle, "fasta")

        print(f"Wrote {len(chrom_dict)} chromosomes to {out_fa}")


if __name__ == "__main__":
    main()
