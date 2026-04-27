#!/usr/bin/env python3

from __future__ import annotations
import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq


STOP = "*"


@dataclass
class ConservedSite:
    position_aa: int
    expected_aa: str
    label: str


def read_fasta(path: str) -> Dict[str, str]:
    d = {}
    for rec in SeqIO.parse(path, "fasta"):
        d[rec.id] = str(rec.seq).upper()
    if not d:
        raise ValueError(f"No FASTA records found in {path}")
    return d


def read_metadata(path: str) -> List[Dict[str, str]]:
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sequence_id", "species", "locus"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Metadata file missing columns: {missing}")
        return list(reader)


def read_conserved_sites(path: Optional[str]) -> List[ConservedSite]:
    if not path:
        return []
    out = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"position_aa", "expected_aa", "label"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Conserved-sites file missing columns: {missing}")
        for row in reader:
            out.append(
                ConservedSite(
                    position_aa=int(row["position_aa"]),
                    expected_aa=row["expected_aa"].upper(),
                    label=row["label"],
                )
            )
    return out


def read_signalp(path: Optional[str], id_col: str, call_col: str) -> Dict[str, bool]:
    if not path:
        return {}
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {id_col, call_col}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"SignalP file missing columns: {missing}")
        for row in reader:
            raw = row[call_col].strip().lower()
            out[row[id_col].strip()] = raw in {"1", "true", "yes", "y", "sp", "present"}
    return out


def ungap(s: str) -> str:
    return s.replace("-", "").replace(".", "")


def translate_cds(cds: str) -> Tuple[str, bool]:
    """
    Returns translated protein and whether CDS length was frame-consistent (%3==0).
    Trims trailing partial codon if present for translation only.
    """
    cds = ungap(cds)
    in_frame = (len(cds) % 3 == 0)
    if len(cds) % 3 != 0:
        cds = cds[: len(cds) - (len(cds) % 3)]
    if not cds:
        return "", in_frame
    return str(Seq(cds).translate(to_stop=False)), in_frame


def strip_terminal_stop(protein: str) -> str:
    return protein[:-1] if protein.endswith(STOP) else protein


def internal_stop_count(protein: str) -> int:
    return strip_terminal_stop(protein).count(STOP)


def compare_translated_to_protein(translated: str, aligned_protein: str) -> Dict[str, float]:
    """
    Compare translated CDS to ungapped aligned protein.
    """
    prot = ungap(aligned_protein)
    translated = strip_terminal_stop(translated)
    prot = strip_terminal_stop(prot)

    if not prot:
        return {
            "protein_len": 0,
            "translated_len": len(translated),
            "coverage": 0.0,
            "identity": 0.0,
            "length_ratio": 0.0,
        }

    compare_len = min(len(prot), len(translated))
    matches = sum(1 for a, b in zip(prot[:compare_len], translated[:compare_len]) if a == b)
    identity = matches / compare_len if compare_len > 0 else 0.0
    coverage = len(translated) / len(prot) if len(prot) > 0 else 0.0
    length_ratio = min(len(translated), len(prot)) / max(len(translated), len(prot)) if max(len(translated), len(prot)) > 0 else 0.0

    return {
        "protein_len": len(prot),
        "translated_len": len(translated),
        "coverage": coverage,
        "identity": identity,
        "length_ratio": length_ratio,
    }


def aa_at_human_position_from_alignment(query_aln_prot: str, human_aln_prot: str, aa_pos_1based: int) -> Optional[str]:
    """
    Get query amino acid at the column corresponding to a human ungapped aa position.
    """
    target_idx = aa_pos_1based - 1
    human_seen = 0
    for col, aa in enumerate(human_aln_prot):
        if aa not in "-.":
            if human_seen == target_idx:
                qaa = query_aln_prot[col]
                return None if qaa in "-." else qaa
            human_seen += 1
    return None


def classify_sequence(
    seq_id: str,
    translated_protein: str,
    aligned_protein: str,
    human_aligned_protein: str,
    human_aa_len: int,
    conserved_sites: List[ConservedSite],
    signalp_calls: Dict[str, bool],
    cds_in_frame: bool,
) -> Dict[str, object]:
    aa_seq = strip_terminal_stop(translated_protein)
    aa_len = len(aa_seq)
    trunc_frac = aa_len / human_aa_len if human_aa_len > 0 else 0.0
    large_truncation = trunc_frac < 0.8

    has_internal_stop = internal_stop_count(translated_protein) > 0
    signalp_ok = signalp_calls.get(seq_id, None)

    conserved_losses = []
    for site in conserved_sites:
        obs = aa_at_human_position_from_alignment(aligned_protein, human_aligned_protein, site.position_aa)
        if obs != site.expected_aa:
            conserved_losses.append(f"{site.label}:{site.expected_aa}{site.position_aa}->{obs or 'NA'}")

    unresolved = False
    unresolved_reasons = []

    # unresolved if too short to assess meaningfully
    if aa_len < 0.3 * human_aa_len:
        unresolved = True
        unresolved_reasons.append("very_short_sequence")

    # unresolved if many gaps around conserved sites
    missing_conserved = sum(1 for site in conserved_sites
                            if aa_at_human_position_from_alignment(aligned_protein, human_aligned_protein, site.position_aa) is None)
    if conserved_sites and missing_conserved > 0.5 * len(conserved_sites):
        unresolved = True
        unresolved_reasons.append("insufficient_conserved_site_coverage")

    disruptive = any([
        has_internal_stop,
        not cds_in_frame,
        large_truncation,
        len(conserved_losses) > 0,
        signalp_ok is False,
    ])

    if unresolved:
        classification = "unresolved"
    elif disruptive:
        classification = "pseudogene"
    else:
        classification = "intact"

    return {
        "sequence_id": seq_id,
        "aa_length": aa_len,
        "human_aa_length": human_aa_len,
        "length_fraction_of_human": round(trunc_frac, 4),
        "cds_length_divisible_by_3": cds_in_frame,
        "internal_stop": has_internal_stop,
        "large_truncation_lt80pct": large_truncation,
        "lost_conserved_residue": len(conserved_losses) > 0,
        "lost_conserved_residue_details": ";".join(conserved_losses),
        "signal_peptide_present": signalp_ok if signalp_ok is not None else "",
        "classification": classification,
        "unresolved_reasons": ";".join(unresolved_reasons),
    }


def support_score(
    translated_protein: str,
    aligned_protein: str,
    cds_in_frame: bool,
) -> Tuple[float, Dict[str, float]]:
    """
    Score support for choosing one representative per species/locus.
    This is NOT an intactness score. It is a consistency/support score.
    """
    cmp = compare_translated_to_protein(translated_protein, aligned_protein)
    internal_stops = internal_stop_count(translated_protein)
    score = 0.0

    # reward CDS/protein agreement
    score += 100.0 * cmp["identity"]
    score += 25.0 * min(cmp["coverage"], 1.0)
    score += 15.0 * cmp["length_ratio"]

    # penalize poor consistency
    if not cds_in_frame:
        score -= 20.0
    score -= 15.0 * internal_stops

    return score, cmp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-aln", required=True, help="Aligned protein FASTA")
    parser.add_argument("--cds-fasta", required=True, help="CDS FASTA for candidate sequences")
    parser.add_argument("--metadata-tsv", required=True, help="TSV with sequence_id, species, locus")
    parser.add_argument("--human-id", required=True, help="Human IFNE ID in protein alignment")
    parser.add_argument("--conserved-tsv", default=None, help="TSV of conserved residues")
    parser.add_argument("--signalp-tsv", default=None, help="Optional SignalP TSV")
    parser.add_argument("--signalp-id-col", default="sequence_id")
    parser.add_argument("--signalp-col", default="has_signal_peptide")
    parser.add_argument("--selected-output", required=True, help="TSV of representative sequences")
    parser.add_argument("--classification-output", required=True, help="TSV of classifications")
    args = parser.parse_args()

    protein_aln = read_fasta(args.protein_aln)
    cds_fasta = read_fasta(args.cds_fasta)
    metadata = read_metadata(args.metadata_tsv)
    conserved_sites = read_conserved_sites(args.conserved_tsv)
    signalp_calls = read_signalp(args.signalp_tsv, args.signalp_id_col, args.signalp_col)

    if args.human_id not in protein_aln:
        raise ValueError(f"Human ID '{args.human_id}' not found in protein alignment")

    human_protein = ungap(protein_aln[args.human_id]).replace(STOP, "")
    human_aa_len = len(human_protein)

    rows = []
    by_group: Dict[Tuple[str, str], List[Dict[str, object]]] = {}

    for row in metadata:
        seq_id = row["sequence_id"]
        species = row["species"]
        locus = row["locus"]

        if seq_id not in cds_fasta or seq_id not in protein_aln:
            continue

        translated, cds_in_frame = translate_cds(cds_fasta[seq_id])
        score, cmp = support_score(translated, protein_aln[seq_id], cds_in_frame)

        rec = {
            "sequence_id": seq_id,
            "species": species,
            "locus": locus,
            "support_score": round(score, 3),
            "protein_identity_to_translated_cds": round(cmp["identity"], 4),
            "protein_coverage_by_translated_cds": round(cmp["coverage"], 4),
            "protein_length_ratio": round(cmp["length_ratio"], 4),
            "cds_in_frame": cds_in_frame,
            "internal_stop_count": internal_stop_count(translated),
        }
        rows.append(rec)
        by_group.setdefault((species, locus), []).append(rec)

    # choose best-supported representative per species/locus
    selected_rows = []
    selected_ids = set()
    for key, candidates in by_group.items():
        candidates_sorted = sorted(
            candidates,
            key=lambda x: (
                x["support_score"],
                x["protein_identity_to_translated_cds"],
                x["protein_coverage_by_translated_cds"],
                -x["internal_stop_count"],
            ),
            reverse=True,
        )
        best = candidates_sorted[0]
        best["selected"] = True
        selected_rows.append(best)
        selected_ids.add(best["sequence_id"])

    # write selected representatives
    with open(args.selected_output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sequence_id", "species", "locus", "support_score",
                "protein_identity_to_translated_cds",
                "protein_coverage_by_translated_cds",
                "protein_length_ratio", "cds_in_frame",
                "internal_stop_count", "selected"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(selected_rows)

    # classify selected representatives
    class_rows = []
    for row in selected_rows:
        seq_id = row["sequence_id"]
        translated, cds_in_frame = translate_cds(cds_fasta[seq_id])

        class_row = classify_sequence(
            seq_id=seq_id,
            translated_protein=translated,
            aligned_protein=protein_aln[seq_id],
            human_aligned_protein=protein_aln[args.human_id],
            human_aa_len=human_aa_len,
            conserved_sites=conserved_sites,
            signalp_calls=signalp_calls,
            cds_in_frame=cds_in_frame,
        )
        class_row["species"] = row["species"]
        class_row["locus"] = row["locus"]
        class_rows.append(class_row)

    with open(args.classification_output, "w", newline="") as fh:
        fieldnames = [
            "sequence_id", "species", "locus",
            "aa_length", "human_aa_length", "length_fraction_of_human",
            "cds_length_divisible_by_3", "internal_stop",
            "large_truncation_lt80pct", "lost_conserved_residue",
            "lost_conserved_residue_details", "signal_peptide_present",
            "classification", "unresolved_reasons"
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(class_rows)


if __name__ == "__main__":
    main()
