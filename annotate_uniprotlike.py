#!/usr/bin/env python3
"""
annotate_uniprotlike.py

Takes:
  1) Helixer-predicted protein FASTA (query)
  2) DIAMOND tabular output with stitle that contains UniProt-like headers

Produces:
  - FASTA with headers like:
      >sp|<YourQueryID>|<UniProtEntryName> ProteinName OS=<Organism> OX=<TaxID> [GN=Gene] PE=<UniProtPE> SV=<UniProtSV>
    If no confident hit:
      >tr|<YourQueryID>|<YourQueryID> Hypothetical protein OS=... OX=... PE=5 SV=1

  - TSV mapping file (query -> UniProt hit details; keeps isoform accessions in hit_acc)
"""

import argparse
import re
import sys

# -----------------------------
# Helpers
# -----------------------------
def normalize_id(raw_id: str) -> str:
    """Remove leading underscores from IDs."""
    return raw_id.lstrip("_")

def wrap_seq(seq: str, width: int) -> str:
    """Wrap sequence to fixed width (e.g. 60 or 80). Use width=0 to disable."""
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

def read_fasta(path):
    """
    Returns:
      seqs: dict[norm_unique_id] = sequence
      id_map: dict[raw_id] = norm_unique_id

    Notes:
      - Leading '_' stripped from IDs.
      - If two raw IDs collapse to same normalized ID, we rename to <id>__dupN.
    """
    seqs = {}
    id_map = {}
    used = set()

    cur_raw = None
    cur_norm = None
    cur = []

    def finalize():
        nonlocal cur_raw, cur_norm, cur
        if cur_raw is None:
            return
        seqs[cur_norm] = "".join(cur)
        cur_raw = None
        cur_norm = None
        cur = []

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                finalize()
                raw_id = line[1:].split()[0]
                base = normalize_id(raw_id)

                norm_id = base
                if norm_id in used:
                    i = 2
                    while f"{base}__dup{i}" in used:
                        i += 1
                    norm_id = f"{base}__dup{i}"
                    print(
                        f"[warn] ID collision after stripping underscores: "
                        f"{raw_id} -> {base}; renamed to {norm_id}",
                        file=sys.stderr
                    )
                used.add(norm_id)
                id_map[raw_id] = norm_id
                cur_raw = raw_id
                cur_norm = norm_id
            else:
                cur.append(line)

    finalize()
    return seqs, id_map

# UniProt-like header parser for DIAMOND stitle
UNIPROT_RE = re.compile(r"^(?P<db>\w+)\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<desc>.+)$")

def parse_uniprot_title(stitle):
    """
    Parse DIAMOND 'stitle' where the subject database is UniProt FASTA.

    Example:
      sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4

    Returns:
      dict with db, acc (keeps isoform suffix like P12345-2), entry, protein_name, gn, pe, sv, is_isoform
      or None if not parseable.
    """
    m = UNIPROT_RE.match(stitle)
    if not m:
        return None

    desc = m.group("desc")
    protein_name = desc.split(" OS=")[0].strip()

    def grab(tag):
        mm = re.search(rf"\b{tag}=([^\s]+)", desc)
        return mm.group(1) if mm else None

    acc = m.group("acc")  # preserves isoform suffix if present
    is_isoform = ("-" in acc and acc.rsplit("-", 1)[-1].isdigit())

    return {
        "db": m.group("db"),
        "acc": acc,
        "entry": m.group("entry"),
        "protein_name": protein_name,
        "gn": grab("GN"),
        "pe": grab("PE"),
        "sv": grab("SV"),
        "is_isoform": is_isoform,
    }

def load_best_hits(tsv_path, id_map, min_pident, min_qcov, max_evalue):
    """
    Load DIAMOND outfmt 6 with columns:
      qseqid pident length qlen qcovhsp evalue bitscore stitle

    Applies filters and keeps best hit per query by bitscore (tie-breaker lower evalue).
    Returns:
      best[qid] = dict(pident, qcov, evalue, bitscore, stitle)
    """
    best = {}
    with open(tsv_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            raw_qseqid = parts[0]
            qseqid = id_map.get(raw_qseqid, normalize_id(raw_qseqid))

            pident = float(parts[1])
            qcov = float(parts[4])     # percent
            evalue = float(parts[5])
            bitscore = float(parts[6])
            stitle = parts[7]

            if pident < min_pident:
                continue
            if qcov < (min_qcov * 100.0):
                continue
            if evalue > max_evalue:
                continue

            if (qseqid not in best) or (bitscore > best[qseqid]["bitscore"]) or (
                bitscore == best[qseqid]["bitscore"] and evalue < best[qseqid]["evalue"]
            ):
                best[qseqid] = {
                    "pident": pident,
                    "qcov": qcov,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "stitle": stitle,
                }
    return best

# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Query proteins FASTA (Helixer output)")
    ap.add_argument("--hits", required=True, help="DIAMOND hits TSV (outfmt 6 with stitle)")
    ap.add_argument("--out_fasta", required=True, help="Output annotated FASTA")
    ap.add_argument("--out_map", required=True, help="Output mapping TSV")
    ap.add_argument("--organism", required=True, help="Organism name to write into OS=")
    ap.add_argument("--taxid", required=True, help="NCBI taxonomy ID to write into OX=")

    ap.add_argument("--min_pident", type=float, default=30.0)
    ap.add_argument("--min_qcov", type=float, default=0.70)  # 0..1
    ap.add_argument("--max_evalue", type=float, default=1e-5)

    ap.add_argument("--wrap", type=int, default=60,
                    help="FASTA sequence line width (60 or 80 common). Use 0 to disable wrapping.")

    ap.add_argument("--keep_uniprot_accession_tag", action="store_true",
                    help="If set, append 'UniProtKB=<accession>' to the FASTA header (keeps isoform suffix like P12345-2).")

    args = ap.parse_args()

    seqs, id_map = read_fasta(args.fasta)
    best = load_best_hits(args.hits, id_map, args.min_pident, args.min_qcov, args.max_evalue)

    # Reverse map norm->raw for reporting (best effort)
    rev = {v: k for k, v in id_map.items()}

    with open(args.out_fasta, "w") as outfa, open(args.out_map, "w") as outmap:
        outmap.write("\t".join([
            "query_id", "raw_query_id",
            "hit_db", "hit_acc", "hit_entry", "hit_is_isoform",
            "hit_pe", "hit_sv",
            "pident", "qcovhsp", "evalue", "bitscore"
        ]) + "\n")

        for qid, seq in seqs.items():
            raw_qid = rev.get(qid, qid)
            hit = best.get(qid)

            if hit:
                u = parse_uniprot_title(hit["stitle"])
                if u:
                    # -----------------------------
                    # REQUESTED HEADER FORMAT:
                    # >db|<YourQueryID>|<UniProtEntryName> ProteinName OS=... OX=... GN=... PE=... SV=...
                    # -----------------------------
                    db = u["db"]                    # sp or tr from UniProt
                    entry = u["entry"]              # e.g. GPB2_YEAST
                    protein_name = u["protein_name"]

                    gn = u.get("gn")                # GN from UniProt hit (if present)
                    pe = u.get("pe") or "4"         # PE from UniProt hit; fallback
                    sv = u.get("sv") or "1"         # SV from UniProt hit; fallback

                    header = f"{db}|{qid}|{entry} {protein_name} OS={args.organism} OX={args.taxid}"
                    if gn:
                        header += f" GN={gn}"
                    header += f" PE={pe} SV={sv}"

                    if args.keep_uniprot_accession_tag:
                        header += f" UniProtKB={u['acc']}"

                    outfa.write(">" + header + "\n")
                    outfa.write(wrap_seq(seq, args.wrap) + "\n")

                    outmap.write("\t".join([
                        qid, raw_qid,
                        u["db"], u["acc"], u["entry"], str(u["is_isoform"]),
                        str(pe), str(sv),
                        f"{hit['pident']:.2f}", f"{hit['qcov']:.2f}",
                        f"{hit['evalue']:.2e}", f"{hit['bitscore']:.1f}",
                    ]) + "\n")
                    continue

            # No confident hit (or unparseable stitle)
            header = f"tr|{qid}|{qid} Hypothetical protein OS={args.organism} OX={args.taxid} PE=5 SV=1"
            outfa.write(">" + header + "\n")
            outfa.write(wrap_seq(seq, args.wrap) + "\n")
            outmap.write("\t".join([
                qid, raw_qid,
                "", "", "", "False",
                "", "",
                "0", "0", "1", "0"
            ]) + "\n")

if __name__ == "__main__":
    main()
