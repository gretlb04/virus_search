import sys
import os
from textwrap import wrap

# ---------- FASTA PARSER ----------

def parse_fasta(path):
    """
    Yield (id, desc, seq) for each record.
    id  = first token after '>'
    desc = rest of header line (may be empty)
    """
    name = None
    desc = ""
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(seq_chunks)
                header = line[1:]
                parts = header.split(None, 1)
                name = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            yield name, desc, "".join(seq_chunks)

# ---------- DNA HELPERS ----------

_complement = str.maketrans("ACGTacgtnN", "TGCAtgcanN")

def revcomp(seq: str) -> str:
    return seq.translate(_complement)[::-1]

# Standard genetic code (nuclear)
CODON_TABLE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TAA": "*", "TAG": "*", "TGA": "*",
}

def translate_nt(seq_nt: str) -> str:
    seq_nt = seq_nt.upper().replace("U", "T")
    aa = []
    for i in range(0, len(seq_nt) - 2, 3):
        codon = seq_nt[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)

# ---------- ORF FINDER ----------

def longest_orf_aa(nt_seq: str, min_aa_len: int = 0) -> str:
    """Return the longest AA ORF across 6 frames (no frame metadata, just sequence)."""
    nt_seq = nt_seq.upper().replace("U", "T")
    best_orf = ""

    for strand_seq in (nt_seq, revcomp(nt_seq)):  # + and - strands
        for frame in range(3):
            aa_seq = translate_nt(strand_seq[frame:])
            for orf in aa_seq.split("*"):
                if len(orf) > len(best_orf) and len(orf) >= min_aa_len:
                    best_orf = orf
    return best_orf

# ---------- MAIN ----------

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 longest_orf.py <input.fasta> <output.fasta>")
        sys.exit(1)

    in_fa = sys.argv[1]
    out_fa = sys.argv[2]

    # Run ID from file name, e.g. SRR6823454.rdrp1.mu.fa -> SRR6823454
    base = os.path.basename(in_fa)
    run_id = base.split(".")[0]  # everything before first dot

    with open(out_fa, "w") as out:
        for name, desc, nt_seq in parse_fasta(in_fa):
            orf = longest_orf_aa(nt_seq, min_aa_len=50)  # tweak cutoff if you want
            if not orf:
                continue
            orf_len = len(orf)
            # Unique ID: runID|contigID|longestORF
            new_id = f"{run_id}|{name}|longestORF"
            # Keep original description as a comment-ish tail
            header_desc = f"lenAA={orf_len}"
            if desc:
                header_desc += f" original={desc}"
            out.write(f">{new_id} {header_desc}\n")
            for chunk in wrap(orf, 60):
                out.write(chunk + "\n")

if __name__ == "__main__":
    main()
