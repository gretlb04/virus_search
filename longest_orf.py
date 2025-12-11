import sys
from textwrap import wrap

# ---------- FASTA PARSER ----------

def parse_fasta(path):
    name = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].split()[0]  # take first token as ID
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            yield name, "".join(seq_chunks)

# ---------- DNA HELPERS ----------

_complement = str.maketrans("ACGTacgtnN", "TGCAtgcanN")

def revcomp(seq: str) -> str:
    return seq.translate(_complement)[::-1]

# Standard genetic code (nuclear)
CODON_TABLE = {
    # Phenylalanine
    "TTT": "F", "TTC": "F",
    # Leucine
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Methionine (start)
    "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine
    "TAT": "Y", "TAC": "Y",
    # Histidine
    "CAT": "H", "CAC": "H",
    # Glutamine
    "CAA": "Q", "CAG": "Q",
    # Asparagine
    "AAT": "N", "AAC": "N",
    # Lysine
    "AAA": "K", "AAG": "K",
    # Aspartic acid
    "GAT": "D", "GAC": "D",
    # Glutamic acid
    "GAA": "E", "GAG": "E",
    # Cysteine
    "TGT": "C", "TGC": "C",
    # Tryptophan
    "TGG": "W",
    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stops
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
    """Return the longest AA ORF across 6 frames. Split on '*' (stops)."""
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

    with open(out_fa, "w") as out:
        for name, nt_seq in parse_fasta(in_fa):
            orf = longest_orf_aa(nt_seq, min_aa_len=50)  # adjust length cutoff if you like
            if not orf:
                continue
            out.write(f">{name}\n")
            for chunk in wrap(orf, 60):
                out.write(chunk + "\n")

if __name__ == "__main__":
    main()
