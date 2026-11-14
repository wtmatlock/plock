import argparse
import csv
import re
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from urllib.parse import quote, unquote


def fail(message):
    raise ValueError(message)


def fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = (line[1:].split() or [""])[0]
            if not name:
                fail(f"{path}: empty FASTA ID")
            if any(name == r[0] for r in records):
                fail(f"{path}: duplicate FASTA ID {name}")
            records.append([name, ""])
        elif line.strip():
            if not records:
                fail(f"{path}: sequence before header")
            records[-1][1] += line.strip()
    return records


def prepare(a):
    # one fasta per plasmid renamed to its filename stem
    Path(a.out).mkdir()
    names = set()
    for path in map(Path, a.fastas):
        name = path.stem
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or name in names:
            fail(f"invalid or duplicate plasmid filename stem: {name}")
        names.add(name)
        records = fasta(path)
        if len(records) != 1 or not records[0][1]:
            fail(f"{path}: expected exactly one non-empty FASTA record")
        Path(a.out, f"{name}.fasta").write_text(f">{name}\n{records[0][1]}\n")


def encode(value):
    return quote(value, safe=" :/")


def decode(value):
    return unquote(value)


def annotate(a):
    # bakta gene names and products written into the ggcaller gffs
    lines = [x for x in Path(a.bakta).read_text().splitlines() if not x.startswith("#")]
    header = lines[0].split("\t") if lines else []
    if not {"ID", "Gene", "Product"} <= set(header):
        fail(f"{a.bakta}: expected Bakta ID, Gene and Product columns")
    col = {x: header.index(x) for x in ("ID", "Gene", "Product")}
    annotations = {}
    for line in lines[1:]:
        row = line.split("\t")
        name = row[col["ID"]]
        if not name or name in annotations:
            fail(f"{a.bakta}: missing or duplicate protein ID {name!r}")
        annotations[name] = (row[col["Gene"]], row[col["Product"]])
    files = sorted(Path(a.gff).glob("*.gff"))
    if not files:
        fail(f"{a.gff}: no GFF files")
    Path(a.out).mkdir()
    seen = set()
    for path in files:
        out, in_fasta = [], False
        for line in path.read_text().splitlines():
            in_fasta = in_fasta or line.startswith("##FASTA")
            fields = line.split("\t")
            if in_fasta or line.startswith("#") or not line.strip() or (len(fields) == 9 and fields[2] != "CDS"):
                out.append(line)
                continue
            if len(fields) != 9:
                fail(f"{path}: malformed GFF record")
            attrs = [x for x in fields[8].split(";") if x]
            ids = [decode(x[3:]) for x in attrs if x.startswith("ID=")]
            if len(ids) != 1 or ids[0] in seen or ids[0] not in annotations:
                fail(f"{path}: missing, duplicate or unannotated CDS ID {ids}")
            seen.add(ids[0])
            gene, product = annotations[ids[0]]
            attrs = [x for x in attrs if not x.startswith(("gene=", "product="))]
            attrs += [f"gene={encode(gene)}"] * bool(gene) + [f"product={encode(product)}"] * bool(product)
            fields[8] = ";".join(attrs)
            out.append("\t".join(fields))
        Path(a.out, path.name).write_text("\n".join(out) + "\n")
    if not seen:
        fail("no CDS features found in ggCaller GFFs")
    print(f"Annotated {len(seen)} ggCaller CDS features; {len(annotations) - len(seen)} additional protein records")


# bacterial code 11
BASES = "TCAG"
CODE = dict(zip((x + y + z for x in BASES for y in BASES for z in BASES),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
REFOUND = re.compile(r"([0-9]+_refound_[0-9]+)_pseudo")


def translate(dna):
    # protein without the final stop or none for a partial codon or a stop inside
    if len(dna) % 3:
        return None
    protein = "".join(CODE.get(dna[i:i + 3], "X") for i in range(0, len(dna), 3))
    protein = protein[:-1] if protein.endswith("*") else protein
    return protein if protein and "*" not in protein else None


def gene_dna(dna, gene_id):
    # panaroo pseudogene ids fall back to the refound gene they came from
    if gene_id in dna:
        return dna[gene_id]
    m = REFOUND.fullmatch(gene_id)
    return dna.get(m.group(1)) if m else None


def core_calls(panaroo, min_presence):
    # dna of each single copy call in the clusters plock analyses with calls far from the usual length left out
    from plock.blocks import LENGTH_TOLERANCE, common_enough, single_copy
    csv.field_size_limit(sys.maxsize)
    with open(Path(panaroo, "gene_data.csv"), newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        ci, cd = header.index("annotation_id"), header.index("dna_sequence")
        dna = {r[ci]: r[cd].upper() for r in reader if len(r) > max(ci, cd)}
    with open(Path(panaroo, "gene_presence_absence.csv"), newline="") as f:
        rows = list(csv.reader(f))
    plasmids, calls = rows[0][3:], {}
    for r in rows[1:]:
        cells = r[3:3 + len(plasmids)]
        if len(cells) < len(plasmids) or not common_enough(cells, min_presence):
            continue
        found = []
        for p, cell in zip(plasmids, cells):
            if single_copy(cell):
                seq = gene_dna(dna, cell)
                if not seq:
                    fail(f"Panaroo gene_data.csv lacks DNA for {cell}")
                found.append((p, cell, seq))
        typical = statistics.median(len(x[2]) for x in found)
        calls[r[0]] = [x for x in found if abs(len(x[2]) - typical) <= LENGTH_TOLERANCE * typical]
    return plasmids, calls


def mafft(path, threads, out):
    with open(out, "w") as f:
        subprocess.run(["mafft", "--auto", "--quiet", "--preservecase", "--thread", str(threads), str(path)],
                       stdout=f, check=True)


def align(a):
    # each core gene's proteins aligned by mafft and its dna threaded onto them by pal2nal
    # written per gene for plock and joined in plasmid order for the tree
    plasmids, calls = core_calls(a.panaroo, a.min_presence)
    if not calls:
        fail("no core genes to align")
    out = Path(a.out, "aligned_gene_sequences")
    out.mkdir(parents=True)
    rows, by_codon, by_dna, dropped = {p: [] for p in plasmids}, 0, 0, 0
    with tempfile.TemporaryDirectory() as tmp:
        protein, dna, aligned, codons = (Path(tmp, x) for x in ("protein.fasta", "dna.fasta", "protein.aln", "codons.fasta"))
        for gene, found in calls.items():
            # calls not in frame are mostly sequences panaroo refound in the plasmid and are left out
            # unless they are most of the gene when its dna is aligned directly instead
            # short ids keep mafft and pal2nal from touching the names and the final stop codon is dropped
            kept = [(p, cell, seq[:3 * len(q)], q) for p, cell, seq in found if (q := translate(seq))]
            if len(kept) >= 2 and 2 * len(kept) > len(found):
                dropped += len(found) - len(kept)
                protein.write_text("".join(f">s{i}\n{x[3]}\n" for i, x in enumerate(kept)))
                dna.write_text("".join(f">s{i}\n{x[2]}\n" for i, x in enumerate(kept)))
                mafft(protein, a.threads, aligned)
                with open(codons, "w") as f:
                    run = subprocess.run(["pal2nal.pl", str(aligned), str(dna), "-output", "fasta", "-codontable", "11"],
                                         stdout=f, stderr=subprocess.PIPE, text=True)
                if run.returncode:
                    fail(f"pal2nal failed for {gene}: {run.stderr.strip()}")
                by_codon += 1
            elif len(found) >= 2:
                kept = [(p, cell, seq, None) for p, cell, seq in found]
                dna.write_text("".join(f">s{i}\n{x[2]}\n" for i, x in enumerate(kept)))
                mafft(dna, a.threads, codons)
                by_dna += 1
            else:
                continue
            seqs = {name: seq.upper() for name, seq in fasta(codons)}
            if sorted(seqs) != sorted(f"s{i}" for i in range(len(kept))) or len({len(x) for x in seqs.values()}) != 1:
                fail(f"no alignment row for every sequence of {gene}")
            width = len(seqs["s0"])
            Path(out, f"{gene}.aln.fas").write_text(
                "".join(f">{x[0]};{x[1]}\n{seqs[f's{i}']}\n" for i, x in enumerate(kept)))
            here = {x[0]: seqs[f"s{i}"] for i, x in enumerate(kept)}
            for p in plasmids:
                rows[p].append(here.get(p, "-" * width))
    # plasmids with no core gene call would be all gaps and are left out of the tree
    empty = [p for p, r in rows.items() if not any(set(x.upper()) & set("ACGT") for x in r)]
    Path(a.out, "core_gene_alignment.aln").write_text(
        "".join(f">{p}\n{''.join(r)}\n" for p, r in rows.items() if p not in empty))
    print(f"Aligned {by_codon} core genes by codon with {dropped} calls out of frame left out, "
          f"and {by_dna} by dna because most of their calls are out of frame")
    if empty:
        print(f"Left {len(empty)} plasmids with no core gene calls out of the tree {' '.join(empty)}")


def parser():
    p = argparse.ArgumentParser(prog="plock_helpers", allow_abbrev=False)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("prepare", allow_abbrev=False)
    s.add_argument("--out", required=True)
    s.add_argument("fastas", nargs="+")
    s = sub.add_parser("annotate", allow_abbrev=False)
    for x in ("--gff", "--bakta", "--out"):
        s.add_argument(x, required=True)
    s = sub.add_parser("align", allow_abbrev=False)
    for x in ("--panaroo", "--out"):
        s.add_argument(x, required=True)
    s.add_argument("--min-presence", type=float, default=0.90)
    s.add_argument("--threads", type=int, default=1)
    return p


def main(argv=None):
    a = parser().parse_args(argv)
    try:
        {"prepare": prepare, "annotate": annotate, "align": align}[a.command](a)
    except (ValueError, OSError, subprocess.CalledProcessError) as e:
        sys.exit(f"plock_helpers: {e}")


if __name__ == "__main__":
    main()
