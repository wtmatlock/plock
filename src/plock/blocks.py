import argparse
import csv
import json
import multiprocessing
import os
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

# one numpy thread per process because work is spread over processes
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(variable, "1")
import numpy as np

from plock import __version__

TREE_Q = 0.05
HITS = 10
LENGTH_TOLERANCE = 0.10
ZERO_LENGTH = 1e-6
STAGES = 5


# settings and logging

def die(message):
    sys.exit(f"plock: {message}")


def need(condition, message):
    if not condition:
        die(message)


def number(x):
    # whole numbers without a decimal point
    return str(int(x)) if x == int(x) else repr(float(x))


def parse_args(argv):
    p = argparse.ArgumentParser(prog="plock", allow_abbrev=False)
    for name in ("panaroo", "gff", "tree"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--alignments")
    p.add_argument("--out", default="results")
    p.add_argument("--min-block-size", type=int, default=5)
    p.add_argument("--min-presence", type=float, default=0.90)
    p.add_argument("--simulations", type=int, default=1000000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--threads", type=int, default=1)
    a = p.parse_args(argv)
    need(a.simulations > 0, "simulations must be positive")
    need(a.min_block_size >= 2, "minimum block size must be at least 2")
    need(0 <= a.min_presence <= 1, "min-presence must be in [0,1]")
    need(a.threads > 0, "threads must be positive")
    return a


def clock(s):
    return f"{s // 3600:02d}:{s // 60 % 60:02d}:{s % 60:02d}"


def span(s):
    if s >= 3600:
        return f"{s // 3600}h {s // 60 % 60:02d}m"
    if s >= 60:
        return f"{s // 60}m {s % 60:02d}s"
    return f"{s}s"


def tenth(done, total):
    # the ten percent step reached
    return 10 if total == 0 else done * 10 // total


class Log:
    # run log to stderr and plock.log with stage lines stamped by time since the start
    def __init__(self, path):
        self.start = time.time()
        self.file = open(path, "w")

    def line(self, text):
        print(text, file=sys.stderr, flush=True)
        self.file.write(text + "\n")
        self.file.flush()

    def stamp(self):
        return f"[{clock(int(time.time() - self.start))}]"

    def stage(self, n, title):
        self.line(f"{self.stamp()} {n}/{STAGES} {title}")

    def info(self, message):
        self.line(" " * 11 + message)

    def warn(self, message):
        self.info(f"warning: {message}")


_task = None


def _call(i):
    return _task(i)


def parallel(task, total, threads, log, label, unit):
    # runs task on 0 to total in forked processes and logs every tenth done
    global _task
    _task = task
    pool = multiprocessing.get_context("fork").Pool(min(threads, total)) if threads > 1 and total > 1 else None
    results = pool.imap(_call, range(total)) if pool else map(_call, range(total))
    since, logged, out = time.time(), 0, []
    for done, result in enumerate(results, 1):
        out.append(result)
        if tenth(done, total) > logged:
            logged = tenth(done, total)
            elapsed = time.time() - since
            tail = f"done in {span(int(elapsed))}" if done >= total else f"eta {span(int(elapsed * (total - done) / done))}"
            log.line(f"{log.stamp()} {label} {logged * 10}% {done}/{total} {unit}, {tail}")
    if pool:
        pool.close()
        pool.join()
    return out


# panaroo calls

def fasta(path):
    name, parts = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = (line[1:].split() or [""])[0], []
                need(name, f"{path}: empty FASTA ID")
            elif line:
                need(name is not None, f"{path}: sequence before header")
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


# a plasmid's call is used only if it is one copy so two copies count as missing

def single_copy(cell):
    return cell != "" and ";" not in cell


def common_enough(cells, min_presence):
    return sum(map(single_copy, cells)) / len(cells) >= min_presence


def cluster_label(gene_names):
    # gene names panaroo lists for a cluster once each joined with a slash
    return "/".join(dict.fromkeys(x.strip() for x in gene_names.split(";") if x.strip()))


def site_patterns(rows):
    # one 0 1 pattern per alignment column and state with gaps as a state and other letters unknown
    # coded 0 or 1 when known and 2 when unknown with repeated patterns kept once
    m = np.frombuffer("".join(rows).encode(), np.uint8).reshape(len(rows), -1)
    seen = np.isin(m, np.frombuffer(b"ACGT-", np.uint8))
    out = [np.zeros((0, len(rows)), np.uint8)]
    for c in b"ACGT-":
        x = m == c
        keep = x.any(0) & (seen & ~x).any(0)
        out.append(np.where(seen[:, keep], x[:, keep], 2).T.astype(np.uint8))
    return np.unique(np.concatenate(out), axis=0)


# the plasmid tree

def parse_newick(text):
    # nodes numbered children before parents with the root last and names on tips only
    parent, length, name, open_nodes = [], [], [], []
    pos = 0

    def skip():
        nonlocal pos
        while pos < len(text) and text[pos].isspace():
            pos += 1

    def peek():
        skip()
        return text[pos] if pos < len(text) else None

    def label():
        nonlocal pos
        skip()
        if pos < len(text) and text[pos] == "'":
            end = text.find("'", pos + 1)
            end = len(text) if end < 0 else end
            value, pos = text[pos + 1:end], end + 1
            return value
        start = pos
        while pos < len(text) and text[pos] not in ":,();":
            pos += 1
        return text[start:pos].strip()

    def finish(kids, value):
        nonlocal pos
        size = 0.0
        if peek() == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",();":
                pos += 1
            raw = text[start:pos].strip()
            try:
                size = float(raw)
            except ValueError:
                raise ValueError(f"tree: invalid branch length '{raw}'")
            if not (size >= 0 and np.isfinite(size)):
                raise ValueError(f"tree: invalid branch length '{raw}'")
        node = len(parent)
        parent.append(-1)
        length.append(size)
        name.append(None if kids else value)
        for kid in kids:
            parent[kid] = node
        return node

    while True:
        if peek() == "(":
            pos += 1
            open_nodes.append([])
            continue
        node = finish([], label())
        while open_nodes:
            open_nodes[-1].append(node)
            c = peek()
            if c == ",":
                pos += 1
                break
            if c != ")":
                raise ValueError(f"tree: expected ',' or ')' at character {pos}")
            pos += 1
            kids = open_nodes.pop()
            node = finish(kids, label())
        else:
            break
    if peek() == ";":
        pos += 1
    if peek() is not None:
        raise ValueError(f"tree: unexpected text after the tree at character {pos}")
    return parent, length, name


class Tree:
    # the plasmid tree used only to place each gene's changes on branches
    def __init__(self, text, plasmids, ignore=()):
        parent, length, name = parse_newick(text)
        if len(plasmids) < 4:
            raise ValueError("tree: at least 4 plasmids are needed")
        # branches too short to hold a substitution are collapsed and so is one branch at a two way root
        # so each split of the plasmids is one branch and one event is not spread over arbitrary branches
        n, root = len(parent), len(parent) - 1
        has_lengths = max(length[:-1], default=0) > ZERO_LENGTH
        drop = [name[i] is None and i != root and has_lengths and length[i] <= ZERO_LENGTH for i in range(n)]
        inner = [i for i in range(n) if parent[i] == root and name[i] is None]
        if parent.count(root) == 2 and inner and not any(drop[i] for i in inner):
            drop[inner[0]] = True
        up = list(parent)
        for i in range(n - 2, -1, -1):
            if drop[up[i]]:
                up[i] = up[up[i]]
        keep = [i for i in range(n) if not drop[i]]
        new = {old: j for j, old in enumerate(keep)}
        parent, name = [new[up[i]] if i != root else -1 for i in keep], [name[i] for i in keep]
        self.collapsed = n - len(keep)
        index = {p: i for i, p in enumerate(plasmids)}
        tip, leaf = [-1] * len(plasmids), [-1] * len(parent)
        for node, tip_name in enumerate(name):
            if tip_name is None or tip_name in ignore:
                continue
            if tip_name not in index:
                raise ValueError(f"tree: tip {tip_name} is not a plasmid in the analysis")
            if tip[index[tip_name]] >= 0:
                raise ValueError(f"tree: tip {tip_name} appears twice")
            tip[index[tip_name]] = node
            leaf[node] = index[tip_name]
        if -1 in tip:
            raise ValueError(f"tree: plasmid {plasmids[tip.index(-1)]} is not in the tree")
        self.parent, self.leaf = np.array(parent), np.array(leaf)
        # tips in node order so the tips below each branch form one run
        rank = np.full(len(parent), -1)
        rank[self.leaf >= 0] = np.arange(len(plasmids))
        lo = np.where(rank >= 0, rank, len(plasmids))
        hi = np.where(rank >= 0, rank + 1, 0)
        for i in range(len(parent) - 1):
            lo[parent[i]] = min(lo[parent[i]], lo[i])
            hi[parent[i]] = max(hi[parent[i]], hi[i])
        self.lo, self.hi = lo[:-1], hi[:-1]
        self.column = rank[tip]

    def changed(self, states, known):
        # branches where each 0 1 pattern changes in a fewest change fitch reconstruction
        # unknown tips are ignored and a node keeps its parent's state whenever its fitch set allows it
        states = np.atleast_2d(states)
        known = np.broadcast_to(known, states.shape)
        n, rows = len(self.parent), states.shape[0]
        c0, c1, kids = (np.zeros((n, rows), np.int32) for _ in range(3))
        zero, one = np.zeros((n, rows), bool), np.zeros((n, rows), bool)
        for i in range(n):
            if self.leaf[i] >= 0:
                s, k = states[:, self.leaf[i]], known[:, self.leaf[i]]
                zero[i], one[i] = k & (s == 0), k & (s == 1)
            else:
                m = np.maximum(c0[i], c1[i])
                zero[i], one[i] = (c0[i] == m) & (kids[i] > 0), (c1[i] == m) & (kids[i] > 0)
            if i + 1 < n:
                up = self.parent[i]
                c0[up] += zero[i]
                c1[up] += one[i]
                kids[up] += zero[i] | one[i]
        state = np.zeros((n, rows), bool)
        state[-1] = one[-1] & ~zero[-1]
        change = np.zeros((n - 1, rows), bool)
        for i in range(n - 2, -1, -1):
            inherited = state[self.parent[i]]
            change[i] = (zero[i] | one[i]) & ~np.where(inherited, one[i], zero[i])
            state[i] = inherited ^ change[i]
        return change.T

    def observable(self, called):
        # branches with called plasmids on both sides which are the only places a gene can be seen to change
        order = np.argsort(self.column)
        cs = np.concatenate([np.zeros((len(called), 1), int), np.cumsum(called[:, order], 1)], 1)
        below = cs[:, self.hi] - cs[:, self.lo]
        return (below > 0) & (below < cs[:, -1:])


# scoring gene pairs

def pair_seed(seed, a, b):
    # seed for one gene pair that is the same whatever else is analysed
    x, y = sorted((a, b))
    mask = (1 << 64) - 1
    h = 0xcbf29ce484222325 ^ (seed * 0x9e3779b97f4a7c15 & mask)
    for byte in x.encode() + b"\0" + y.encode():
        h = ((h ^ byte) * 0x100000001b3) & mask
    return h | 1


def draw(weight, k, size, rng):
    # k distinct items per row chosen one at a time in proportion to weight
    # cheap draws with replacement give each row its first distinct items and rows still short
    # are finished with exponential keys over the other items which continues the same process
    cum = np.cumsum(weight)
    width = 2 * k + 8
    draws = np.minimum(np.searchsorted(cum, rng.random((size, width)) * cum[-1], side="right"), len(weight) - 1)
    order = np.argsort(draws, axis=1, kind="stable")
    ranked = np.take_along_axis(draws, order, 1)
    new = np.ones(ranked.shape, bool)
    new[:, 1:] = ranked[:, 1:] != ranked[:, :-1]
    first = np.sort(np.where(new, order, width), axis=1)[:, :k]
    picked = np.take_along_axis(draws, np.minimum(first, width - 1), 1)
    short = np.flatnonzero(first[:, -1] == width)
    if len(short):
        keys = rng.standard_exponential((len(short), len(weight))) / weight
        found = first[short] < width
        keys[np.nonzero(found)[0], picked[short][found]] = -1
        picked[short] = np.argpartition(keys, k - 1, axis=1)[:, :k]
    return picked


def co_change(a, b, weight, simulations, seed):
    # branches where both genes changed and how often as many arise by chance when the gene with fewer
    # changes has them redrawn in proportion to weight which is one more than how many other genes changed there
    # simulations stop once HITS reach the observed count and p is then HITS over the number run
    if a.sum() < b.sum():
        a, b = b, a
    k, observed = int(b.sum()), int((a & b).sum())
    rng = np.random.default_rng(seed)
    chunk = max(1, int(2e7 // (2 * k + 8)))
    done = reached = total = 0
    size = min(chunk, 100)
    while done < simulations and k:
        size = min(size, simulations - done)
        found = a[draw(weight, k, size, rng)].sum(1)
        hits = np.flatnonzero(found >= observed)
        if reached + len(hits) >= HITS:
            stop = int(hits[HITS - reached - 1]) + 1
            return observed, float(total + found[:stop].sum()) / (done + stop), HITS / (done + stop)
        reached, done, total, size = reached + len(hits), done + size, total + found.sum(), min(2 * size, chunk)
    return observed, float(total) / max(done, 1), (reached + 1) / (done + 1)


def benjamini_hochberg(p, tests):
    # adjusted values in input order with untested hypotheses counted as p of 1
    tests = max(tests, len(p))
    order = sorted(range(len(p)), key=lambda i: p[i])
    q, floor = [1.0] * len(p), 1.0
    for rank in range(len(order) - 1, -1, -1):
        floor = min(floor, p[order[rank]] * tests / (rank + 1))
        q[order[rank]] = floor
    return q


# blocks

def louvain(n, edges):
    node, e, size = list(range(n)), list(edges), n
    while True:
        m2 = sum(2 * x[2] for x in e)
        if m2 <= 0:
            break
        adj, k = [[] for _ in range(size)], [0.0] * size
        for u, v, w in e:
            if u == v:
                k[u] += 2 * w
                continue
            adj[u].append((v, w))
            adj[v].append((u, w))
            k[u] += w
            k[v] += w
        com, tot, moved = list(range(size)), k[:], True
        while moved:
            moved = False
            for u in range(size):
                cu, wc = com[u], {}
                for v, w in adj[u]:
                    wc[com[v]] = wc.get(com[v], 0.0) + w
                tot[cu] -= k[u]
                best, top = cu, wc.get(cu, 0.0) - tot[cu] * k[u] / m2
                for c, w in wc.items():
                    gain = w - tot[c] * k[u] / m2
                    if c != cu and (gain > top + 1e-12 or (abs(gain - top) <= 1e-12 and best != cu and c < best)):
                        best, top = c, gain
                com[u] = best
                tot[best] += k[u]
                moved |= best != cu
        # relabel then carry labels down to the original nodes and repeat on the merged graph
        lab = {}
        for u in range(size):
            lab.setdefault(com[u], len(lab))
        node = [lab[com[x]] for x in node]
        if len(lab) == size:
            break
        e = aggregate_edges(e, [lab[com[u]] for u in range(size)])
        size = len(lab)
    return node


def aggregate_edges(edges, labels):
    agg = {}
    for u, v, w in edges:
        key = tuple(sorted((labels[u], labels[v])))
        agg[key] = agg.get(key, 0.0) + w
    return [(a, b, w) for (a, b), w in sorted(agg.items())]


def blocks_of(n, edges, min_block):
    # blocks of at least min_block genes numbered largest first with minus one for unassigned
    lab = louvain(n, edges)
    count = Counter(lab)
    keep = sorted((x for x in count if count[x] >= min_block), key=lambda x: (-count[x], x))
    rank = {x: i + 1 for i, x in enumerate(keep)}
    return [rank.get(x, -1) for x in lab]


def circular_gap(a, b, length):
    # bases between two genes the shorter way round with overlaps at zero and unknown length taken as linear
    first, second = (a, b) if a[0] <= b[0] else (b, a)
    inner = max(0, second[0] - first[1] - 1)
    end = max(first[1], second[1])
    if length < end:
        return inner
    return min(inner, length - end + max(0, first[0] - 1))


def middle(v):
    return float(np.median(v)) if len(v) else None


# inputs

def read_gffs(directory):
    # cds coordinates and each plasmid's length from its gff header or else its embedded sequence
    pos, plen = {}, {}
    if not Path(directory).is_dir():
        return pos, plen
    for path in Path(directory).iterdir():
        embedded = None
        for line in path.read_text().splitlines():
            if embedded is not None:
                embedded += 0 if line.startswith(">") else len(line.strip())
                continue
            if line.startswith("##FASTA"):
                embedded = 0
                continue
            if line.startswith("##sequence-region "):
                q = line.split()[1:]
                if len(q) >= 3 and q[2].isdigit():
                    plen[q[0]] = int(q[2])
                continue
            p = line.split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            ids = [a[3:] for a in p[8].split(";") if a.startswith("ID=")]
            if ids and p[3].isdigit() and p[4].isdigit():
                pos[ids[0]] = (int(p[3]), int(p[4]), p[6][:1] or "+")
        if embedded and path.stem not in plen:
            plen[path.stem] = embedded
    return pos, plen


def embed(value):
    # json that is safe inside a script element
    return json.dumps(value, separators=(",", ":")).replace("<", "\\u003c")


class Pair(NamedTuple):
    # a gene pair with the branches each changed on where both can be seen and the branches they share
    a: int
    b: int
    changes_a: int
    changes_b: int
    co: int
    expected: float
    p: float
    q: float = None


@dataclass
class Genes:
    # the analysed panaroo clusters
    plasmids: list    # plasmid names in panaroo column order
    clusters: list    # every panaroo cluster row
    names: list       # analysed cluster names
    seqs: list        # aligned dna of each analysed gene in each plasmid or empty when not one copy
    ids: dict         # annotation id of each analysed gene and plasmid index
    rare: int = 0     # clusters below the presence cut off
    multi_calls: int = 0
    atypical: int = 0  # calls far from their gene's usual length
    unaligned: int = 0  # common clusters without a panaroo alignment


def read_panaroo(directory, min_presence, alignments=None):
    # analysed clusters with each single copy call's row in its gene's alignment
    with open(f"{directory}/gene_presence_absence.csv", newline="") as f:
        rows = list(csv.reader(f))
    n = len(rows[0]) - 3
    genes = Genes(rows[0][3:], [r for r in rows[1:] if len(r) >= 3 + n], [], [], {})
    column = {p: i for i, p in enumerate(genes.plasmids)}
    for r in genes.clusters:
        cells = r[3:3 + n]
        if not common_enough(cells, min_presence):
            genes.rare += 1
            continue
        path = Path(alignments or Path(directory, "aligned_gene_sequences"), f"{r[0]}.aln.fas")
        if not path.exists():
            genes.unaligned += 1
            continue
        genes.multi_calls += sum(";" in x for x in cells)
        # headers are plasmid then annotation id and mafft marks reversed sequences with _R_
        records = {}
        for header, seq in fasta(path):
            records.setdefault(header.split(";")[0].removeprefix("_R_"), []).append(seq.upper())
        alleles, ids = [""] * n, {}
        for p, seqs in records.items():
            pi = column.get(p)
            if pi is not None and single_copy(cells[pi]) and len(seqs) == 1:
                alleles[pi], ids[pi] = seqs[0], cells[pi]
        # a truncated or fused call is a poor allele and is treated as missing
        length = {pi: len(alleles[pi]) - alleles[pi].count("-") for pi in ids}
        typical = np.median(list(length.values())) if ids else 0
        for pi in [pi for pi in ids if abs(length[pi] - typical) > LENGTH_TOLERANCE * typical]:
            alleles[pi] = ""
            del ids[pi]
            genes.atypical += 1
        genes.ids.update({(len(genes.names), pi): cell for pi, cell in ids.items()})
        genes.names.append(r[0])
        genes.seqs.append(alleles)
    return genes


def keep_plasmids(genes, keep):
    # the same genes over a subset of plasmids
    new, n = {old: i for i, old in enumerate(keep)}, len(genes.plasmids)
    genes.clusters = [r[:3] + [r[3 + i] for i in keep] + r[3 + n:] for r in genes.clusters]
    genes.plasmids = [genes.plasmids[i] for i in keep]
    genes.seqs = [[s[i] for i in keep] for s in genes.seqs]
    genes.ids = {(gi, new[pi]): cell for (gi, pi), cell in genes.ids.items() if pi in new}
    return genes


def median_gap(a, b, genes, pos, lengths):
    # median bases between two genes over the plasmids carrying both
    found = []
    for pi, length in enumerate(lengths):
        x, y = pos.get(genes.ids.get((a, pi))), pos.get(genes.ids.get((b, pi)))
        if x and y:
            found.append(circular_gap(x[:2], y[:2], length))
    return middle(found)


# outputs

def write_gene_pairs(path, pairs, genes, gaps):
    # every pair of analysed genes from the strongest evidence down
    with open(path, "w") as f:
        f.write("gene_a\tgene_b\tchanges_a\tchanges_b\tco_changes\texpected_co_changes\ttree_p\ttree_q\t"
                "median_distance_bp\n")
        for x in sorted(pairs, key=lambda x: (x.p, x.a, x.b)):
            d = gaps[(x.a, x.b)]
            f.write(f"{genes.names[x.a]}\t{genes.names[x.b]}\t{x.changes_a}\t{x.changes_b}\t{x.co}\t"
                    f"{x.expected:.2f}\t{x.p:.3g}\t{x.q:.3g}\t{'NA' if d is None else number(d)}\n")


def write_block_fastas(out, final, genes, pos):
    # block genes joined in each plasmid's own gene order with genes lacking coordinates last
    unplaced = 0
    for b in sorted(set(final[final > 0].tolist())):
        with open(f"{out}/block{b}.fasta", "w") as f:
            for pi, p in enumerate(genes.plasmids):
                order = []
                for gi in np.flatnonzero(final == b):
                    if genes.seqs[gi][pi]:
                        at = pos.get(genes.ids.get((gi, pi)))
                        unplaced += at is None
                        order.append((at is None, at[0] if at else 0, gi))
                order.sort()
                if order:
                    f.write(f">{p} genes={','.join(genes.names[x[2]] for x in order)} block={b}\n"
                            f"{''.join(genes.seqs[x[2]][pi] for x in order).replace('-', '')}\n")
    return unplaced


def write_blocks_table(path, genes, final, called, changes):
    # each analysed gene with its block and how often it changed
    with open(path, "w") as f:
        f.write("gene\tblock\tplasmids_called\tbranches_changed\n")
        for i, name in enumerate(genes.names):
            f.write(f"{name}\t{final[i] if final[i] > 0 else 'unassigned'}\t{called[i].sum()}\t{changes[i]}\n")


def write_report(path, params, genes, final, pairs, pos, plen):
    # the html report with every panaroo cluster and each plasmid's genes and chords between analysed genes
    gene_index = {x: i for i, x in enumerate(genes.names)}
    cluster_of_gene, clusters, plasmid_genes = {}, [], [[] for _ in genes.plasmids]
    for ci, r in enumerate(genes.clusters):
        gi = gene_index.get(r[0])
        if gi is not None:
            cluster_of_gene[gi] = ci
        clusters.append({"n": r[0], "l": cluster_label(r[1]), "b": -2 if gi is None else int(final[gi])})
        for pi, cell in enumerate(r[3:3 + len(genes.plasmids)]):
            for x in cell.split(";"):
                if x in pos:
                    st, en, strand = pos[x]
                    plasmid_genes[pi].append([ci, st, en, -1 if strand == "-" else 1])
    links = [[cluster_of_gene[x.a], cluster_of_gene[x.b], round(float(-np.log10(x.p)), 2), float(f"{x.q:.3g}")]
             for x in pairs if x.p <= TREE_Q and x.a in cluster_of_gene and x.b in cluster_of_gene]
    plasmids = [{"p": p, "len": plen.get(p, 0), "genes": sorted(plasmid_genes[pi], key=lambda x: x[1])}
                for pi, p in enumerate(genes.plasmids)]
    html = (Path(__file__).resolve().with_name("plot.html").read_text()
            .replace("__PARAMS__", embed(params))
            .replace("__CLUSTERS__", embed(clusters)[1:-1])
            .replace("__LINKS__", embed(links)[1:-1])
            .replace("__PLASMIDS__", embed(plasmids)[1:-1]))
    Path(path).write_text(html)


def plural(n, word):
    return f"{n} {word}{'' if n == 1 else 's'}"


# the run

def main(argv=None):
    c = parse_args(sys.argv[1:] if argv is None else argv)
    os.makedirs(c.out, exist_ok=True)
    log = Log(f"{c.out}/plock.log")
    log.line(f"{log.stamp()} plock {__version__}")
    log.info(f"seed {c.seed}, {plural(c.threads, 'thread')}, up to {c.simulations} tree simulations per pair "
             f"stopping once {HITS} reach the observed co-changes")
    log.info(f"gene presence >= {number(c.min_presence)}, blocks of >= {c.min_block_size} genes")

    # genes
    genes = read_panaroo(c.panaroo, c.min_presence, c.alignments)
    # plasmids without an analysed call carry no evidence and are left out
    empty = [p for pi, p in enumerate(genes.plasmids) if genes.seqs and not any(s[pi] for s in genes.seqs)]
    if empty:
        genes = keep_plasmids(genes, [pi for pi, p in enumerate(genes.plasmids) if p not in set(empty)])
    plasmids, n, g = genes.plasmids, len(genes.plasmids), len(genes.names)
    called = np.array([[a != "" for a in s] for s in genes.seqs], bool).reshape(g, n)
    log.stage(1, "Genes")
    log.info(f"{len(genes.clusters)} gene clusters in {n + len(empty)} plasmids")
    if empty:
        log.warn(f"{plural(len(empty), 'plasmid')} with no analysed gene calls left out {' '.join(empty)}")
    log.info(f"{g} analysed; {genes.rare} below {number(c.min_presence * 100)}% single-copy presence")
    if genes.multi_calls:
        log.info(f"{genes.multi_calls} multi-copy calls in analysed genes treated as missing")
    if genes.atypical:
        log.info(f"{genes.atypical} calls more than {number(LENGTH_TOLERANCE * 100)}% from their gene's median length "
                 f"treated as missing")
    if genes.unaligned:
        message = (f"{genes.unaligned} common clusters have no alignment in "
                   f"{c.alignments or c.panaroo + '/aligned_gene_sequences'}; make them with plock_helpers align "
                   f"and the same --min-presence")
        need(g > 0, message)
        log.warn(message)
    if g < c.min_block_size:
        log.warn("too few genes to form a block")
    try:
        tree = Tree(Path(c.tree).read_text(), plasmids, set(empty))
    except (OSError, ValueError) as e:
        die(e)

    pos, plen = read_gffs(c.gff)
    # each gene's variable alignment columns as patterns over all plasmids with uncalled plasmids unknown
    patterns = []
    for gi in range(g):
        rows = np.flatnonzero(called[gi])
        full = np.full((0, n), 2, np.uint8)
        if len(rows):
            local = site_patterns([genes.seqs[gi][pi] for pi in rows])
            full = np.full((len(local), n), 2, np.uint8)
            full[:, rows] = local
        patterns.append(full)
    # plasmids with the same state at every analysed pattern add no new evidence
    every = np.concatenate(patterns)
    haplotypes = len(np.unique(every.T, axis=0)) if len(every) else 1
    log.info(f"{len(every)} variable site patterns; {haplotypes} distinct core-gene haplotypes")

    # where on the tree each gene changed
    log.stage(2, "Changes")
    branches = len(tree.parent) - 1
    log.info(f"tree: {n} tips, {branches} branches after collapsing {tree.collapsed} too short to hold a substitution "
             f"or at a two way root")
    changed = np.array([tree.changed(s == 1, s < 2).any(0) if len(s) else np.zeros(branches, bool)
                        for s in patterns], bool).reshape(g, branches)
    seen = tree.observable(called)
    # a branch where many other genes changed is expected to hold co-changes so it is weighted by that count
    activity = changed.sum(0)
    log.info(f"{changed.sum()} gene changes on {(activity > 0).sum()} branches; "
             f"a median gene changed on {number(np.median(changed.sum(1)))} branches")

    # every pair of genes tested
    log.stage(3, "Testing")
    todo = [(a, b) for a in range(g) for b in range(a + 1, g)]

    def score(i):
        a, b = todo[i]
        on = seen[a] & seen[b] & (activity > 0)
        x, y = changed[a, on], changed[b, on]
        weight = activity[on] - x - y + 1
        r = co_change(x, y, weight, c.simulations, pair_seed(c.seed, genes.names[a], genes.names[b]))
        return Pair(a, b, int(x.sum()), int(y.sum()), *r)

    pairs = parallel(score, len(todo), c.threads, log, "testing", "pairs")
    pairs = [x._replace(q=q) for x, q in zip(pairs, benjamini_hochberg([x.p for x in pairs], len(pairs)))]
    lengths = [plen.get(p, 0) for p in plasmids]
    gaps = {(x.a, x.b): median_gap(x.a, x.b, genes, pos, lengths) for x in pairs}
    beyond = sum(x.q <= TREE_Q for x in pairs)
    log.info(f"{len(pairs)} gene pairs tested; {beyond} with tree q <= {TREE_Q}")

    # blocks are genes joined by pairs that changed together more often than chance
    log.stage(4, "Blocks")
    edges = [(x.a, x.b, float(-np.log10(x.p))) for x in pairs if x.q <= TREE_Q]
    final = np.array(blocks_of(g, edges, c.min_block_size), int)
    sizes = Counter(final[final > 0].tolist())
    kept = sorted(sizes)
    inside = {b: [x for x in pairs if final[x.a] == final[x.b] == b] for b in kept}

    def apart(within):
        return middle([d for d in (gaps[(x.a, x.b)] for x in within) if d is not None])

    def kb(d):
        return "NA" if d is None else f"{d / 1000:.1f} kb"

    listing = "; ".join(f"block {b} has {sizes[b]} genes, {sum(x.q <= TREE_Q for x in inside[b])} of "
                        f"{len(inside[b])} pairs with tree q <= {TREE_Q}, median distance {kb(apart(inside[b]))}"
                        for b in kept)
    block_genes = int((final > 0).sum())
    log.info(f"Louvain: {plural(len(kept), 'block')} of >= {c.min_block_size} genes{': ' if kept else ''}{listing}; "
             f"{g - block_genes} of {g} genes unassigned")

    # outputs
    log.stage(5, "Output")
    write_gene_pairs(f"{c.out}/gene_pairs.tsv", pairs, genes, gaps)
    unplaced = write_block_fastas(c.out, final, genes, pos)
    if unplaced:
        log.warn(f"{unplaced} block-gene calls have no ggCaller coordinates; placed after positioned genes in block FASTAs")
    write_blocks_table(f"{c.out}/blocks.tsv", genes, final, called, changed.sum(1))
    # run statistics written to summary.tsv and shown in the report from the same values
    summary = [("plasmids", n), ("gene_clusters", len(genes.clusters)), ("core_single_copy_gene_clusters_analysed", g),
               ("distinct_core_gene_haplotypes", haplotypes), ("gene_changes", int(changed.sum())),
               ("gene_pairs_tested", len(pairs)), ("gene_pairs_tree_q_le_0.05", beyond),
               ("blocks", len(kept)), ("genes_in_blocks", block_genes)]
    Path(f"{c.out}/summary.tsv").write_text("statistic\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in summary))
    block_stats = [{"b": b, "genes": sizes[b], "pairs": len(inside[b]), "q": sum(x.q <= TREE_Q for x in inside[b]),
                    "co": middle([x.co for x in inside[b]]), "p": middle([x.p for x in inside[b]]),
                    "distance": apart(inside[b])} for b in kept]
    # significant and tested pairs within and between blocks with unassigned genes as one group
    groups = kept + ([-1] if (final < 0).any() else [])
    at = {b: i for i, b in enumerate(groups)}
    between = [[[0, 0] for _ in groups] for _ in groups]
    for x in pairs:
        i, j = sorted((at[final[x.a]], at[final[x.b]]))
        between[j][i][0] += int(x.q <= TREE_Q)
        between[j][i][1] += 1
    params = {"plasmids": n, "analysed": g, "clusters": len(genes.clusters), "haplotypes": haplotypes,
              "changes": int(changed.sum()), "tested": len(pairs), "beyond": beyond, "block_genes": block_genes,
              "blocks": len(kept), "block_stats": block_stats,
              "between": {"groups": [int(b) for b in groups], "cells": between}}
    write_report(f"{c.out}/plasmids.html", params, genes, final, pairs, pos, plen)
    log.info(f"wrote plasmids.html, summary.tsv, gene_pairs.tsv, blocks.tsv, {len(kept)} block "
             f"FASTA{'' if len(kept) == 1 else 's'} and plock.log to {c.out}")
    log.line(f"{log.stamp()} done in {span(int(time.time() - log.start))}")


if __name__ == "__main__":
    main()
