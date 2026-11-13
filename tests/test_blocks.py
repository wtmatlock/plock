import math
import random

import numpy as np
import pytest

from plock import blocks as plock
from plock.blocks import Tree


def names(n):
    return [f"p{i}" for i in range(n)]


def small_tree():
    # eight plasmids in two clades of four with all branches of length 1
    return Tree("(((p0:1,p1:1):1,(p2:1,p3:1):1):1,((p4:1,p5:1):1,(p6:1,p7:1):1):1);", names(8))


def changes(tree, states, known):
    # fewest changes of each pattern from the reconstruction
    counts = tree.changed(states, known).sum(1)
    return int(counts[0]) if np.ndim(states) == 1 else counts


ALL = np.ones(8, bool)


# a second copy in one plasmid makes that plasmid missing and not the whole gene
def test_multi_copy_calls_count_as_missing():
    cells = ["a1", "a2;a3", "a4", "a5"]
    assert not plock.single_copy("a2;a3") and not plock.single_copy("") and plock.single_copy("a1")
    assert plock.common_enough(cells, 0.75)
    assert not plock.common_enough(cells, 0.9)


def write_panaroo(d, presence, alignments):
    # a panaroo folder with its presence table and core gene alignments
    (d / "aligned_gene_sequences").mkdir(parents=True)
    (d / "gene_presence_absence.csv").write_text("\n".join(presence) + "\n")
    for gene, records in alignments.items():
        (d / "aligned_gene_sequences" / f"{gene}.aln.fas").write_text("".join(f">{h}\n{s}\n" for h, s in records))


# mafft's reversed mark is ignored and a plasmid with two records or a truncated or fused call counts as missing
# while a common cluster without an alignment is counted and skipped
def test_alignments_are_read_with_poor_calls_missing(tmp_path):
    write_panaroo(tmp_path, ["Gene,Non-unique Gene name,Annotation,p0,p1,p2,p3,p4",
                             "g,,,a0,a1,a2,a3,a4", "h,,,b0,b1,b2,b3,b4", "r,,,c0,,,,"],
                  {"g": [("p0;a0", "ACGTA" * 19 + "ACG-A"), ("_R_p1;a1", "acgta" * 20), ("p2;a2", "A" * 100),
                         ("p2;x", "A" * 100), ("p3;a3", "A" * 30 + "-" * 70), ("p4;a4", "C" * 100)]})
    genes = plock.read_panaroo(tmp_path, 0.75)
    assert genes.names == ["g"] and genes.rare == 1 and genes.unaligned == 1
    assert genes.seqs[0][0] == "ACGTA" * 19 + "ACG-A" and genes.seqs[0][1] == "ACGTA" * 20
    assert genes.seqs[0][2] == "" and genes.seqs[0][3] == ""
    assert genes.atypical == 1 and sorted(genes.ids) == [(0, 0), (0, 1), (0, 4)]


# each state in a column is one pattern against the other states with gaps as a state
# and other letters unknown while fixed columns and repeats are dropped
def test_site_patterns_split_each_state_from_the_rest():
    pats = plock.site_patterns(["AAG-", "CAGT", "ACGT", "NCGT"]).tolist()
    assert sorted(pats) == sorted([[1, 0, 1, 2], [0, 1, 0, 2], [1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 0, 0], [0, 1, 1, 1]])
    assert len(plock.site_patterns(["AAAA", "AAAA"])) == 0


# the label is the gene names panaroo lists once each or nothing
def test_cluster_labels_are_the_listed_gene_names():
    assert plock.cluster_label(";mobI") == "mobI"
    assert plock.cluster_label("parA;") == "parA"
    assert plock.cluster_label(";insD2;tnp;insD2") == "insD2/tnp"
    assert plock.cluster_label("") == ""


def test_durations_and_progress_steps():
    assert (plock.span(59), plock.span(61), plock.span(3725)) == ("59s", "1m 01s", "1h 02m")
    assert [plock.tenth(d, 118) for d in (0, 11, 12, 118)] == [0, 0, 1, 10]


# monotone and order preserving with untested hypotheses counted as p of 1
def test_benjamini_hochberg_counts_untested_hypotheses():
    q = plock.benjamini_hochberg([0.04, 0.01, 0.03, 0.5], 4)
    assert q == pytest.approx([0.0533333333, 0.04, 0.0533333333, 0.5])
    assert plock.benjamini_hochberg([0.001, 0.002], 100) == pytest.approx([0.1, 0.1])
    assert plock.benjamini_hochberg([], 10) == []


def random_tree(n, rng):
    # a random binary tree built by merging random lineages
    live = [(f"p{i}", 0.0) for i in range(n)]
    while len(live) > 1:
        a = live.pop(rng.randrange(len(live)))
        b = live.pop(rng.randrange(len(live)))
        h = max(a[1], b[1]) + 0.02 + 0.05 * rng.random()
        live.append((f"({a[0]}:{h - a[1]},{b[0]}:{h - b[1]})", h))
    return live[0][0] + ";"


def evolve(text, expected_changes, rng):
    # a 0 1 pattern from independent flips along branches in proportion to length
    parent, length, name = plock.parse_newick(text)
    total = sum(length[:-1])
    state = [0] * len(parent)
    for i in range(len(parent) - 2, -1, -1):
        state[i] = state[parent[i]] ^ (rng.random() < 1 - math.exp(-expected_changes * length[i] / total))
    tips = {x: state[i] for i, x in enumerate(name) if x is not None}
    return np.array([tips[f"p{i}"] for i in range(len(tips))], np.uint8)


# report json cannot close the script element and round trips
def test_json_is_safe_inside_a_script_element():
    text = plock.embed({"n": "gene\"\\\n</script>"})
    assert "<" not in text
    assert plock.json.loads(text) == {"n": "gene\"\\\n</script>"}


# gaps inside and across the origin and overlapping and with unknown length
def test_circular_gap_takes_the_shorter_way_round():
    assert plock.circular_gap((1, 10), (21, 30), 100) == 10
    assert plock.circular_gap((5, 10), (91, 95), 100) == 9
    assert plock.circular_gap((1, 50), (40, 60), 100) == 0
    assert plock.circular_gap((21, 30), (1, 10), 0) == 10


def test_unknown_options_and_missing_values_are_rejected():
    base = [x for name in ("panaroo", "gff", "tree") for x in (f"--{name}", "x")]
    assert plock.parse_args(base + ["--out", "x", "--seed", "1"]).seed == 1
    for bad in (["--sede", "1"], ["--seed"], ["--out", "--seed", "1"], ["stray"], ["--sim", "5"], base[:2]):
        with pytest.raises(SystemExit):
            plock.parse_args(bad if bad == base[:2] else base + bad)


def test_empty_network_produces_no_blocks():
    assert plock.blocks_of(2, [], 2) == [-1, -1]


# quoted names and internal labels and branch lengths and tips that do not match
def test_newick_names_lengths_and_mismatches():
    ids = ["a b", "c", "d", "e"]
    t = Tree("( ('a b':0.5, c:1e-3)90:2 ,\n (d:3, e:1) );", ids)
    # a root branch without a length is collapsed and the other root branch is then kept
    assert len(t.parent) == 6
    # with both root branches long one of them is collapsed so the split they share is one branch
    assert len(Tree("(('a b':1,c:1):1,(d:1,e:1):1);", ids).parent) == 6
    assert (t.leaf >= 0).sum() == 4
    Tree("('a b':1,c:1,d:1,e:1)", ids)
    for bad, tips in (("(a:1,b:1,c:1,d:1);", ids), ("('a b':1,c:1,d:1);", ids), ("('a b':1,c:1,d:1,e:x);", ids),
                      ("('a b':1,c:1,d:1,d:1);", ids), ("('a b':1,c:1,d:1);", ids[:3]), ("('a b':1,c:1", ids)):
        with pytest.raises(ValueError):
            Tree(bad, tips)


# rooting does not change counts or placements and polytomies count correctly
def test_user_supplied_trees():
    ids = names(6)
    unrooted = Tree("(p0:1,p1:1,((p2:1,p3:1):1,(p4:1,p5:1):1):1);", ids)
    rooted = Tree("((p0:1,p1:1)95:0.5,((p2:1,p3:1)100:1,(p4:1,p5:1)80:1)99:0.5);", ids)
    every = np.ones(6, bool)
    for states in ([0, 0, 1, 1, 1, 1], [0, 1, 0, 1, 0, 1], [1, 1, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1]):
        assert changes(unrooted, states, every) == changes(rooted, states, every)
        assert (unrooted.changed(states, every) == rooted.changed(states, every)).all()
    star = Tree("(p0:1,p1:1,p2:1,p3:1,p4:1,p5:1);", ids)
    assert changes(star, [0, 0, 1, 1, 1, 0], every) == 3
    assert changes(star, [0, 0, 1, 1, 0, 0], every) == 2
    # a tree without lengths keeps its shape
    bare = Tree("((p0,p1),(p2,p3),(p4,p5));", ids)
    assert bare.collapsed == 0 and changes(bare, [1, 1, 0, 0, 0, 0], every) == 1
    with pytest.raises(ValueError):
        Tree("((p0:1,p1:1),(p2:1,p3:1),(p4:1,p5:1,p6:1));", ids)
    with pytest.raises(ValueError):
        Tree("((p0:1,p1:1),(p2:1,p3:1),p4:1);", ids)


# fitch counts with an unknown tip and each change placed on one branch with a pattern and its complement alike
def test_changes_are_counted_and_placed():
    t = small_tree()
    assert changes(t, [0, 0, 0, 0, 1, 1, 1, 1], ALL) == 1
    assert changes(t, [0, 1, 0, 1, 0, 1, 0, 1], ALL) == 4
    assert changes(t, [0] * 8, ALL) == 0
    unknown = ALL.copy()
    unknown[3] = False
    assert changes(t, [0, 0, 0, 1, 1, 1, 1, 1], unknown) == 1
    split = t.changed(np.array([[0, 0, 0, 0, 1, 1, 1, 1], [1, 1, 1, 1, 0, 0, 0, 0]]), ALL)
    assert (split[0] == split[1]).all()
    cherry = t.changed([0, 0, 1, 1, 0, 0, 0, 0], ALL)[0]
    assert t.hi[cherry] - t.lo[cherry] == 2


# branches too short to hold a substitution are merged so one event is one branch
def test_zero_length_branches_are_collapsed():
    ids = names(6)
    t = Tree("((p0:1,p1:1):0,(p2:1,p3:1):1,(p4:1,p5:1):1);", ids)
    assert t.collapsed == 1 and len(t.parent) == 9
    assert changes(t, [1, 1, 0, 0, 0, 0], np.ones(6, bool)) == 2


# a gene can only be seen to change where it is called on both sides of a branch
def test_branches_are_observable_only_with_calls_on_both_sides():
    t = small_tree()
    called = np.array([ALL, [True] * 4 + [False] * 4])
    seen = t.observable(called)
    assert seen[0].all()
    inside = (t.lo >= t.column[4]) & (t.hi <= t.column[7] + 1)
    assert not seen[1][inside].any() and seen[1][~inside].sum() > 0


# items are chosen one at a time in proportion to weight and never twice in a row
def test_draws_are_distinct_and_follow_weights():
    w = np.array([3.0, 1, 1, 1, 1, 1])
    rng = np.random.default_rng(1)
    for k in range(1, 7):
        picked = plock.draw(w, k, 200, rng)
        assert all(len(set(row)) == k for row in picked.tolist())
    first = np.bincount(plock.draw(w, 1, 40000, rng).ravel(), minlength=6)
    assert abs(first[0] / first[1:].mean() - 3) < 0.15
    # the heavy item is among two picks with chance 36 in 56 and among five with chance 55 in 56
    assert abs((plock.draw(w, 2, 40000, rng) == 0).any(1).mean() - 36 / 56) < 0.01
    assert abs((plock.draw(w, 5, 40000, rng) == 0).any(1).mean() - 55 / 56) < 0.004


# sharing every change is rare by chance while sharing none is never surprising
# and a pair stops after ten simulations reach it or runs to the cap
def test_co_change_counts_shared_branches_against_chance():
    w = np.ones(40)
    a = np.zeros(40, bool)
    a[:6] = True
    co, expected, p = plock.co_change(a, a.copy(), w, 2000, 1)
    assert co == 6 and expected < 2 and p == 1 / 2001
    b = np.zeros(40, bool)
    b[20:23] = True
    assert plock.co_change(a, b, w, 2000, 1)[2] == 1.0
    b[:2] = True
    co, expected, p = plock.co_change(a, b, w, 100000, 1)
    assert co == 2 and 0.05 < p < 0.6 and plock.HITS / p == int(plock.HITS / p)
    assert plock.co_change(a, b, w, 100000, 1) == plock.co_change(b, a, w, 100000, 1)
    assert all(type(v) in (int, float) for v in plock.co_change(a, b, w, 100000, 1))


# genes that changed independently stay apart and a pair that changed on the same branches is found
def test_independent_genes_rarely_and_a_shared_history_always_passes():
    rng = random.Random(5)
    text = random_tree(64, rng)
    tree = Tree(text, names(64))
    every = np.ones(64, bool)
    shared = evolve(text, 8.0, rng)
    while changes(tree, shared, every) < 6:
        shared = evolve(text, 8.0, rng)
    genes = [[evolve(text, 4.0, rng) for _ in range(3)] for _ in range(10)]
    genes[0][0] = genes[1][0] = shared
    changed = np.array([tree.changed(np.array(s), every).any(0) for s in genes])
    activity = changed.sum(0)
    on = activity > 0
    pairs = [(a, b) for a in range(10) for b in range(a + 1, 10)]
    p = [plock.co_change(changed[a, on], changed[b, on], activity[on] - changed[a, on] - changed[b, on] + 1, 20000, 1)[2]
         for a, b in pairs]
    q = plock.benjamini_hochberg(p, len(pairs))
    assert q[0] <= 0.05, f"planted pair was not found p {p[0]} q {q[0]}"
    assert sum(x <= 0.05 for x in q[1:]) <= 1


# a pair's simulations do not depend on gene order
def test_pair_seed_ignores_order():
    assert plock.pair_seed(1, "x", "y") == plock.pair_seed(1, "y", "x")
    assert plock.pair_seed(1, "x", "y") != plock.pair_seed(2, "x", "y")
    assert plock.pair_seed(1, "x", "y") != plock.pair_seed(1, "x", "z")


# network aggregation keeps edge and self loop weights
def test_aggregation_counts_each_edge_and_self_loop_once():
    edges = [(0, 0, 3.0), (0, 1, 2.0), (1, 2, 4.0), (2, 2, 1.0)]
    result = plock.aggregate_edges(edges, [0, 0, 1])
    assert result == [(0, 0, 5.0), (0, 1, 4.0), (1, 1, 1.0)]
    assert sum(e[2] for e in edges) == sum(e[2] for e in result)


def balanced(tips):
    # a balanced tree with branches of length 1
    while len(tips) > 1:
        tips = [f"({tips[i]}:1,{tips[i + 1]}:1)" for i in range(0, len(tips), 2)]
    return tips[0] + ";"


def write_inputs(d, n=16):
    # plasmids with six genes of two alleles where genes 0 to 2 and genes 3 to 5 change together
    rng = random.Random(2)

    def dna(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    groups = [[i % 2 for i in range(n)], [i // 2 % 2 for i in range(n)]]
    alleles = []
    for j in range(6):
        s = dna(30)
        alleles.append([s, s[:15] + {"A": "C", "C": "G", "G": "T", "T": "A"}[s[15]] + s[16:]])
    spacers = [dna(20) for _ in range(7)]
    (d / "gff").mkdir()
    presence = ["Gene,Non-unique Gene name,Annotation," + ",".join(names(n))]
    cells, records = [[] for _ in range(6)], {f"g{j}": [] for j in range(6)}
    for i in range(n):
        seq, gff = spacers[0], [f"##gff-version 3", f"##sequence-region p{i} 1 0"]
        for j in range(6):
            allele = alleles[j][groups[j // 3][i]]
            gff.append(f"p{i}\tggcaller\tCDS\t{len(seq) + 1}\t{len(seq) + 30}\t.\t+\t0\tID=p{i}_{j}")
            records[f"g{j}"].append((f"p{i};p{i}_{j}", allele))
            cells[j].append(f"p{i}_{j}")
            seq += allele + spacers[j + 1]
        gff[1] = f"##sequence-region p{i} 1 {len(seq)}"
        (d / "gff" / f"p{i}.gff").write_text("\n".join(gff + ["##FASTA", f">p{i}", seq]) + "\n")
    presence += [f"g{j},{'gene' + str(j) if j else ''},," + ",".join(cells[j]) for j in range(6)]
    write_panaroo(d / "panaroo", presence, records)
    (d / "tree.nwk").write_text(balanced(names(n)))


def run_main(d, out):
    plock.main(["--panaroo", str(d / "panaroo"), "--gff", str(d / "gff"), "--tree", str(d / "tree.nwk"),
                "--out", str(out), "--simulations", "5000", "--min-block-size", "2"])
    return {line.split("\t")[0]: line.split("\t")[1:] for line in (out / "blocks.tsv").read_text().splitlines()[1:]}


# the whole tool on a small made up data set finds the two groups and writes every output
def test_end_to_end_on_a_small_data_set(tmp_path):
    write_inputs(tmp_path)
    blocks = run_main(tmp_path, tmp_path / "out")
    assert blocks["g0"][0] == blocks["g1"][0] == blocks["g2"][0] != "unassigned"
    assert blocks["g3"][0] == blocks["g4"][0] == blocks["g5"][0] != "unassigned"
    assert blocks["g0"][0] != blocks["g3"][0]
    out = tmp_path / "out"
    pairs = (out / "gene_pairs.tsv").read_text().splitlines()
    assert len(pairs) == 16
    summary = dict(line.split("\t") for line in (out / "summary.tsv").read_text().splitlines())
    assert summary["blocks"] == "2" and summary["genes_in_blocks"] == "6" and summary["gene_pairs_tested"] == "15"
    assert summary["gene_pairs_tree_q_le_0.05"] == "6"
    html = (out / "plasmids.html").read_text()
    assert "__PARAMS__" not in html and '"l":"gene1"' in html
    assert (out / "block1.fasta").read_text().count(">") == 16
