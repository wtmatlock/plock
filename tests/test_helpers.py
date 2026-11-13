import shutil

import pytest

from plock import helpers


def test_attribute_escaping_and_decoding():
    assert helpers.encode("abc;é=1 /") == "abc%3B%C3%A9%3D1 /"
    assert helpers.decode("abc%3B%C3%A9%3D1 /") == "abc;é=1 /"


def test_invalid_options_are_rejected():
    for bad in (["weights"], ["prepare", "--out"], ["annotate", "--gff", "x", "--bak", "y", "--out", "z"]):
        with pytest.raises(SystemExit):
            helpers.parser().parse_args(bad)


# ids are taken from the filename and files with more than one record are refused
def test_prepare_normalizes_ids_and_rejects_multiple_records(tmp_path):
    source = tmp_path / "sample.fa"
    source.write_text(">original description\nAcgT\n")
    helpers.main(["prepare", "--out", str(tmp_path / "out"), str(source)])
    assert (tmp_path / "out" / "sample.fasta").read_text() == ">sample\nAcgT\n"
    source.write_text(">one\nAC\n>two\nGT\n")
    with pytest.raises(SystemExit):
        helpers.main(["prepare", "--out", str(tmp_path / "bad"), str(source)])


# bakta names replace the old ones and the embedded fasta is kept
def test_annotate_replaces_fields_and_preserves_fasta(tmp_path):
    (tmp_path / "gff").mkdir()
    (tmp_path / "gff" / "sample.gff").write_text(
        "##gff-version 3\ns\tg\tCDS\t1\t3\t.\t+\t0\tID=a%3Bb;gene=old;product=old\n##FASTA\n>s\nACG\n")
    (tmp_path / "proteins.tsv").write_text("# comment\nID\tGene\tProduct\na;b\tnew\talpha;beta\n")
    helpers.main(["annotate", "--gff", str(tmp_path / "gff"), "--bakta", str(tmp_path / "proteins.tsv"),
                  "--out", str(tmp_path / "out")])
    assert (tmp_path / "out" / "sample.gff").read_text() == (
        "##gff-version 3\ns\tg\tCDS\t1\t3\t.\t+\t0\tID=a%3Bb;gene=new;product=alpha%3Bbeta\n##FASTA\n>s\nACG\n")


# the final stop is dropped and a partial codon or a stop inside means no protein
def test_translation_uses_whole_codons_without_inner_stops():
    assert helpers.translate("ATGAAATAA") == "MK"
    assert helpers.translate("ATGNNNTGG") == "MXW"
    assert helpers.translate("ATGAAAT") is None
    assert helpers.translate("ATGTAAAAA") is None
    assert helpers.translate("TAA") is None


def write_panaroo(d, dna, presence):
    d.mkdir()
    (d / "gene_data.csv").write_text("annotation_id,dna_sequence\n" + "".join(f"{k},{v}\n" for k, v in dna.items()))
    (d / "gene_presence_absence.csv").write_text("\n".join(presence) + "\n")


# only single copy calls in common clusters are aligned and calls far from the usual length are left out
def test_core_calls_follow_plock_rules(tmp_path):
    gene = "ATG" + "GCT" * 30 + "TAA"
    write_panaroo(tmp_path / "p", {"a0": gene, "a1": gene, "a2": gene[:30] + "TAA", "a3": gene, "b0": gene},
                  ["Gene,Non-unique Gene name,Annotation,p0,p1,p2,p3", "g,,,a0,a1,a2,a3;x", "h,,,b0,,,"])
    plasmids, calls = helpers.core_calls(tmp_path / "p", 0.75)
    assert plasmids == ["p0", "p1", "p2", "p3"] and list(calls) == ["g"]
    assert [c[0] for c in calls["g"]] == ["p0", "p1"]


# a codon alignment keeps indels in whole codons and leaves out the odd call out of frame
# while a gene mostly out of frame is aligned as dna with every call kept
@pytest.mark.skipif(not (shutil.which("mafft") and shutil.which("pal2nal.pl")), reason="mafft and pal2nal are not installed")
def test_alignment_by_codon(tmp_path):
    base = "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTTAA"
    stop = "ATGTAA" + base[6:]
    dna = {"a0": base, "a1": base[:15] + base[18:], "a2": base.replace("AAA", "AAG", 1), "a3": base[:-3] + "TGA",
           "b0": base, "b1": base, "b2": base, "b3": stop, "c0": base, "c1": stop, "c2": stop, "c3": stop}
    write_panaroo(tmp_path / "p", dna, ["Gene,Non-unique Gene name,Annotation,p0,p1,p2,p3",
                                        "g,,,a0,a1,a2,a3", "h,,,b0,b1,b2,b3", "k,,,c0,c1,c2,c3"])
    helpers.main(["align", "--panaroo", str(tmp_path / "p"), "--out", str(tmp_path / "msa")])
    genes = tmp_path / "msa" / "aligned_gene_sequences"
    rows = dict(helpers.fasta(genes / "g.aln.fas"))
    assert sorted(rows) == ["p0;a0", "p1;a1", "p2;a2", "p3;a3"]
    assert len({len(x) for x in rows.values()}) == 1 and rows["p1;a1"].count("-") == 3
    assert rows["p0;a0"] == base[:-3]
    h = dict(helpers.fasta(genes / "h.aln.fas"))
    assert sorted(h) == ["p0;b0", "p1;b1", "p2;b2"]
    k = dict(helpers.fasta(genes / "k.aln.fas"))
    assert sorted(k) == ["p0;c0", "p1;c1", "p2;c2", "p3;c3"] and k["p0;c0"].replace("-", "") == base
    joined = dict(helpers.fasta(tmp_path / "msa" / "core_gene_alignment.aln"))
    assert list(joined) == ["p0", "p1", "p2", "p3"] and len({len(x) for x in joined.values()}) == 1
