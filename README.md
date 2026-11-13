# plock

plock finds blocks of core plasmid genes that co-vary. It takes in a sample of plasmid sequences, then calls, clusters and aligns the core genes. It builds a plasmid tree from those alignments and places each gene's changes on its branches. Two genes that change on the same branches more often than chance are linked. Linked genes are clustered into blocks.

> Reconciling each gene tree with the plasmid tree (for example with AleRax) is more thorough but scales poorly. plock is a faster alternative for large samples.

> Oversampled groups are handled by the tree. plock counts changes on tree branches, so identical plasmids sit on zero-length branches with no changes, meaning they add nothing to the test.

## What you need

- [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/), with Docker running.
- A directory of plasmid sequences, one FASTA file per plasmid. Plasmids are assumed circular. The filename without its extension is the plasmid ID.

## Setup and run

Build the Docker images from the plock directory:

```sh
docker build -t plock/annotation:0.1 -f docker/annotation/Dockerfile .
docker build -t plock/panaroo:0.1 docker/panaroo
docker build -t plock/tree:0.1 -f docker/tree/Dockerfile .
docker build -t plock/plock:0.1 -f docker/plock/Dockerfile .
```

Run the tests:

```sh
docker run --rm -v "$PWD":/plock -w /plock plock/plock:0.1 python -m pytest
```

Run plock:

```sh
nextflow run main.nf --plasmids my_plasmids/
```

On Apple Silicon, try the experimental profile:

```sh
nextflow run main.nf -profile apple_silicon --plasmids my_plasmids/
```

It runs the annotation image under x86-64 emulation, uses an ARM build of DIAMOND and runs ggCaller on one CPU. Emulation is slow. If a tool fails, use an x86-64 Linux host.

The first run downloads Bakta's light database (about 1 GB) to `.plock-cache/bakta-1.12`. Later runs reuse it.

- `--bakta_full` uses the full database (about 32 GB). It names more genes.
- `--bakta_db /path/to/db` uses an existing Bakta 1.12 database.
- `--no_bakta` skips Bakta. Clusters keep Panaroo's `group_` names.

plock builds its own tree. To use yours, pass `--tree mytree.nwk`. It must be Newick with branch lengths in substitutions per site, and tip names must match the plasmid IDs.

## Results

Results go to `results/` by default. Start with `plasmids.html`.

| File | Contents |
|---|---|
| `plasmids.html` | Interactive plot of blocks on each plasmid |
| `summary.tsv` | Run statistics |
| `blocks.tsv` | Block of each gene |
| `gene_pairs.tsv` | Statistics for each gene pair |
| `block1.fasta`, etc. | Block genes joined per plasmid, unaligned |
| `plock.log` | Log |
| `raw/` | Outputs of each tool |

Some genes have no coordinates, mostly sequences Panaroo refound. They go after the positioned genes in block FASTAs and may be missing from the plot. Gene names come from Panaroo.

## Options

For example:

```sh
nextflow run main.nf --plasmids my_plasmids/ --outdir my_results --cpus 4 --memory '8 GB'
```

| Option | Default | Meaning |
|---|---|---|
| `--outdir` | `results` | Output directory |
| `--cpus` | `8` | CPUs per task |
| `--memory` | `16 GB` | Memory per task |
| `--single_job` | `true` | Run one task at a time. `false` allows concurrent tasks |
| `--no_bakta` | off | Skip Bakta |
| `--bakta_full` | off | Use Bakta's full database |
| `--bakta_db` | downloaded | Existing Bakta database |
| `--bakta_cache` | `.plock-cache/bakta-1.12` | Where to cache the downloaded database |
| `--tree` | built | Newick tree to use instead of building one |
| `--min_presence` | `0.90` | Fraction of plasmids a gene must be called in |
| `--min_block_size` | `5` | Fewest genes in a block |
| `--simulations` | `1000000` | Most redraws per gene pair |
| `--seed` | `1` | Seed for IQ-TREE and redraws |
| `--k` | `31` | ggCaller k-mer length |

## How it works

1. ggCaller predicts genes, Bakta names them and Panaroo clusters them.
2. plock keeps clusters called once in at least `--min_presence` of plasmids. A plasmid counts as missing for a gene if it has two copies, or if its call is more than 10% longer or shorter than the median. Copies cannot be matched between plasmids, and truncated or fused calls would add false gap changes.
3. Each gene is aligned by codon. Its proteins (code 11) are aligned with MAFFT and the DNA is threaded back with PAL2NAL, so gaps cover whole codons. Calls with a partial codon or an internal stop are left out. These are almost all Panaroo refound sequences. If they are most of a gene's calls, its DNA is aligned with MAFFT instead.
4. Each variable column gives one pattern per base or gap, marking which plasmids have it. Other letters are unknown.
5. IQ-TREE 3 builds a tree from the joined gene alignments. ClonalFrameML re-estimates its branch lengths without recombination. Branches too short to hold a substitution are collapsed, as is one branch at a two-way root. The topology averages the genes' histories, so a gene's changes are its departures from that average.
6. Fitch parsimony places each pattern's changes on the tree. A gene changed on a branch if any of its patterns did.
7. For each gene pair, plock counts the branches where both changed. Only branches with both genes called on both sides are used. It then redraws the changes of the gene with fewer changes, without replacement, weighted by one plus the number of other genes that changed on each branch. Branches where many genes change, such as a large import, count for less. `tree_p` is the share of redraws with at least as many shared branches.
8. Pairs with `tree_q` ≤ 0.05 form a network weighted by −log10(*p*). Louvain splits it into blocks of at least `--min_block_size` genes.

## Reading the outputs

### `summary.tsv`

- `distinct_core_gene_haplotypes` is the number of distinct plasmid rows across all gene patterns. Plasmids with the same row carry the same core gene variants.
- `gene_changes` is the number of gene changes on the tree, summed over genes.

### `gene_pairs.tsv`

Rows are sorted from strongest evidence down.

- `changes_a`, `changes_b`: branches where each gene changed, among branches where both can be seen.
- `co_changes`: branches where both changed.
- `expected_co_changes`: mean shared branches in the redraws.
- `tree_p`: chance of at least `co_changes` shared branches if the genes changed independently. Redraws stop once 10 reach `co_changes`, giving 10/*L* for *L* redraws (Besag and Clifford). Otherwise *p* = (*r* + 1)/(*S* + 1), with *r* of *S* = `--simulations` redraws reaching it. Pairs with no shared branches get 1.
- `tree_q`: Benjamini–Hochberg adjusted `tree_p` over all gene pairs.
- `median_distance_bp`: median distance between the genes, the shorter way round the plasmid.

### `blocks.tsv`

- `block`: the gene's block, or `unassigned`.
- `plasmids_called`: plasmids with the gene called once.
- `branches_changed`: branches where the gene changed.

### Reading a block

- A block is a set of genes that change on the same branches more often than chance. Often, it is a stretch of DNA that moves as a unit, so its genes are close together. Significant pairs that are further apart are potentially more interesting.
- Genes with few changes have lower power.
- The test is conservative for large blocks, because their own changes raise the weight of the branches they share.
- Blocks come from one partition with no uncertainty estimate. Small blocks may be unstable.

## Dependencies

Please cite these tools if you use plock.

| Tool | Version | Purpose |
|---|---|---|
| Nextflow | 25.04.2 (tested) | Runs the workflow |
| ggCaller | 1.5.2 | Predicts genes |
| Bakta | 1.12.2 | Names genes |
| Panaroo | 1.8.1 | Clusters genes |
| MAFFT | 7.526 | Aligns proteins |
| PAL2NAL | 14.1 | Threads DNA onto protein alignments |
| IQ-TREE | 3.1.4 | Builds the tree |
| ClonalFrameML | 1.20 | Estimates branch lengths without recombination |
