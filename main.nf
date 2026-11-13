#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def baktaName() { params.bakta_full.toString().toBoolean() ? 'db' : 'db-light' }

process PREPARE {
    input: path fastas
    output: path 'plasmids'
    script:
    """
    plock_helpers prepare --out plasmids ${fastas}
    """
}

process GGCALLER {
    input: path plasmids
    output: tuple path('ggcaller/GFF'), path('ggcaller/gene_calls.faa')
    script:
    """
    realpath ${plasmids}/*.fasta > refs.txt
    micromamba run -n ggcaller ggcaller --refs refs.txt --kmer ${params.k} --out ggcaller --threads ${task.cpus}
    """
}

process BAKTA_DB {
    publishDir "${params.bakta_cache}", mode: 'copy', overwrite: false
    output: path "${baktaName()}"
    script:
    """
    micromamba run -n bakta curl -fL --retry 50 --retry-all-errors --retry-delay 10 -C - \
        -o ${baktaName()}.tar.xz https://zenodo.org/record/14916843/files/${baktaName()}.tar.xz
    micromamba run -n bakta bakta_db install --db-file ${baktaName()}.tar.xz --output .
    rm ${baktaName()}.tar.xz
    """
}

process BAKTA {
    input:
    tuple path(gff), path(proteins)
    path db
    output: tuple path('annotated_gff'), path('bakta')
    script:
    """
    micromamba run -n bakta bakta_proteins --debug \
    --db "${db}" --output bakta --prefix proteins \
    --threads ${task.cpus} ${proteins}
    plock_helpers annotate --gff ${gff} --bakta bakta/proteins.tsv --out annotated_gff
    """
}

process PANAROO {
    input: tuple path(gff), path(bakta)
    output: tuple path('panaroo'), path(gff)
    script:
    """
    micromamba run -n panaroo panaroo -i ${gff}/*.gff -o panaroo --clean-mode moderate --threads ${task.cpus}
    """
}

process ALIGN {
    input: tuple path(panaroo), path(gff)
    output: path 'msa'
    script:
    """
    micromamba run -n tree plock_helpers align --panaroo ${panaroo} --min-presence ${params.min_presence} \
        --threads ${task.cpus} --out msa
    """
}

process TREE {
    input: path msa
    output:
    path 'clonal_frame_tree.nwk', emit: tree
    path 'clonalframeml.importation_status.txt'
    script:
    """
    micromamba run -n tree iqtree3 -s ${msa}/core_gene_alignment.aln -m GTR+G --fast \
        --seed ${params.seed} -T ${task.cpus} --prefix start
    micromamba run -n tree ClonalFrameML start.treefile ${msa}/core_gene_alignment.aln clonalframeml \
        -num_threads ${task.cpus}
    cp clonalframeml.labelled_tree.newick clonal_frame_tree.nwk
    """
}

process BLOCKS {
    input:
    tuple path(panaroo), path(gff)
    path msa
    path tree
    output: path 'out/*'
    script:
    """
    plock --panaroo ${panaroo} --gff ${gff} --alignments ${msa}/aligned_gene_sequences --tree ${tree} --out out \
        --min-block-size ${params.min_block_size} --min-presence ${params.min_presence} \
        --simulations ${params.simulations} --seed ${params.seed} --threads ${task.cpus}
    """
}

workflow {
    if (!(params.single_job.toString() in ['true', 'false']))
        error "single_job must be true or false"
    if (!(params.no_bakta.toString() in ['true', 'false']))
        error "--no_bakta takes no value"
    if (!(params.bakta_full.toString() in ['true', 'false']))
        error "--bakta_full takes no value"
    if ((params.k as int) < 1 || (params.k as int) % 2 == 0)
        error "k must be a positive odd integer"
    if ((params.simulations as int) < 1)
        error "simulations must be positive"
    fa = Channel.fromPath("${params.plasmids}/*.{fa,fasta,fna}", checkIfExists: true)
        .map { p ->
            if (!(p.name ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*\.(fa|fasta|fna)/))
                error "Invalid plasmid filename: ${p.name}"
            p
        }.collect()
    prepared = PREPARE(fa)
    calls = GGCALLER(prepared)
    if (!params.no_bakta.toString().toBoolean()) {
        database = params.bakta_db ? Channel.value(file(params.bakta_db, checkIfExists: true)) :
            (file("${params.bakta_cache}/${baktaName()}").exists() ?
                Channel.value(file("${params.bakta_cache}/${baktaName()}")) : BAKTA_DB())
        annotations = BAKTA(calls, database)
    } else {
        annotations = calls
    }
    clusters = PANAROO(annotations)
    msa = ALIGN(clusters)
    tree = params.tree ? Channel.value(file(params.tree, checkIfExists: true)) : TREE(msa).tree
    BLOCKS(clusters, msa, tree)
}
