/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    /*
    Here, cpus is what we call a DIRECTIVE. Directives are used to specify the execution
    requirements of a process. For example, the cpus directive specifies the number of CPUs required to execute the process.
    Directives can be added under the process declaration.
    */
    cpus 2

    // Defining process input as a path variable called transcriptome
    input:
    path transcriptome

    // Defining process output as a path variable called 'salmon_index'
    output:
    path 'salmon_index'

    // Defining process script (what to do with input)
    /*
    The INDEX process (using the salmon tool) creates salmon_index, an indexed transcriptome 
    that is passed as an output to the index_ch channel.
    */
    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

workflow {
    index_ch = INDEX(params.transcriptome_file)
    index_ch.view()
}
