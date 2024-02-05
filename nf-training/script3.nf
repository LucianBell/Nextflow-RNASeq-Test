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
    .stripIndent()

// The fromFilePairs channel factory takes a glob pattern as input and returns a channel of tuples.
// Here, we are calling the function fromFilePairs, passing reads as the input, and storing the output at
// the read_pairs_ch (since the output is a CHANNEL OF TUPLES)
/*
                                                        set is used to determine the channel where
                                                        the output will be stored
        inside the () in the fromFilePairs is where we
        put the available  (checkIfExists: true)
*/
Channel.fromFilePairs(params.reads, checkIfExists: true).set{ read_pairs_ch }
read_pairs_ch.view()
