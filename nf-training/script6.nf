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

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

/*
 * QUANTIFICATION process, reciving the binary index created by INDEX and
 * the RNA-Seq read pair fastq files. This being a tuple composed of two elements 
 * (a value: sample_id and a list of paths to the fastq reads: reads). 
*/

process QUANTIFICATION {
    cpus 4
    memory '4 GB'

    tag "Salmon on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

/*
* In this FASQC process we're implementing a quality control step to our workflow.
* (CPU and memory directives are here to speed up the process)
*/

process FASTQC {
    cpus 4
    memory '4 GB'
    
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

// This step collects the outputs from the QUANTIFICATION and FASTQC processes to create a final report using the MultiQC tool.

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'
    
    // It creates the final report in the results folder in the current work directory.
    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    /*
    * IMPORTANT: note the use of the mix and collect operators chained together to
    * gather the outputs of the QUANTIFICATION and FASTQC processes as a single input.
    *
    * We can combine operators to manipulate channels. Here, mix combine 2 channels and
    * the collect operator is used to complete the channel contents as a sigle element. 
    */
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}
