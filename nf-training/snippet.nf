ch = Channel.of("1","2","3")
ch.view() 

//ch1 = Channel.of(1, 2, 3)
//ch2 = Channel.of(1)
//
//process SUM {
//    input:
//    val x
//    val y
//
//    output:
//    stdout
//
//    script:
//    """
//    echo \$(($x+$y))
//    """
//}
//
//workflow {
//    SUM(ch1, ch2.first()).view()
//}

Channel
    .fromPath('./data/ggal/**.fq', hidden: true)
    .view()

// This script is interacting with the National Center for Biotechnology Information (NCBI) 
// to collect data from the NCBI Sequence Read Archive (SRA). 
params.ncbi_api_key = '4d3bfe7f5b2a080dc8b08800a72451026108'

params.accession = ['ERR908507', 'ERR908506']

process FASTQC {
    input:
    tuple val(sample_id), path(reads_file)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}

workflow {
    // Here, we use the Channel.fromSRA channel factory to query the NCBI SRA archive.
    reads = Channel.fromSRA(params.accession, apiKey: params.ncbi_api_key)
    FASTQC(reads)
}
