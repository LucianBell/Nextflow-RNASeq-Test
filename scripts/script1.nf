// Defines three workflow input parameters and uses the groovy println command to print one of these to the console.

// Params -> A Nextflow special object used to pass parameters to the workflow

// params.reads -> declaring, with params, the variable "reads"
// $projectDir -> pre defined variable that represents the project directory
// {1, 2} -> notation for set expansion (meaning it expects to get a gut_1.fq and a gut_2.fq file)
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"

// params.transcriptome_file -> declaring, with params, the variable "transcriptome_file"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"

// params.multiqc -> declaring, with params, the variable "multiqc"
params.multiqc = "$projectDir/multiqc"

params.outdir = "results"

// log.info -> groovy function to print multiline string statements
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads: $params.reads
    transcriptome_file: $params.transcriptome_file
    multiqc: $params.multiqc
    outdir: $params.outdir
""".stripIndent(true)
