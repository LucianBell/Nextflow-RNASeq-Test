// Defining params for the project

// $projectDir gives you the directory of your project

// Param reads
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
// Param transcriptome_file
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
// Param multiqc
params.multiqc = "$projectDir/multiqc"
// Param outdir
params.outdir = "results"

// Multi-line print
log.info """
            PARAMS
*******************************
Reads: ${params.reads}
transcriptome_file: ${params.transcriptome_file}
multiqc: ${params.multiqc}
outdir: ${params.outdir}
""".stripIndent(true)
