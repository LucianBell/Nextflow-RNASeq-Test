// Main workflow file

include { script1 } from "$projectDir/scripts/script1.nf"

workflow {
    process myProcess {
        // Use process defined in script1.nf
        script:
        """
        script1_process
        """
    }
}