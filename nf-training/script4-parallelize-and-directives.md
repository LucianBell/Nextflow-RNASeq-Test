# Explaining how nextflow uses parallelization

- **Concept ->** Parallelization is the technique of dividing a large computational task into smaller sub-tasks that can be executed concurrently on multiple processors or cores
- Nextflow parallelizes the execution of your workflow simply by providing multiple sets of input data to your script.

# Explaning nextflow's directives
- **Concept ->** Directives are optional settings that affect the execution of the current process.
- They must be entered at the top of the process body, before any other declaration blocks (input, output, etc), and have the following syntax:

    `
    directive with simple value -> name value
    `

    `
    directive with list value -> name arg1, arg2, arg3
    `

    `
    directive with map value -> name key1: val1, key2: val2
    `

    `
    directive with value and options -> name arg, opt1: val1, opt2: val2
    `
- Some directives are generally available to all processes, while **others depend on the executor currently defined**.

## Accelerator directives
- The accelerator directive allows you to request hardware accelerators (e.g. GPUs) for the task execution.

script example:
---
    process foo {
        accelerator 4, type: 'nvidia-tesla-k80'

        script:
        """
        your_gpu_enabled --command --line
        """
    }
---

## publishDir directive
- The publishDir directive allows you to publish the process output files to a specified folder.

- The publishDir directive can be specified more than once in order to publish output files to different target directories based on different rules.


script example:
---
    process foo {
        publishDir '/data/chunks', mode: 'copy', overwrite: false

        output:
        path 'chunk_*'

        '''
        printf 'Hola' | split -b 1 - chunk_
        '''
    }
---

- **Files are copied into the specified directory in an asynchronous manner**, so they may not be immediately available in the publish directory at the end of the process execution. For this reason, **downstream processes should not try to access output files through the publish directory, but through channels.**


