---
title: Cache and resume
description: Foundational Nextflow Training Workshop
---

# Execution cache and resume

The Nextflow caching mechanism works by assigning a unique ID to each task which is used to create a separate execution directory where the tasks are executed and the results stored.

The task unique ID is generated as a 128-bit hash value composing the task input values, files and command string.

The workflow work directory is organized as shown below:

```txt
work/
├── 12
│   └── 1adacb582d2198cd32db0e6f808bce
│       ├── genome.fa -> /data/../genome.fa
│       └── index
│           ├── hash.bin
│           ├── header.json
│           ├── indexing.log
│           ├── quasi_index.log
│           ├── refInfo.json
│           ├── rsd.bin
│           ├── sa.bin
│           ├── txpInfo.bin
│           └── versionInfo.json
├── 19
│   └── 663679d1d87bfeafacf30c1deaf81b
│       ├── ggal_gut
│       │   ├── aux_info
│       │   │   ├── ambig_info.tsv
│       │   │   ├── expected_bias.gz
│       │   │   ├── fld.gz
│       │   │   ├── meta_info.json
│       │   │   ├── observed_bias.gz
│       │   │   └── observed_bias_3p.gz
│       │   ├── cmd_info.json
│       │   ├── libParams
│       │   │   └── flenDist.txt
│       │   ├── lib_format_counts.json
│       │   ├── logs
│       │   │   └── salmon_quant.log
│       │   └── quant.sf
│       ├── ggal_gut_1.fq -> /data/../ggal_gut_1.fq
│       ├── ggal_gut_2.fq -> /data/../ggal_gut_2.fq
│       └── index -> /data/../asciidocs/day2/work/12/1adacb582d2198cd32db0e6f808bce/index
```

!!! info

    You can create these plots using the `tree` function if you have it installed.

### How resume works

The `-resume` command-line option allows the continuation of a workflow execution from the last step that was completed successfully:

```bash
nextflow run <script> -resume
```

In practical terms, the workflow is executed from the beginning. However, before launching the execution of a process, Nextflow uses the task unique ID to check if the work directory already exists and that it contains a valid command exit state with the expected output files.

If this condition is satisfied the task execution is skipped and previously computed results are used as the process results.

The first task for which a new output is computed invalidates all downstream executions in the remaining DAG.

### Work directory

The task work directories are created in the folder `work` in the launching path by default. It is recommended that this is a **scratch** storage area that can be cleaned up once the computation is completed.

!!! note

    It is recommended that you store the final output(s) in a different location using one or more [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) directives.

!!! warning

    Make sure to delete your work directory occasionally, else your machine/environment may be filled with unused files.

A different location for the execution work directory can be specified using the command line option `-w`:

```bash
nextflow run <script> -w /some/scratch/dir
```

!!! warning

    If you delete or move the workflow work directory, it will prevent the use of the resume feature in subsequent runs.

The hash code for input files is computed using:

-   The complete file path
-   The file size
-   The last modified timestamp

Therefore, just **touching** a file will invalidate the related task execution.

### How to organize _in-silico_ experiments

It’s good practice to organize each **experiment** in its own folder. The main experiment input parameters should be specified using a Nextflow config file. This makes it simple to track and replicate an experiment over time.

!!! note

    In the same experiment, the same workflow can be executed multiple times, however, launching two (or more) Nextflow instances in the same directory concurrently should be avoided.

The `nextflow log` command lists the executions run in the current folder:

```console
$ nextflow log
```

```console title="Output"
TIMESTAMP            DURATION  RUN NAME          STATUS  REVISION ID  SESSION ID                            COMMAND
2019-05-06 12:07:32  1.2s      focused_carson    ERR     a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run hello
2019-05-06 12:08:33  21.1s     mighty_boyd       OK      a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run rnaseq-nf -with-docker
2019-05-06 12:31:15  1.2s      insane_celsius    ERR     b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf
2019-05-06 12:31:24  17s       stupefied_euclid  OK      b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf -resume -with-docker
```

You can use either the **session ID** or the **run name** to recover a specific execution:

```bash
nextflow run rnaseq-nf -resume mighty_boyd
```

### Execution provenance

The `log` command, when provided with a **run name** or **session ID**, can return many useful bits of information about a workflow execution that can be used to create a provenance report.

By default, it will list the work directories used to compute each task:

```console
nextflow log tiny_fermat
```

```console title="Output"
/data/.../work/7b/3753ff13b1fa5348d2d9b6f512153a
/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
/data/.../work/82/ba67e3175bd9e6479d4310e5a92f99
/data/.../work/e5/2816b9d4e7b402bfdd6597c2c2403d
/data/.../work/3b/3485d00b0115f89e4c202eacf82eba
```

The `-f` (fields) option can be used to specify which metadata should be printed by the `log` command:

```console
$ nextflow log tiny_fermat -f 'process,exit,hash,duration'
```

```console title="Output"
index    0   7b/3753ff  2.0s
fastqc   0   c1/56a36d  9.3s
fastqc   0   f7/659c65  9.1s
quant    0   82/ba67e3  2.7s
quant    0   e5/2816b9  3.2s
multiqc  0   3b/3485d0  6.3s
```

The complete list of available fields can be retrieved with the command:

```bash
nextflow log -l
```

The `-F` option allows the specification of filtering criteria to print only a subset of tasks:

```console
$ nextflow log tiny_fermat -F 'process =~ /fastqc/'
```

```console title="Output"
/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
```

This can be useful to locate specific task work directories.

Finally, the `-t` option enables the creation of a basic custom provenance report, showing a template file in any format of your choice:

```html
<div>
    <h2>${name}</h2>
    <div>
        Script:
        <pre>${script}</pre>
    </div>

    <ul>
        <li>Exit: ${exit}</li>
        <li>Status: ${status}</li>
        <li>Work dir: ${workdir}</li>
        <li>Container: ${container}</li>
    </ul>
</div>
```

!!! exercise

    Save the above snippet in a file named `template.html`. Then run this command (using the correct id for your run, e.g., `tiny_fermat`):

    ```bash
    nextflow log tiny_fermat -t template.html > prov.html
    ```

    Finally, open the `prov.html` file with a browser.

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How the workflow execution cache works
    2. How to use the `-resume` command line option
    3. How to organize _in-silico_ experiments in different folders
    4. How to create and customize a basic provenance report

## Troubleshooting resume

Being able to resume workflows is a key feature of Nextflow, but it doesn't always work as you expect. In this section you will learn common reasons why Nextflow may be ignoring your cached results.

!!! tip

    To learn more details about the resume mechanism and how to troubleshoot please refer to the following three blog posts:

     1. [Demystifying Nextflow resume](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html)
     2. [Troubleshooting Nextflow resume](https://www.nextflow.io/blog/2019/troubleshooting-nextflow-resume.html)
     3. [Analyzing caching behavior of pipelines](https://nextflow.io/blog/2022/caching-behavior-analysis.html)

### Input file changed

Make sure that there’s no change in your input file(s). Don’t forget the task unique hash is computed by taking into account the complete file path, the last modified timestamp and the file size. If any of this information has changed, the workflow will be re-executed even if the input content is the same.

### A process modifies an input

A process should never alter input files, otherwise the `resume` for future executions will be invalidated for the same reason explained in the previous point.

### Inconsistent file attributes

Some shared file systems, such as [NFS](https://en.wikipedia.org/wiki/Network_File_System), may report an inconsistent file timestamp (i.e. a different timestamp for the same file) even if it has not been modified. To prevent this problem use the [lenient cache strategy](https://www.nextflow.io/docs/latest/process.html#cache).

### Race condition in global variable

Nextflow is designed to simplify parallel programming without taking care about race conditions and the access to shared resources. One of the few cases in which a race condition can arise is when using a global variable with two (or more) operators

```groovy linenums="1" title="snippet.nf"
Channel
    .of(1, 2, 3)
    .map { it -> X = it; X += 2 }
    .view { "ch1 = $it" }

Channel
    .of(1, 2, 3)
    .map { it -> X = it; X *= 2 }
    .view { "ch2 = $it" }
```

The problem in this snippet is that the `X` variable in the closure definition is defined in the global scope. Therefore, since operators are executed in parallel, the `X` value can be overwritten by the other `map` invocation.

The correct implementation requires the use of the `def` keyword to declare the variable **local**.

```groovy linenums="1" title="snippet.nf"
Channel
    .of(1, 2, 3)
    .map { it -> def X = it; X += 2 }
    .println { "ch1 = $it" }

Channel
    .of(1, 2, 3)
    .map { it -> def X = it; X *= 2 }
    .println { "ch2 = $it" }
```

### Non-deterministic input channels

While dataflow channel ordering is guaranteed – data is read in the same order in which it’s written in the channel – be aware that there is no guarantee that the elements will maintain their order in the process _output_ channel. This is because a process may spawn multiple tasks, which can run in parallel. For example, the operation on the second element may end sooner than the operation on the first element, changing the output channel order.

In practical terms, consider the following snippet:

```groovy linenums="1" title="snippet.nf"
process FOO {
    input:
    val x

    output:
    tuple val(task.index), val(x)

    script:
    """
    sleep \$((RANDOM % 3))
    """
}

workflow {
    channel.of('A', 'B', 'C', 'D') | FOO | view
}
```

Just like you saw at the beginning of this tutorial with HELLO WORLD or WORLD HELLO, the output of the snippet above can be:

```console title="Output"
[0, A]
[3, C]
[4, D]
[2, B]
[1, A]
```

..and that order will likely be different every time the workflow is run.

Imagine that you now have two processes like this, whose output channels are acting as input channels to a third process. Both channels will be independently random, so the third process must not expect them to retain a paired sequence. If it does assume that the first element in the first process output channel is related to the first element in the second process output channel, there will be a mismatch.

A common solution for this is to use what is commonly referred to as a _meta map_. A groovy object with sample information is passed out together with the file results within an output channel as a tuple. This can then be used to pair samples from separate channels together for downstream use. For example, instead of putting just `/some/path/myoutput.bam` into a channel, you could use `['SRR123', '/some/path/myoutput.bam']` to make sure the processes are not incurring into a mismatch. Check the example below:

```groovy linenums="1" title="snippet.nf"
// For example purposes only.
// These would normally be outputs from upstream processes.
Channel
    .of(
        [[id: 'sample_1'], '/path/to/sample_1.bam'],
        [[id: 'sample_2'], '/path/to/sample_2.bam']
    )
    .set { bam }

// NB: sample_2 is now the first element, instead of sample_1
Channel
    .of(
        [[id: 'sample_2'], '/path/to/sample_2.bai'],
        [[id: 'sample_1'], '/path/to/sample_1.bai']
    )
    .set { bai }

// Instead of feeding the downstream process with these two channels separately, you can
// join them and provide a single channel where the sample meta map is implicitly matched:
bam
    .join(bai)
    | PROCESS_C
```

If meta maps are not possible, an alternative is to use the [`fair`](https://nextflow.io/docs/edge/process.html#fair) process directive. When this directive is specified, Nextflow will guarantee that the order of outputs will match the order of inputs (not the order in which the tasks run, only the order of the output channel).

!!! warning

     Depending on your situation, using the `fair` directive will lead to a decrease in performance.

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. Common reasons why Nextflow may be ignoring your cached results
