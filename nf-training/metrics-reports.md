# Metrics and Reports in Nextflow

- Nextflow can produce multiple reports and charts providing several runtime metrics and execution information. These can be enabled by using the following command line options:

        The -with-report option enables the creation of the workflow execution report.

        The -with-trace option enables the creation of a tab separated value (TSV) file containing runtime information for each executed task.

        The -with-timeline option enables the creation of the workflow timeline report showing how processes were executed over time. This may be useful to identify the most time consuming tasks and bottlenecks.

        Finally, the -with-dag option enables the rendering of the workflow execution direct acyclic graph representation. The dag needs to be given a name (-with-dag dag.png). Note: This feature requires the installation of Graphviz on your computer. See here for further details.