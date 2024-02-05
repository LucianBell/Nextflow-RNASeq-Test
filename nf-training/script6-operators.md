# Nextflow Operators

- In Nextflow, operators are used to **define the relationships between processes**
    - They basically specify how the output of one process should be connected to the input of another process.
- Every operator, with the exception of set and subscribe, **produces one or more new channels**

## Most common operators:
    Filtering: filter, randomSample, take, unique

    Reduction: collect, groupTuple, reduce

    Parsing text data: splitCsv, splitJson, splitText

    Combining channels: combine, concat, join, mix

    Forking channels: branch, multiMap

    Maths: count, max, min, sum

    Other: ifEmpty, map, set, view
---
