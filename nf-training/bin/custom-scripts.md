# Using custom scripts with Nextflow
- Real-world workflows use a lot of custom user scripts. 
- Nextflow allows you to consistently use and manage these scripts. 
    - Simply put them in a **directory named bin in the workflow project root.**
    - They will be automatically added to the workflow execution PATH.
- **REMEMBER:** You will need to give execute permission to the files in the bin folder
