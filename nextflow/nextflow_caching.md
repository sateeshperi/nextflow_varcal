---
layout: main
title: Workflow caching and checkpointing
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_caching
---
{% include _nextflow_nextflow_caching_toc.html %}


<hr>
<center>This is part 6 of 8 of a <a href="/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

A key features of workflow management systems, like Nextflow, is re-entrancy which is the ability to restart a pipeline after an error from the last successful process. Re-entrancy can also skip time consuming sucessfully completed steps, such as index creation, when adding more data to a pipeline. Nextflow achieves re-entrancy by automatically keeping track of all the processes executed in your pipeline via caching and checkpointing.

## Resume

To restart from the last successfully executed process we add the command line option -resume to the Nextflow command.

For example, the command below would resume the wc.nf script from the last successful process.

```bash
cd ~/nextflow_tutorial
nextflow run word_count.nf --input 'data/untrimmed_fastq/*.fastq.gz' -resume
```
We can see in the output that the results from the process NUM_LINES has been retrieved from the cache.

>Output
>```bash
>N E X T F L O W  ~  version 21.04.3
>Launching `word_count.nf` [chaotic_knuth] - revision: 65cc76d410
>[dd/c19734] process > NUM_LINES (5) [100%] 6 of 6, cached: 6 ✔
>SRR2589044_2.fastq.gz4428360
>SRR2584863_1.fastq.gz6213036
>SRR2589044_1.fastq.gz4428360
>SRR2584866_1.fastq.gz11073592
>SRR2584866_2.fastq.gz11073592
>SRR2584863_2.fastq.gz6213036
>```

## How does resume work?

The mechanism works by assigning a unique ID to each task. This unique ID is used to create a separate execution directory, within the work directory, where the tasks are executed and the results stored. A task’s unique ID is generated as a 128-bit hash number obtained from a composition of the task’s:

*   Inputs values
*   Input files
*   Command line string
*   Container ID
*   Conda environment
*   Environment modules
*   Any executed scripts in the bin directory

When we resume a workflow Nextflow uses this unique ID to check if:

1.  The working directory exists
2.  It contains a valid command exit status
3.  It contains the expected output files.

If these conditions are satisfied, the task execution is skipped and the previously computed outputs are applied. When a task requires recomputation, ie. the conditions above are not fulfilled, the downstream tasks are automatically invalidated.

Therefore, if you modify some parts of your script, or alter the input data using `-resume`, will only execute the processes that are actually changed.

The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

## The Work directory

By default the pipeline results are cached in the directory work where the pipeline is launched.

We can use the Bash tree command to list the contents of the work directory. Note: By default tree does not print hidden files (those beginning with a dot .). Use the -a to view all files.

```bash
tree -a work
```

>Output
>```bash
>work
>├── 25
>│   └── a0d9448a5073401b9b1324d4f2bc5e
>│       ├── .command.begin
>│       ├── .command.err
>│       ├── .command.log
>│       ├── .command.out
>│       ├── .command.run
>│       ├── .command.sh
>│       ├── .exitcode
>│       └── SRR2584863_1.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz
>├── 3b
>│   └── 1d775f7d53553b8ec092c33373f481
>│       ├── .command.begin
>│       ├── .command.err
>│       ├── .command.log
>│       ├── .command.out
>│       ├── .command.run
>│       ├── .command.sh
>│       ├── .exitcode
>│       └── SRR2589044_1.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz
>├── 57
>│   └── 0f0c12f9fdb72115a7b2ed0bb078dd
>│       ├── .command.begin
>│       ├── .command.err
>│       ├── .command.log
>│       ├── .command.out
>│       ├── .command.run
>│       ├── .command.sh
>│       ├── .exitcode
>│       └── SRR2589044_2.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz
>├── b4
>│   └── 31ea937ef26c77f1ab688b40d1de15
>│       ├── .command.begin
>│       ├── .command.err
>│       ├── .command.log
>│       ├── .command.out
>│       ├── .command.run
>│       ├── .command.sh
>│       ├── .exitcode
>│       └── SRR2584863_2.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz
>├── ee
>│   └── 9f0a30bc02d4235033e25097aa2353
>│       ├── .command.begin
>│       ├── .command.err
>│       ├── .command.log
>│       ├── .command.out
>│       ├── .command.run
>│       ├── .command.sh
>│       ├── .exitcode
>│       └── SRR2584866_1.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz
>└── f6
>    └── bd1ee73d3fac620ade747ae749164a
>        ├── .command.begin
>        ├── .command.err
>        ├── .command.log
>        ├── .command.out
>        ├── .command.run
>        ├── .command.sh
>        ├── .exitcode
>        └── SRR2584866_2.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz
>
>12 directories, 48 file
>```

## Task execution directory

Within the work directory there are multiple task execution directories. There is one directory for each time a process is executed. These task directories are identified by the process execution hash. For example the task directory fa/cd3e49b63eadd6248aa357083763c1 would be location for the process identified by the hash fa/cd3e49 .

The task execution directory contains:

*   **`.command.sh`**: The command script.

*   **`.command.run`**: The command wrapped used to run the job.
   
*   **`.command.out`**: The complete job standard output.
   
*   **`.command.err`**: The complete job standard error.
   
*   **`.command.log`**: The wrapper execution output.
   
*   **`.command.begin`**: A file created as soon as the job is launched.
   
*   **`.exitcode`**: A file containing the task exit code.
   
*   Any task input files (symlinks)
   
*   Any task output files


## Specifying another work directory

Depending on your script, this work folder can take a lot of disk space. You can specify another work directory using the command line option -w. Note Using a different work directory will mean that any jobs will need to re-run from the beginning.

```bash
mkdir second_work_dir 
nextflow run word_count.nf --input 'data/untrimmed_fastq/*.fastq.gz' -w second_work_dir -resume
```

## Clean the work directory

If you are sure you won’t resume your pipeline execution, clean this folder periodically using the command nextflow clean.

```bash
nextflow clean [run_name|session_id] [options]
```

*   You need to specify the options -n to print names of file to be removed without deleting them, or -f to force the removal of the files. Note If you want to removes only temporary files but retains execution log entries and metadata you need to add the option -k.

*   If you want to clean the temporary files for multiple runs you can use the options, -before, -after or -but before the run name.

For example, the command below would remove all the temporary files and log entries for runs before the run gigantic_minsky.

```bash
nextflow clean -f -before gigantic_minsky
```


> Key Points
>*  Nextflow automatically keeps track of all the processes executed in your pipeline via checkpointing.
>*  Nextflow caches intermediate data in task directories within the work directory.
>*  Nextflow caching and checkpointing allows re-entrancy into a workflow after a pipeline error or using new data, skipping steps that have been successfully executed. - Re-entrancy is enabled using the -resume option.

---

<h5><a href="/nextflow/nextflow_modules" style="float: left"><b>Back to:</b>NextFlow Modules</a>

<a href="/nextflow/nextflow_install" style="float: right"><b>Next:</b>NextFlow Install</a></h5>