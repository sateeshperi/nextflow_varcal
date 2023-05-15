---
layout: main
title: Workflow caching and checkpointing
categories: [nextflow]
tags: [cluster, nextflow, workflow, bioinformatics, tutorial]
permalink: /nextflow/nextflow_caching
---

{% include _nextflow_nextflow_caching_toc.html %}

<hr>
<center>This is part 6 of 8 of a <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

One of the essential features of workflow management systems, such as Nextflow, is their re-entrancy capability. This feature allows a pipeline to be restarted from the last successful process in case of an error, thereby saving a significant amount of time by skipping processes that have already been completed successfully. This is particularly beneficial in cases like adding more data to a pipeline, where index creation and other time-consuming steps can be skipped.

This re-entrancy is achieved in Nextflow through a combination of caching and checkpointing, allowing it to keep track of all executed processes automatically.

## Resume

To restart a pipeline from the last successful process, we use the `-resume` option in the Nextflow command.

For example, to resume the `word_count.nf` script from the last successful process, the command would be:

```bash
cd /workspace/nextflow_tutorial
nextflow run word_count.nf --input 'data/untrimmed_fastq/*.fastq.gz' -resume
```

In the output, we can see that the results from the `NUM_LINES` process have been retrieved from the cache.

> Output
>
> ```bash
> N E X T F L O W  ~  version 21.04.3
> Launching `word_count.nf` [chaotic_knuth] - revision: 65cc76d410
> [dd/c19734] process > NUM_LINES (5) [100%] 6 of 6, cached: 6 ✔
> SRR2589044_2.fastq.gz4428360
> SRR2584863_1.fastq.gz6213036
> SRR2589044_1.fastq.gz4428360
> SRR2584866_1.fastq.gz11073592
> SRR2584866_2.fastq.gz11073592
> SRR2584863_2.fastq.gz6213036
> ```

## How does resume work?

Nextflow assigns a unique ID to each task, which is used to create a separate execution directory within the work directory. This directory is where tasks are executed and results are stored. The unique ID is a 128-bit hash number derived from a combination of the task’s:

- Input values
- Input files
- Command line string
- Container ID
- Conda environment
- Environment modules
- Any executed scripts in the bin directory

Upon resuming a workflow, Nextflow uses this unique ID to verify:

1.  Whether the working directory exists
2.  Whether it contains a valid command exit status
3.  Whether it contains the expected output files

If these conditions are satisfied, the task execution is skipped and the previously computed outputs are used. If a task needs to be recomputed, i.e., if the conditions above are not met, the downstream tasks are automatically invalidated.

Therefore, if you modify some parts of your script or change the input data and use `-resume`, only the processes that have actually changed will be executed. The execution of the unchanged processes will be skipped, and the cached results will be used instead. This feature is incredibly helpful when testing or modifying a part of your pipeline without having to re-execute it from scratch.

## The Work directory

By default the pipeline results are cached in the directory work where the pipeline is launched.

We can use the Bash tree command to list the contents of the work directory. Note: By default tree does not print hidden files (those beginning with a dot .). Use the -a to view all files.

```bash
tree -a work
```

> Output
>
> ```bash
> work
> ├── 25
> │   └── a0d9448a5073401b9b1324d4f2bc5e
> │       ├── .command.begin
> │       ├── .command.err
> │       ├── .command.log
> │       ├── .command.out
> │       ├── .command.run
> │       ├── .command.sh
> │       ├── .exitcode
> │       └── SRR2584863_1.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz
> ├── 3b
> │   └── 1d775f7d53553b8ec092c33373f481
> │       ├── .command.begin
> │       ├── .command.err
> │       ├── .command.log
> │       ├── .command.out
> │       ├── .command.run
> │       ├── .command.sh
> │       ├── .exitcode
> │       └── SRR2589044_1.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz
> ├── 57
> │   └── 0f0c12f9fdb72115a7b2ed0bb078dd
> │       ├── .command.begin
> │       ├── .command.err
> │       ├── .command.log
> │       ├── .command.out
> │       ├── .command.run
> │       ├── .command.sh
> │       ├── .exitcode
> │       └── SRR2589044_2.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz
> ├── b4
> │   └── 31ea937ef26c77f1ab688b40d1de15
> │       ├── .command.begin
> │       ├── .command.err
> │       ├── .command.log
> │       ├── .command.out
> │       ├── .command.run
> │       ├── .command.sh
> │       ├── .exitcode
> │       └── SRR2584863_2.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz
> ├── ee
> │   └── 9f0a30bc02d4235033e25097aa2353
> │       ├── .command.begin
> │       ├── .command.err
> │       ├── .command.log
> │       ├── .command.out
> │       ├── .command.run
> │       ├── .command.sh
> │       ├── .exitcode
> │       └── SRR2584866_1.fastq.gz -> nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz
> └── f6
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
> 12 directories, 48 file
> ```

In Nextflow, every process execution is stored in a separate directory within the "work" directory. Each of these task directories is identified by a unique execution hash. For instance, the task directory "fa/cd3e49b63eadd6248aa357083763c1" would be associated with the process identified by the hash "fa/cd3e49".

Each task execution directory contains:

- **`.command.sh`**: The actual command script.
- **`.command.run`**: The job's command wrapper.
- **`.command.out`**: The complete standard output of the job.
- **`.command.err`**: The complete standard error of the job.
- **`.command.log`**: The output of the wrapper execution.
- **`.command.begin`**: A file created as soon as the job is launched.
- **`.exitcode`**: A file containing the task's exit code.
- Task input files (symlinks)
- Task output files

## Specifying Another Work Directory

Given that this "work" folder can consume significant disk space depending on your script, you can specify another directory for this purpose using the `-w` command-line option. Note, however, that using a different work directory will necessitate re-running any jobs from the start.

```bash
mkdir second_work_dir
nextflow run word_count.nf --input 'data/untrimmed_fastq/*.fastq.gz' -w second_work_dir -resume
```

## Cleaning the Work Directory

If you're certain that you won't resume your pipeline execution, you can periodically clean this folder using the `nextflow clean` command.

```bash
nextflow clean [run_name|session_id] [options]
```

Here, you can use the `-n` option to print names of files to be removed without actually deleting them, or `-f` to forcefully remove the files. Note that if you want to remove only temporary files while retaining execution log entries and metadata, you need to use the `-k` option.

To clean the temporary files for multiple runs, you can use the `-before`, `-after` or `-but` options before the run name.

For instance, the command below would remove all the temporary files and log entries for runs executed before the run "gigantic_minsky".

```bash
nextflow clean -f -before gigantic_minsky
```

> Key Points
>
> - Nextflow automatically tracks all processes executed in your pipeline via checkpointing.
> - Nextflow caches intermediate data in task directories within the work directory.
> - Nextflow's caching and checkpointing facilitate re-entrancy into a workflow after a pipeline error or when using new data, thus skipping steps that have been successfully executed. This re-entrancy is enabled using the `-resume` option.

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_reporting" style="float: left"><b>Back to:</b>NextFlow Reporting</a>

<a href="/nextflow_varcal/nextflow/nextflow_nfcore_variantcall" style="float: right"><b>Next:</b>nf-core/variantcall</a></h5>
