---
layout: main
title: NF-Core
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_nfcore
---
{% include _nextflow_nextflow_nfcore_toc.html %}


<hr>
<center>This is part 2 of 14 of a <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Learning Objectives
* Understand what nf-core is and how it relates to Nextflow.
* Use the nf-core helper tool to find nf-core pipelines.
* Run a nf-core pipeline (nf-core/bactmap) using a test dataset.
* Run a nf-core pipeline (nf-core/bactmap) using a test dataset on a High Performance Computing (HPC) environment.

<br>

## What is nf-core?

**[nf-core](https://nf-co.re/)** is a community effort to create a set of best-practice pipelines built using Nextflow. These pipelines are developed and maintained by experts in their respective fields, ensuring they are reliable, robust, and follow established guidelines. nf-core pipelines are designed to be easily portable and compatible with various computing environments, such as local machines, HPC clusters, and cloud platforms. Pipelines are governed by a **set of guidelines**, enforced by **community code reviews** and **automatic code testing**.

<img src="images/what_is_nf-core.png" alt="drawing" width="600"/>

<br>

## What is nf-core tools?

* nf-core provides a suite of helper tools aim to help users run and develop pipelines. These allow you to list all available pipelines and versions, with information about what versions you're running locally. There are also commands to help downloading pipelines for use offline.


### nf-core tools sub-commands

You can use the `--help` option to see the range of nf-core tools sub-commands.

```bash
nf-core --help
```

>```bash
                                          >,--./,-.
          >___     __   __   __   ___     /,-._.--~\
    >|\ | |__  __ /  ` /  \ |__) |__         }  {
    >| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          >`._,._,'
>
    >nf-core/tools version 2.8 - https://nf-co.re
>
>
>                                                                                                    
 >Usage: nf-core [OPTIONS] COMMAND [ARGS]...                                                         
>                                                                                                    
 >nf-core/tools provides a set of helper tools for use with nf-core Nextflow pipelines.              
 >It is designed for both end-users running pipelines and also developers creating new pipelines.    
>                                                                                                    
>╭─ Options ────────────────────────────────────────────────────────────────────────────────────────╮
>│ --version                        Show the version and exit.                                      │
>│ --verbose        -v              Print verbose output to the console.                            │
>│ --hide-progress                  Don't show progress bars.                                       │
>│ --log-file       -l  <filename>  Save a verbose log to a file.                                   │
>│ --help           -h              Show this message and exit.                                     │
>╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
>╭─ Commands for users ─────────────────────────────────────────────────────────────────────────────╮
>│ list        List available nf-core pipelines with local info.                                    │
>│ launch      Launch a pipeline using a web GUI or command line prompts.                           │
>│ download    Download a pipeline, nf-core/configs and pipeline singularity images.                │
>│ licences    List software licences for a given workflow (DSL1 only).                             │
>╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
>╭─ Commands for developers ────────────────────────────────────────────────────────────────────────╮
>│ create            Create a new pipeline using the nf-core template.                              │
>│ lint              Check pipeline code against nf-core guidelines.                                │
>│ modules           Commands to manage Nextflow DSL2 modules (tool wrappers).                      │
>│ subworkflows      Commands to manage Nextflow DSL2 subworkflows (tool wrappers).                 │
>│ schema            Suite of tools for developers to manage pipeline schema.                       │
>│ bump-version      Update nf-core pipeline version number.                                        │
>│ sync              Sync a pipeline TEMPLATE branch with the nf-core template.                     │
>╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
>```

## Listing Available nf-core Pipelines

The easiest way to see all available nf-core pipelines is by using the `nf-core list` command, which displays all pipelines from the nf-core GitHub repository. The output provides information on the latest version number, release date, and whether you have the latest version if the pipeline has been pulled locally using Nextflow.

To list all available nf-core pipelines, run:

```bash
nf-core list
```

>```bash
                                          >,--./,-.
          >___     __   __   __   ___     /,-._.--~\
    >|\ | |__  __ /  ` /  \ |__) |__         }  {
    >| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          >`._,._,'
>
    >nf-core/tools version 2.8 - https://nf-co.re
>
>
>┏━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
>┃ Pipeline Name          ┃ Stars ┃ Latest Release ┃      Released ┃ Last Pulled ┃ Have latest release? ┃
>┡━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
>│ smrnaseq               │    49 │          2.2.1 │    2 days ago │           - │ -                    │
>│ scrnaseq               │    84 │          2.3.0 │   1 weeks ago │           - │ -                    │
>│ funcscan               │    30 │          1.1.0 │   2 weeks ago │           - │ -                    │
>│ rnafusion              │    95 │          2.3.4 │   3 weeks ago │           - │ -                    │
>│ rnaseq                 │   610 │         3.11.2 │   3 weeks ago │ 3 weeks ago │ No (v3.11.1)         │
>│ demultiplex            │    26 │          1.2.0 │   3 weeks ago │           - │ -                    │
>│ differentialabundance  │    19 │          1.2.0 │   3 weeks ago │           - │ -                    │
>│ mhcquant               │    21 │          2.4.1 │  1 months ago │           - │ -                    │
>│ viralintegration       │     8 │          0.1.0 │  2 months ago │           - │ -                    │
>│ quantms                │     9 │          1.1.1 │  2 months ago │           - │ -                    │
>│ viralrecon             │    93 │          2.6.0 │  2 months ago │           - │ -                    │
>│ airrflow               │    24 │            3.0 │  2 months ago │           - │ -                    │
>│ epitopeprediction      │    25 │          2.2.1 │  2 months ago │           - │ -                    │
>│ isoseq                 │    12 │          1.1.4 │  2 months ago │           - │ -                    │
>│ taxprofiler            │    50 │          1.0.0 │  2 months ago │           - │ -                    │
>│ nanoseq                │   110 │          3.1.0 │  2 months ago │           - │ -                    │
>│ cutandrun              │    41 │            3.1 │  2 months ago │           - │ -                    │
>│ circdna                │    13 │          1.0.2 │  2 months ago │           - │ -                    │
>│ ampliseq               │   114 │          2.5.0 │  2 months ago │           - │ -                    │
>│ mag                    │   126 │          2.3.0 │  2 months ago │           - │ -                    │
>│ nascent                │     8 │          2.1.1 │  3 months ago │           - │ -                    │
>│ phyloplace             │     3 │          1.0.0 │  3 months ago │           - │ -                    │
>│ proteinfold            │    21 │          1.0.0 │  3 months ago │           - │ -                    │
>│ crisprseq              │     8 │            1.0 │  3 months ago │           - │ -                    │
>│ hic                    │    51 │          2.0.0 │  4 months ago │           - │ -                    │
>│ sarek                  │   239 │          3.1.2 │  4 months ago │           - │ -                    │
>│ fetchngs               │    78 │            1.9 │  5 months ago │           - │ -                    │
>│ methylseq              │   104 │          2.3.0 │  5 months ago │ 1 weeks ago │ Yes (v2.3.0)         │
>│ atacseq                │   134 │            2.0 │  5 months ago │           - │ -                    │
>│ eager                  │    91 │          2.4.6 │  6 months ago │           - │ -                    │
>│ coproid                │     7 │          1.1.1 │  6 months ago │           - │ -                    │
>│ hgtseq                 │    17 │          1.0.0 │  7 months ago │           - │ -                    │
>│ hlatyping              │    42 │          2.0.0 │  7 months ago │           - │ -                    │
>│ chipseq                │   145 │          2.0.0 │  7 months ago │           - │ -                    │
>│ rnavar                 │    16 │          1.0.0 │ 11 months ago │           - │ -                    │
>│ mnaseseq               │     9 │          1.0.0 │ 12 months ago │           - │ -                    │
>│ hicar                  │     3 │          1.0.0 │   1 years ago │           - │ -                    │
>│ bamtofastq             │     9 │          1.2.0 │   1 years ago │           - │ -                    │
>│ bacass                 │    42 │          2.0.0 │   2 years ago │           - │ -                    │
>│ bactmap                │    41 │          1.0.0 │   2 years ago │           - │ -                    │
>│ metaboigniter          │    10 │          1.0.1 │   2 years ago │           - │ -                    │
>│ diaproteomics          │    10 │          1.2.4 │   2 years ago │           - │ -                    │
>│ clipseq                │    13 │          1.0.0 │   2 years ago │           - │ -                    │
>│ pgdb                   │     3 │          1.0.0 │   2 years ago │           - │ -                    │
>│ dualrnaseq             │    12 │          1.0.0 │   2 years ago │           - │ -                    │
>│ cageseq                │     9 │          1.0.2 │   2 years ago │           - │ -                    │
>│ proteomicslfq          │    29 │          1.0.0 │   3 years ago │           - │ -                    │
>│ imcyto                 │    20 │          1.0.0 │   3 years ago │           - │ -                    │
>│ slamseq                │     4 │          1.0.0 │   3 years ago │           - │ -                    │
>│ callingcards           │     1 │            dev │             - │           - │ -                    │
>│ circrna                │    27 │            dev │             - │           - │ -                    │
>│ fastquorum             │     8 │            dev │             - │           - │ -                    │
>│ genomeannotator        │     9 │            dev │             - │           - │ -                    │
>│ genomeassembler        │    12 │            dev │             - │           - │ -                    │
>│ gwas                   │    12 │            dev │             - │           - │ -                    │
>│ lncpipe                │    25 │            dev │             - │           - │ -                    │
>│ metapep                │     3 │            dev │             - │           - │ -                    │
>│ metatdenovo            │     3 │            dev │             - │           - │ -                    │
>│ nanostring             │     2 │            dev │             - │           - │ -                    │
>│ pangenome              │    23 │            dev │             - │           - │ -                    │
>│ radseq                 │     0 │            dev │             - │           - │ -                    │
>│ raredisease            │    37 │            dev │             - │           - │ -                    │
>│ rnadnavar              │     0 │            dev │             - │           - │ -                    │
>│ rnasplice              │     3 │            dev │             - │           - │ -                    │
>│ scflow                 │    19 │            dev │             - │           - │ -                    │
>│ spatialtranscriptomics │    19 │            dev │             - │           - │ -                    │
>│ spinningjenny          │     0 │            dev │             - │           - │ -                    │
>│ variantcatalogue       │     4 │            dev │             - │           - │ -                    │
>└────────────────────────┴───────┴────────────────┴───────────────┴─────────────┴──────────────────────┘
>```

---

## Running nf-core Pipelines

## Usage Instructions and Documentation

*   General documentation and instructions for Nextflow and nf-core can be found on the [nf-core website](https://nf-co.re/). Each pipeline has pipeline-specific documentation bundled within the `/docs` folder. You can read this documentation locally, on GitHub, or on the nf-core website.

*   Every pipeline has its own webpage on the nf-core website, such as [nf-co.re/bactmap](https://nf-co.re/bactmap).

*   Alongside this documentation, each pipeline provides a basic command line reference. You can access this by running the pipeline with the `--help` flag, like this:

```bash
nextflow run nf-core/viralintegration -r 0.1.0 --help
```

>```bash
>N E X T F L O W  ~  version 23.04.1
>Pulling nf-core/viralintegration ...
> downloaded from https://github.com/nf-core/viralintegration.git
>Launching `https://github.com/nf-core/viralintegration` [drunk_cantor] DSL2 - revision: 88a9d1708a [0.1.0]
>
>
>------------------------------------------------------
>                                        ,--./,-.
>        ___     __   __   __   ___     /,-._.--~'
>  |\ | |__  __ /  ` /  \ |__) |__         }  {
>  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
>                                        `._,._,'
>  nf-core/viralintegration v0.1.0-g88a9d17
>------------------------------------------------------
>Typical pipeline command:
>
>  nextflow run nf-core/viralintegration --input samplesheet.csv --genome GRCh37 -profile docker
>
>Input/output options
>  --input                       [string]  Path to comma-separated file containing information about the samples in the experiment.
>  --viral_fasta                 [string]  Path to fasta file with viral genomes [default: 
>                                          https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/VIRUS_INSERTION_FINDING_LIB_SUPPLEMENT/virus_db.fasta] 
>  --outdir                      [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud 
>                                          infrastructure. 
>  --email                       [string]  Email address for completion summary.
>  --multiqc_title               [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.
>  --min_reads                   [integer] Minimum reads [default: 5]
>  --max_hits                    [integer] Maximum hits [default: 50]
>  --remove_duplicates           [boolean] To remove duplicate viral hits or not [default: true]
>
>Reference genome options
>  --genome                      [string]  Name of iGenomes reference.
>  --fasta                       [string]  Path to FASTA genome file.
>  --gtf                         [string]  Path to GTF annotation file.
>
>Generic options
>  --multiqc_methods_description [string]  Custom MultiQC yaml file containing HTML including a methods description.
>
>!! Hiding 24 params, use --show_hidden_params to show them !!
>------------------------------------------------------
>If you use nf-core/viralintegration for your analysis please cite:
>
>* The nf-core framework
>  https://doi.org/10.1038/s41587-020-0439-x
>
>* Software dependencies
>  https://github.com/nf-core/viralintegration/blob/master/CITATIONS.md
>------------------------------------------------------
>```

<br>

## Running Pipelines with Test Profile

*   The `test` config profile specifies URLs for test data and all required parameters. This allows you to test any nf-core pipeline using the following command:

```bash
nextflow run nf-core/viralintegration -r 0.1.0 -profile test,docker --outdir test_results
```

---

>**If you encounter an error saying `Failed to pull singularity image` with `status:255`, follow the instructions below to create a `/tmp` folder in your `$HOME` directory and re-run the `nextflow run` command:**
>
>```bash
>mkdir -p $HOME/tmp && export TMPDIR=$HOME/tmp
>```

---

sample output should looks like this:

>```bash
>N E X T F L O W  ~  version 23.04.1
>Launching `https://github.com/nf-core/viralintegration` [clever_lagrange] DSL2 - revision: 88a9d1708a [0.1.0]
>
>
>------------------------------------------------------
>                                        ,--./,-.
>        ___     __   __   __   ___     /,-._.--~'
>  |\ | |__  __ /  ` /  \ |__) |__         }  {
>  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
>                                        `._,._,'
>  nf-core/viralintegration v0.1.0-g88a9d17
>------------------------------------------------------
>Core Nextflow options
>  revision                  : 0.1.0
>  runName                   : clever_lagrange
>  containerEngine           : docker
>  launchDir                 : /home/speri/github/nextflow_varcal
>  workDir                   : /home/speri/github/nextflow_varcal/work
>  projectDir                : /home/speri/.nextflow/assets/nf-core/viralintegration
>  userName                  : speri
>  profile                   : test,docker
>  configFiles               : /home/speri/.nextflow/assets/nf-core/viralintegration/nextflow.config
>
>Input/output options
>  input                     : /home/speri/.nextflow/assets/nf-core/viralintegration/assets/samplesheet.csv
>  viral_fasta               : https://raw.githubusercontent.com/broadinstitute/CTAT-VirusIntegrationFinder/master/testing/HPV16.fa
>  outdir                    : test_results
>  remove_duplicates         : true
>
>Reference genome options
>  genome                    : GRCh38
>  fasta                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralintegration/testdata/GRCh38_chr18.fa
>  gtf                       : https://raw.githubusercontent.com/nf-core/test-datasets/viralintegration/testdata/GRCh38_chr18.gtf
>
>Institutional config options
>  config_profile_name       : Test profile
>  config_profile_description: Minimal test dataset to check pipeline function
>
>Max job request options
>  max_cpus                  : 2
>  max_memory                : 6.GB
>  max_time                  : 6.h
>
>!! Only displaying parameters that differ from the pipeline defaults !!
>------------------------------------------------------
>If you use nf-core/viralintegration for your analysis please cite:
>
>* The nf-core framework
>  https://doi.org/10.1038/s41587-020-0439-x
>
>* Software dependencies
>  https://github.com/nf-core/viralintegration/blob/master/CITATIONS.md
>------------------------------------------------------
>executor >  local (45)
>[ed/a2ba7c] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:INPUT_CHECK:SAMPLESHEET_CHECK (samplesheet.csv)            [100%] 1 of 1 ✔
>[07/65a563] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:FASTQC (redherring_T1)                                     [100%] 3 of 3 ✔
>[e9/933cd3] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:STAR_GENOMEGENERATE_HOST (GRCh38_chr18.fa)                 [100%] 1 of 1 ✔
>[15/2b5059] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:STAR_ALIGN_HOST (test_T1)                                  [100%] 3 of 3 ✔
>[25/bf8d6c] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:TRIMMOMATIC (test_T1)                                      [100%] 3 of 3 ✔
>[4d/14ba97] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:POLYA_STRIPPER (test_T1)                                   [100%] 3 of 3 ✔
>[ce/a01528] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:CAT_FASTA (GRCh38_chr18.fa)                                [100%] 1 of 1 ✔
>[49/ae91ad] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:STAR_GENOMEGENERATE_PLUS (GRCh38_chr18_plus_viraldb.fasta) [100%] 1 of 1 ✔
>[82/b7b59c] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:STAR_ALIGN_PLUS (test_T1)                                  [100%] 3 of 3 ✔
>[b3/049d18] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:SAMTOOLS_SORT_PLUS (test_T1)                               [100%] 3 of 3 ✔
>[34/1d56e9] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:SAMTOOLS_INDEX_PLUS (test_T1)                              [100%] 3 of 3 ✔
>[53/26ccb7] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:INSERTION_SITE_CANDIDATES (test_T1)                        [100%] 3 of 3 ✔
>[b2/0a7f19] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:ABRIDGED_TSV (test_T1)                                     [100%] 3 of 3 ✔
>[fb/10e9f7] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:VIRUS_REPORT (test_T1)                                     [100%] 3 of 3 ✔
>[be/3e71ed] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:EXTRACT_CHIMERIC_GENOMIC_TARGETS (test_T1)                 [100%] 3 of 3 ✔
>[0f/d0b6ac] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:STAR_ALIGN_VALIDATE (test_T1)                              [100%] 1 of 1 ✔
>[5e/2e4a1a] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:SAMTOOLS_SORT_VALIDATE (test_T1)                           [100%] 1 of 1 ✔
>[e5/9585bf] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:SAMTOOLS_INDEX_VALIDATE (test_T1)                          [100%] 1 of 1 ✔
>[95/8c1f8c] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:REMOVE_DUPLICATES (test_T1)                                [100%] 1 of 1 ✔
>[2f/59a8b2] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:CHIMERIC_CONTIG_EVIDENCE_ANALYZER (test_T1)                [100%] 1 of 1 ✔
>[76/21d361] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:SUMMARY_REPORT (test_T1)                                   [100%] 1 of 1 ✔
>[45/4dc88a] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:CUSTOM_DUMPSOFTWAREVERSIONS (1)                            [100%] 1 of 1 ✔
>[b9/8e6aff] process > NFCORE_VIRALINTEGRATION:VIRALINTEGRATION:MULTIQC                                                    [100%] 1 of 1 ✔
>-[nf-core/viralintegration] Pipeline completed successfully-
>Completed at: 14-May-2023 15:30:41
>Duration    : 5m 26s
>CPU hours   : 0.2
>Succeeded   : 45
>```

Running with the test profile is a great way to confirm that you have Nextflow configured properly for your system before attempting to run with real data. This helps ensure that the pipeline, dependencies, and execution environment are set up correctly, making it easier to troubleshoot any issues that may arise when running with your actual data.

The `results` folder will contain:

```bash
test_results/
├── abridged
├── cat
├── chimeric
├── extract
├── fastqc
├── insertion
├── multiqc
├── pipeline_info
├── polya
├── remove
├── samtools
├── star
├── summary
├── trimmomatic
└── virus
```

---

* Many of the techniques and resources described above require an active internet connection at runtime - pipeline files, configuration profiles, and software containers are all dynamically fetched when the pipeline is launched. This can be a problem for people using secure computing resources that do not have connections to the internet.

* To help with this, the `nf-core download` command automates the fetching of required files for running nf-core pipelines offline. The command can download a specific release of a pipeline with `-r/--release`. By default, the pipeline will download the pipeline code and the institutional nf-core/configs files.

```bash
nf-core download nf-core/viralintegration -r 0.1.0
```

The command prompt will ask whether to download the Singularity containers or none. Choose **`singularity`**

```bash
? Download software container images: (Use arrow keys)
   none
 » singularity
```

Next step, choose **`n`** to define a shared Singularity cache directory.

```bash
? Define $NXF_SINGULARITY_CACHEDIR for a shared Singularity image download folder? [y/n]: n
```

Next, for the compression step, choose **`none`**

```bash
? Choose compression type: (Use arrow keys)
 » none
   tar.gz
   tar.bz2
   zip
```

This will download the necessary files, containers, and configurations to run the specified nf-core pipeline offline. Be sure to transfer the downloaded files to the offline environment and configure the appropriate paths before executing the pipeline.

```bash

```

Now you can run the nf-core/viralintegration pipeline with the provided test data using the test profile

```bash
cd nf-core-viralintegration-1.0.0/workflow/
nextflow run main.nf -profile test,singularity
```

<br>

## Troubleshooting

If you encounter issues running your pipeline, you can use the nf-core website to [troubleshoot common mistakes and issues](https://nf-co.re/usage/troubleshooting).

## Extra resources and getting help

* If you still have an issue with running the pipeline, feel free to contact the nf-core community via the Slack channel. The nf-core Slack organization has channels dedicated to each pipeline, as well as specific topics (e.g. #help, #pipelines, #tools, #configs, and much more).

* The nf-core Slack can be found at `https://nfcore.slack.com` (NB: no hyphen in nfcore!). To join, you will need an invite, which you can get at `https://nf-co.re/join/slack`.

* You can also get help by opening an issue in the respective pipeline repository on GitHub, asking for assistance.

* If you have problems that are directly related to Nextflow and not our pipelines or the nf-core framework tools, check out the [Nextflow Gitter channel](https://gitter.im/nextflow-io/nextflow) or the [Google group](https://groups.google.com/forum/#!forum/nextflow).

---

> Citation: If you use nf-core tools in your work, please cite the nf-core publication as follows:
>> The nf-core framework for community-curated bioinformatics pipelines. Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen. Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x. ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

## Quick Recap
>* nf-core is a community-led project to develop a set of best-practice pipelines built using the Nextflow workflow management system.
>* nf-core tools is a suite of helper tools that aims to help people run and develop pipelines.
>* nf-core pipelines can be found using the nf-core helper tool `--list` option or from the nf-core website.
>* nf-core pipelines can be run using the `nextflow run nf-core/<pipeline>` syntax, or launched and parameters configured using the nf-core helper tool launch option.
---

<h5><a href="/nextflow_varcal/nextflow/nextflow_intro" style="float: left"><b>Back to:</b>Introduction</a>

<a href="/nextflow_varcal/nextflow/nextflow_slurm" style="float: right"><b>Next:</b>NF-Core @ HPC</a></h5>
