---
layout: main
title: nf-core module creation
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_nfcore_module_create
---
{% include _nextflow_nextflow_nfcore_module_create_toc.html %}


<hr>
<center>This is part of <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>


## Overview

* Initial Environment Setup (`conda`)
* Initialize a new `nf-core` pipeline
* Install or create `nf-core` modules needed for the pipeline
* Create subworkflows for different pipeline steps/phases
* Create a primary workflow
* Modify `config` files for the pipeline components
* Test the pipeline
* Deployment, distribution, and support

<br>

# 1.   Initial environment setup (optional; recommended for local setup)

## *(requires conda/mamba to be installed prior to beginning)*

<br>

## Create a Conda environment for Nextflow
```Bash
# Create a new empty environment
mamba create -n nextflow
mamba install -c conda-forge -c bioconda nextflow=21.10.6 nf-core=2.2 graphviz openjdk=8.0.312 git=2.35.0

# Activate the environment
conda activate nextflow

# Create a clean, sharable copy of this conda environment
conda env export | grep -v "prefix" > env.nextflow.yml
```


# Creating a New `nf-core` Pipeline from Scratch

This uses templates standardized by `nf-core tools`, which are available from the `nf-core` package in our conda environment

```Bash
# This will activate an interactive prompt to initialize the pipeline name, author, and git repo
nf-core create

# $ Workflow Name: quaisar
# $ Description: nf-core demo workflow
# $ Author: hseabolt

cd nf-core-quaisar

# This is the new repo for this pipeline -- many files and template files
tree .

.
├── assets
│   ├── email_template.html
│   ├── email_template.txt
│   ├── multiqc_config.yaml
│   ├── nf-core-quaisar_logo_light.png
│   ├── samplesheet.csv
│   ├── schema_input.json
│   └── sendmail_template.txt
├── bin
│   └── check_samplesheet.py
├── CHANGELOG.md
├── CITATIONS.md
├── CODE_OF_CONDUCT.md
├── conf
│   ├── base.config
│   ├── igenomes.config
│   ├── modules.config
│   ├── test.config
│   └── test_full.config
├── docs
│   ├── images
│   │   ├── mqc_fastqc_adapter.png
│   │   ├── mqc_fastqc_counts.png
│   │   ├── mqc_fastqc_quality.png
│   │   ├── nf-core-quaisar_logo_dark.png
│   │   └── nf-core-quaisar_logo_light.png
│   ├── output.md
│   ├── README.md
│   └── usage.md
├── lib
│   ├── nfcore_external_java_deps.jar
│   ├── NfcoreSchema.groovy
│   ├── NfcoreTemplate.groovy
│   ├── Utils.groovy
│   ├── WorkflowMain.groovy
│   └── WorkflowQuaisar.groovy
├── LICENSE
├── main.nf
├── modules
│   ├── local
│   │   └── samplesheet_check.nf
│   └── nf-core
│       └── modules
│           ├── custom
│           │   └── dumpsoftwareversions
│           │       ├── main.nf
│           │       ├── meta.yml
│           │       └── templates
│           │           └── dumpsoftwareversions.py
│           ├── fastqc
│           │   ├── main.nf
│           │   └── meta.yml
│           └── multiqc
│               ├── main.nf
│               └── meta.yml
├── modules.json
├── nextflow.config
├── nextflow_schema.json
├── README.md
├── subworkflows
│   └── local
│       └── input_check.nf
└── workflows
    └── quaisar.nf
```


Edit your modules, sub-workflows, workflows, etc. and then don't forget to add

```Bash
git add --all
git commit -m "Initial commit of quaisar nf-core pipeline"

# Prior to pushing to an online repo, make sure to create the repo online and don't add any files to it
git remote add origin <github repo URL>
git push origin master
```


# nf-core Modules 101

* nf-core has many modules already available for a variety of bioinformatics software, which can be downloaded and installed directly into the new pipeline with `nf-core` tools.  To check which modules are already created (and peer-reviewed!), use the following command(s):

<br>

## Checking which modules are available in all of nf-core

```Bash
# Can scroll thru this on the terminal or pipe to grep to check a specific module
nf-core modules list remote

# Grep example
nf-core modules list remote | grep "bwa/align"
```
<br>

## Install an available module direct from nf-core

* nf-core provides an easy interface for module installs. Most softwares consist of a single command plus command-line parameter arguments -- for these, the modules can be installed directly.  However, some softwares have multiple sub-commands, eg. `bwa index`, `bwa align`, `bwa mem`, etc.  In these cases, `nf-core` requires that each sub-command be created as a standalone module within the parent module, which can also be directly installed.  Beware that just installing the parent module (ie. `bwa` in this case) through `nf-core` does not give you direct access to sub-modules.

```bash
# Install nf-core modules from nf-core remote (master) repo
nf-core modules install samtools_stats

# Install modules that have specific sub-modules 
nf-core modules install bwa/index
```

# Creating new modules for nf-core
## TL;DR i just want to write code

* Join nf-core Slack (https://nf-co.re/join).
* `Fork` the `nf-core/modules` repo (https://github.com/nf-core/modules) to your own Github repos, then `git clone` it down locally from your `fork`.
* Create a new `branch` for the new module within your local clone of `nf-core/modules`.
* Raise an issue on the main nf-core Github (https://github.com/nf-core) for the new module to alert others that you are working on this module (`Issues` --> `New Module`.  Add yourself to the `Assignees`).
* Check that a `conda` recipe, or `Docker`/`Singularity` image already available in `Bioconda`, `Biocontainers`, `quay.io`, or similar for the software you want to create a module for.
* Use `nf-core` tools to create, edit, and test the module.
* `git push` your module's `branch` to your `nf-core/modules` fork.
* Open a `pull request` on the main `nf-core/modules` repo.
* Request in the `nf-core` Slack (channel `request-review` ) that someone review the `pull request` and approve the `merge`.

<br>

## Fork the nf-core modules Github repo

Navigate to https://github.com/nf-core/modules and click `Fork` in the top right corner of the web browser. This will create a forked copy of the entire `nf-core/modules` repo to your own Github repo.  

Next, we need to let other `nf-core` developers know that we are working on creating new modules so that we do not duplicate efforts or step over the work of others.  To do this, in the main `nf-core/modules` Github repo, open an `Issue` (https://github.com/nf-core/modules/issues) with `New Issue` --> `New Module`.  


*(Don't forget to search the open modules to see if someone is already working on the module you want to add!)*

* Read over the `markdown` text in the dialogue box for new modules and fill out/edit appropriately.  
* Assign yourself in the `Assignees` in the top right corner.
* When you are finished, then click `Submit New Issue`.

<br>

Now, we need to clone the repo down to local storage so that we can branch/modify it.

```bash
# From somewhere that you want to clone the fork of the modules repo from your Git.
# DO NOT do this inside the workflow that you are creating
git clone https://github.com/hseabolt/modules.git
cd modules

# The modules directory should look something like this
tree .
.
├── docs
│   └── images
│       └── nfcore-modules_logo.png
├── LICENSE
├── main.nf
├── modules
│   ├── abacas
│   │   ├── main.nf
│   │   └── meta.yml
│   ├── abricate
│   │   ├── run
│   │   │   ├── main.nf
│   │   │   └── meta.yml
│   │   └── summary
│   │       ├── main.nf
│   │       └── meta.yml
│   ├── adapterremoval
│   │   ├── main.nf
│   │   └── meta.yml
...

# Create a new branch for the module you want to create
# Here I am doing Nonpareil, which is a software used for metagenomics analyses.
git checkout -b nonpareil

# IMPORTANT: before modifying or creating ANY new files, always make check to make sure that you are not working in the master branch!!!
nf-core modules create nonpareil --author @hseabolt --label process_low --meta

# The new directories for the new module are created in ./modules and ./tests
#    ./modules/nonpareil/main.nf
#    ./modules/nonpareil/meta.yml
#    ./tests/modules/nonpareil/main.nf
#    ./tests/modules/nonpareil/test.yml
#    ./tests/modules/nonpareil/nextflow.config
#    ./tests/config/pytest_modules.yml
```

<br>

### Where do I write my actual code?

<br>

At this point, you are ready to write the module code itself in `./modules/main.nf`, which includes a number of `TODO` statements that must all be addressed for the module creation to be completed successfully and merged into `nf-core`.
* The main code for running the software going in the `./modules/nonpareil/main.nf` file, including the links to the `Singularity` containers, inputs and outputs, and any other details.
* It is very important to include the code block at the bottom to capture the `Version` number of the module software (only the version number, not a whole line of info).
* If there are any special parameter args (e.g. `-t 1 -p 2`, etc.), these should be included in the `./tests/modules/nonpareil/nextflow.config` and **NOT** in the module code in `main.nf`.

<br>

An example `main.nf` to run the module:
```Groovy
process NONPAREIL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::nonpareil=3.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil%3A3.4.1--r41h9f5acd7_1' :
        'quay.io/biocontainers/nonpareil' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqtk \\
        -s $reads \\
        $args \\
		| gzip --no-name > ${prefix}.fastq.gz
			
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(echo \$(nonpareil -V) | cut -f2 | sed 's/v//')
    END_VERSIONS
    """
}
```

<br>

### The accompanying `meta.yml` file describes the new module in detail

<br>

The `./modules/nonpareil/meta.yml` file serves as a documentation companion to `./modules/nonpareil/main.nf` and describes the module's required inputs and outputs (which **must** match those inside the `main.nf`!) and also includes details on the module's authors, sources, and descriptions.  The descriptors can usually be found on the original Github or similar where the source code for the software is hosted.

An example `meta.yml` file for this module (`./modules/nonpareil/meta.yml`)

```yml
name: nonpareil
description: Estimate metagenomic coverage and sequence diversity
keywords:
  - diversity
  - metagenomics
  - nonpareil
  - kmer
tools:
  - nonpareil:
      description: Nonpareil uses the redundancy of the reads in metagenomic datasets to estimate the average coverage and predict the amount of sequences that will be required to achieve 'nearly complete coverage'.
      homepage: https://github.com/lmrodriguezr/nonpareil
      documentation: https://nonpareil.readthedocs.io/en/latest/index.html
      tool_dev_url: https://github.com/lmrodriguezr/nonpareil
      licence: ['Artistic Licence 2.0']

input:
  - meta:
      type: map
      description: |
        Groovy Map containing sample information
        e.g. [ id:'test', single_end:false ]
  - reads:
      type: file
      description: List of input FastQ files of size 1 and 2 for single-end and paired-end data,respectively.
      pattern: "*.{fastq.gz}"

output:
  - meta:
      type: map
      description: |
        Groovy Map containing sample information
        e.g. [ id:'test', single_end:false ]
  - versions:
      type: file
      description: File containing software versions
      pattern: "versions.yml"
  - reads:
      type: file
      description: Subsampled FastQ files, 1 for single-end data or 2 for paired-end data.
      pattern: "*.{fastq.gz}"

authors:
  - "@hseabolt"
```

### How do I locally test my module?

<br>

* Once the module code `main.nf` and `meta.yml` have been written and you are happy with it, it is time to test our module code locally.

* Next, identify some suitable testing data for your module.  `nf-core` has a sizable set of available data in multiple common formats used in bioinformatics available in the `nf-core/test-datasets` repo on their main Github. These can be easily inserted into `nf-core` module testing using a meta map without having to download any data directly.
* Check the list of available data in the file `./modules/tests/config/test_data.config`.

<br>
 
Edit the `./tests/modules/nonpareil/main.nf` to include the test data you identified for your module and ensure that the data and args match what is expected in the main code `./modules/nonpareil/main.nf`.  For example, if you wanted to include a gzipped FASTQ file from the available `nf-core` data sets (e.g. `test_1_fastq_gz` on line 53), your `./modules/tests/modules/nonpareil/main.nf` file might look like this:

```Groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NONPAREIL } from '../../../../modules/nonpareil/main.nf'

// Test with single-end data
workflow test_seqtk_sample_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]

    NONPAREIL ( input, 50 )
}

// Test with paired-end data
workflow test_seqtk_sample_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    NONPAREIL ( input, 50 )
}
```


If you want to include/test any optional parameter arguments for the module, these can be specified in the `./modules/tests/modules/nonpareil/nextflow.config` file:  An example of this using an arg `-T kmer` and setting the `prefix` parameter that is fed into the main code might look something like this:

<br>

```Groovy
process {

    publishDir = { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" }

    withName: NONPAREIL {
        ext.args = '-T kmer' // extend this string to include all args, do not use separate strings
        ext.prefix = { "${meta.id}.np" }
    }

}
```


Now you are ready to begin testing your module for `nf-core`. 


```bash
# Start to test the module (in the root folder of nf-core/modules)
nf-core modules create-test-yml nonpareil

# Will create an interactive prompt on the terminal.
# Just press ENTER / leave blank most things.
# For testing profile: use Singularity profile

# If the testing is successful, a new file will be created: tests/modules/nonpareil/test.yml

# Now run local linting which will give you more specific details on the tests run.
nf-core lint nonpareil --dir .

# If all tests and linting are successful (or in some cases, safely ignorable), then we are ready to commit the code to git and initiate merge requests
git status
git add --all  # or instead of --all, add the files listed by git status
git commit -m "Created new nf-core module Nonpareil"
git push -u origin nonpareil
```

If you want to manually check the testing output files for correctness (eg. FASTA files created by the module), `cd` into the `work` directory specified in the terminal output from `create-test-yml` and verify the output files are what you expect.

<br>

### How do I merge my local module code with the main `nf-core` repo for others to use?

<br>

* From here, you will need to create a `Pull Request` (`PR`) across forks with your local copy of `nf-core/modules`.

* On the `nf-core/modules` Github page, navigate to the `Pull requests` tab, then choose `New pull request`.  Fill out the dialogue box for the new pull request, get the `Issue Number` from the original `New Module Issue` that you opened previously, and include this to let others know which `Issue` you are closing out.  Double check that you are merging `hseabolt/modules`->`nonpareil` into `nf-core/modules`->`master`.  

*  Once a new `PR` is created, `nf-core` Github kicks off some automated integration and testing functionality. You must wait for all of these checks to complete and pass prior to proceeding.     
    * If any tests do not pass, you must address these before your new code will be merged.  Most often this is small things like whitespace or minor formatting issues.  ( As you address these, you must re-run `git add`, `git commit`, `git push` to your local `fork`, which will automatically feed to the open `PR` with `nf-core` and automatically re-trigger the integration tests.  
    * You do not need to create a new `PR`.  
    * Once module tests pass, add `Labels` (`New Module`, `Ready_For_Review`) on the right side of the `PR` dashboard in Git).

<br>

### Request peer review from other `nf-core` developers:

<br>

* Copy the `PR`'s URL link and paste in the `nf-core` Slack channel `request-review` and politely request someone to peer-review your new module.  Once a reviewer agrees to review, they will be able to make comments or ask questions through Github that you must address satisfactorily.  This may require updates to the code or simple responses to questions.  Once this is done and the reviewer(s) are satisfied, they will approve the `pull request` and that's a wrap.  Your module is officially included in `nf-core/modules`!  

<br>

## Congratulations, you are officially a contributor to `nf-core`!  

<br>

Anyone wanting to use your module now can be installed directly with

```Bash
nf-core modules install nonpareil
```


# Documentation and Links 
* Join nf-core (https://nf-co.re/join)
* nf-core Documentation (https://nf-co.re/docs/usage/introduction)  
* nf-core Github (https://github.com/nf-core)
    * modules (https://github.com/nf-core/modules)
    * test datasets (https://github.com/nf-core/test-datasets)
* List of Singularity containers available in Galaxy (https://depot.galaxyproject.org/singularity/)
* nf-core Troubleshooting (https://nf-co.re/docs/usage/troubleshooting)


---

<h5><a href="/nextflow_varcal/nextflow/nextflow_nfcore_variantcall" style="float: left"><b>Back to:</b>nf-core/variantcall</a>

<a href="/nextflow_varcal/nextflow/index" style="float: right"><b>Next:</b>Table of Contents</a></h5>
