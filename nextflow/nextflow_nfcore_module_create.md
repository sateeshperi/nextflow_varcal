---
layout: main
title: nf-core module creation
categories: [nextflow]
tags: [cluster, nextflow, workflow, bioinformatics, tutorial]
permalink: /nextflow/nextflow_nfcore_module_create
---

{% include _nextflow_nextflow_nfcore_module_create_toc.html %}

<hr>
<center>This is part of <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Overview

- Initial Environment Setup (`conda`)
- Initialize a new `nf-core` pipeline
- Install or create `nf-core` modules needed for the pipeline
- Create subworkflows for different pipeline steps/phases
- Create a primary workflow
- Modify `config` files for the pipeline components
- Test the pipeline
- Deployment, distribution, and support

<br>

# 1. Initial environment setup (optional; recommended for local setup)

## _(requires conda/mamba to be installed prior to beginning)_

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

## Creating a New `nf-core` Pipeline from Scratch

Here are the steps to create a new `nf-core` pipeline from scratch:

1. **Initial Environment Setup (Optional; Recommended for Local Setup)**

   Before you begin, you need to set up your environment. This requires `conda` or `mamba` to be installed. First, create a Conda environment for Nextflow:

   ```bash
   # Create a new empty environment
   mamba create -n nextflow
   mamba install -c conda-forge -c bioconda nextflow=21.10.6 nf-core=2.2 graphviz openjdk=8.0.312 git=2.35.0

   # Activate the environment
   conda activate nextflow

   # Create a clean, sharable copy of this conda environment
   conda env export | grep -v "prefix" > env.nextflow.yml
   ```

2. **Initialize a New `nf-core` Pipeline**

   Next, you can initialize your new pipeline. This will use templates standardized by `nf-core tools`, which are available from the `nf-core` package in your Conda environment:

   ```bash
   # This will activate an interactive prompt to initialize the pipeline name, author, and git repo
   nf-core create

   # Workflow Name: quaisar
   # Description: nf-core demo workflow
   # Author: hseabolt
   ```

   Then, navigate to your new pipeline's directory:

   ```bash
   cd nf-core-quaisar
   ```

3. **Install or Create `nf-core` Modules Needed for the Pipeline**

   `nf-core` has many modules already available for a variety of bioinformatics software. You can download and install these directly into your new pipeline using `nf-core` tools. Use the following command to check which modules are already created and peer-reviewed:

   ```bash
   # Scroll through this on the terminal or pipe to grep to check a specific module
   nf-core modules list remote

   # Grep example
   nf-core modules list remote | grep "bwa/align"
   ```

4. **Create Subworkflows for Different Pipeline Steps/Phases**

   This is where you would define the series of steps that form the pipeline. Each subworkflow can contain one or more modules.

5. **Create a Primary Workflow**

   The primary workflow ties together all the subworkflows into a cohesive pipeline.

6. **Modify Config Files for the Pipeline Components**

   You can adjust settings for your pipeline by modifying the appropriate configuration files.

7. **Test the Pipeline**

   It's good practice to thoroughly test your pipeline before deploying it.

8. **Deployment, Distribution, and Support**

   Once your pipeline is tested and ready, you can deploy and distribute it. You may also want to provide support for users of your pipeline.

Remember to track your changes with Git:

```bash
git add --all
git commit -m "Initial commit of quaisar nf-core pipeline"

# Prior to pushing to an online repo, make sure to create the repo online and don't add any files to it
git remote add origin <github repo URL>
git push origin master
```

## Installing an Available Module Directly from nf-core

nf-core provides a simple interface for module installations. Most software consists of a single command plus command-line parameter arguments -- these modules can be installed directly. However, some software has multiple sub-commands, such as `bwa index`, `bwa align`, `bwa mem`, etc. In these cases, nf-core requires each sub-command to be created as a standalone module within the parent module, which can also be directly installed. Be aware that installing the parent module (i.e., `bwa` in this case) through nf-core does not give you direct access to its sub-modules.

```bash
# Install nf-core modules from nf-core remote (master) repo
nf-core modules install samtools_stats

# Install modules that have specific sub-modules
nf-core modules install bwa/index
```

## Creating New Modules for nf-core

If you just want to start writing code, here's a brief rundown of the steps:

1. **Join nf-core Slack**: Go to https://nf-co.re/join to join the community.

2. **Fork and Clone the `nf-core/modules` Repo**: Fork the `nf-core/modules` repo (https://github.com/nf-core/modules) to your own GitHub repositories, then clone it locally from your fork.

3. **Create a New Branch**: Create a new branch for the new module within your local clone of `nf-core/modules`.

4. **Raise an Issue on the Main nf-core GitHub**: Alert others that you are working on this module by raising an issue on the main nf-core GitHub (https://github.com/nf-core). Go to `Issues` --> `New Module`. Add yourself to the `Assignees`.

5. **Check for Existing Conda Recipe, Docker, or Singularity Image**: Ensure a Conda recipe, Docker, or Singularity image already exists in Bioconda, Biocontainers, quay.io, or similar for the software you want to create a module for.

6. **Use nf-core Tools to Create, Edit, and Test the Module**: The nf-core tools package provides useful commands for module creation and testing.

7. **Push Your Module's Branch to Your nf-core/modules Fork**: After you're done editing and testing, push your module's branch to your fork of `nf-core/modules`.

8. **Open a Pull Request on the Main nf-core/modules Repo**: Request to merge your changes into the main nf-core/modules repository.

9. **Request a Review in the nf-core Slack**: In the `request-review` channel, ask for someone to review the pull request and approve the merge.

## Forking the nf-core Modules GitHub Repo

Start by navigating to https://github.com/nf-core/modules and click `Fork` in the top right corner of the web browser. This will create a forked copy of the entire `nf-core/modules` repo in your own GitHub repository.

Next, we should let other `nf-core` developers know that we are working on creating new modules. This way, we avoid duplicating efforts or interfering with the work of others. In the main `nf-core/modules` GitHub repo, open an `Issue` (https://github.com/nf-core/modules/issues) with `New Issue` --> `New Module`.

_(Remember to search the open modules to see if someone is already working on the module you intend to add!)_

Read over the `markdown` text in the dialogue box for new modules and fill out/edit appropriately. Assign yourself in the `Assignees` in the top right corner. When you are finished, click `Submit New Issue`.

Now, we need to clone the repo locally so that we can branch/modify it.

```bash
# From a directory where you want to clone the fork of the modules repo from your Git.
# DO NOT do this inside the workflow that you are creating
git clone https://github.com/<your-username>/modules.git
cd modules

# The modules directory should look something like this
ls
# You should see: docs, LICENSE, main.nf, modules...

# Create a new branch for the module you want to create
# Here we are creating "nonpareil", a software used for metagenomics analyses.
git checkout -b nonpareil

# IMPORTANT: Before modifying or creating ANY new files, always make sure that you are not working in the master branch!!!
nf-core modules create nonpareil --author @<your-username> --label process_low --meta

# The new directories for the new module are created in ./modules and ./tests
#    ./modules/nonpareil/main.nf
#    ./modules/nonpareil/meta.yml
#    ./tests/modules/nonpareil/main.nf
#    ./tests/modules/nonpareil/test.yml
#    ./tests/modules/nonpareil/nextflow.config
#    ./tests/config/pytest_modules.yml
```

Now you can start editing these files to add your new module!

### Where Do I Write My Actual Code?

At this point, you are ready to write the actual module code. This involves modifying the `main.nf` file, which includes a number of `TODO` statements that must all be addressed for the module creation to be successful and ready to be merged into `nf-core`.

Your main code for running the software goes in the `./modules/nonpareil/main.nf` file. This includes references to the `Singularity` containers, inputs, outputs, and any other details. It's critical to include the code block at the bottom that captures the `Version` number of the module software (only the version number, not a whole line of info).

If there are any special parameter arguments (e.g., `-t 1 -p 2`, etc.), these should be included in the `./tests/modules/nonpareil/nextflow.config` and **NOT** in the module code in `main.nf`.

Here's an example `main.nf` to run the module:

```groovy
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
    nonpareil \\
        -s $reads \\
        $args \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(nonpareil -V | cut -d ' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
```

In this script, `nonpareil` is the command to run the software, `-s $reads` refers to the input, `$args` contains optional arguments, and `-o ${prefix}` specifies the output prefix. The version of `nonpareil` is then captured and written into a `versions.yml` file.

### The Accompanying `meta.yml` File Describes the New Module in Detail

The `./modules/nonpareil/meta.yml` file acts as a companion to the `main.nf` file, providing documentation that describes the module's purpose, authors, sources, and details about its inputs and outputs. It is important that the input and output specifications in the `meta.yml` file match those defined in the `main.nf` file.

Typically, you can find the information required for the `meta.yml` file from the original GitHub repository or other resources where the source code for the software is hosted.

Here's an example of what a `meta.yml` file might look like for this module (`./modules/nonpareil/meta.yml`)

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
      licence: ["Artistic Licence 2.0"]

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

After you've written the module code `main.nf` and `meta.yml` and you're satisfied with it, it's time to test the module code locally. Here's how you can go about it.

#### Step 1: Identify suitable testing data for your module

`nf-core` provides a large set of test data in various common bioinformatics formats. These datasets are available in the `nf-core/test-datasets` repository on their main GitHub. You can easily use this data for `nf-core` module testing using a meta map, without the need to download any data directly.

You can check the list of available data in the `./modules/tests/config/test_data.config` file.

#### Step 2: Edit the test module

Next, you need to modify the `./tests/modules/nonpareil/main.nf` file to include the test data you identified for your module. Ensure that the data and arguments match what the main code in `./modules/nonpareil/main.nf` expects.

For instance, suppose you want to include a gzipped FASTQ file from the available `nf-core` datasets (e.g., `test_1_fastq_gz`). Your `./modules/tests/modules/nonpareil/main.nf` file might look like this:

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

In this example, we're testing our `NONPAREIL` module with both single-end and paired-end data. The `input` for each workflow is a tuple containing a meta map and one or two input files, depending on the type of data. The `checkIfExists: true` option ensures that the workflow will only run if the specified input file exists. Finally, `NONPAREIL ( input, 50 )` runs our `NONPAREIL` module with the specified input and an additional argument of `50`.

<br>

### Incorporating optional parameters for the module

If you want to incorporate or test any optional parameters for your module, you can specify them in the `./modules/tests/modules/nonpareil/nextflow.config` file. An example of this, using an argument `-T kmer` and setting the `prefix` parameter that is fed into the main code, might look like this:

```Groovy
process {
    publishDir = { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" }

    withName: NONPAREIL {
        ext.args = '-T kmer' // extend this string to include all args, do not use separate strings
        ext.prefix = { "${meta.id}.np" }
    }
}
```

### Testing the module

Now you're ready to test your `nf-core` module. In the root folder of `nf-core/modules`, you can run the following commands:

```bash
# Start to test the module
nf-core modules create-test-yml nonpareil

# This will create an interactive prompt on the terminal.
# For most options, just press ENTER / leave blank.
# For testing profile: use Singularity profile
```

If the testing is successful, a new file will be created: `tests/modules/nonpareil/test.yml`.

You can then run local linting, which will give you more specific details about the tests that were run:

```bash
# Now run local linting
nf-core lint nonpareil --dir .
```

If all tests and linting are successful (or in some cases, safely ignorable), then you're ready to commit the code to git and initiate merge requests:

```bash
# Check the status of the repository
git status

# Stage changes for commit
git add --all  # or instead of --all, add the files listed by git status

# Commit the changes
git commit -m "Created new nf-core module Nonpareil"

# Push the changes to your remote repository
git push -u origin nonpareil
```

If you want to manually check the testing output files for correctness (e.g., FASTA files created by the module), you can `cd` into the `work` directory specified in the terminal output from `create-test-yml` and verify that the output files are what you expect.

<br>

### Merging your local module code with the main `nf-core` repository

To share your module with others, you'll need to create a `Pull Request` (`PR`) across forks from your local copy of `nf-core/modules`.

1. On the `nf-core/modules` Github page, navigate to the `Pull requests` tab, then choose `New pull request`.
2. Fill out the dialogue box for the new pull request. Make sure to get the `Issue Number` from the original `New Module Issue` that you opened previously, and include this to let others know which `Issue` you are closing out.
3. Double check that you are merging `yourusername/modules`->`nonpareil` into `nf-core/modules`->`master`.

When a new `PR` is created, the `nf-core` Github will automatically run some integration tests and checks. You must wait for all of these checks to complete and pass prior to proceeding. If any tests do not pass, you need to address these issues before your new code can be merged. These are often small things like whitespace or minor formatting issues.

As you address these issues, you need to re-run `git add`, `git commit`, `git push` to your local fork, which will automatically update the open `PR` with `nf-core` and re-trigger the integration tests. You do not need to create a new `PR`.

Once all module tests pass, add `Labels` (`New Module`, `Ready_For_Review`) on the right side of the `PR` dashboard in Github.

### Requesting peer review from other `nf-core` developers

After the tests pass, you should request a review from other `nf-core` developers:

1. Copy the `PR`'s URL link and paste it in the `nf-core` Slack channel `#request-review` and politely ask for someone to review your new module.
2. Once a reviewer agrees to review, they will be able to make comments or ask questions through Github that you must address satisfactorily. This may require updates to the code or simple responses to questions.
3. Once your reviewer(s) are satisfied, they will approve the `pull request`. Congratulations! Your module is officially included in `nf-core/modules`!

Anyone wanting to use your module now can install it directly with:

```bash
nf-core modules install nonpareil
```

## Congratulations, you are officially a contributor to `nf-core`!

### Documentation and Links

- [Join nf-core](https://nf-co.re/join)
- [nf-core Documentation](https://nf-co.re/docs/usage/introduction)
- [nf-core Github](https://github.com/nf-core)
  - [modules](https://github.com/nf-core/modules)
  - [test datasets](https://github.com/nf-core/test-datasets)
- [List of Singularity containers available in Galaxy](https://depot.galaxyproject.org/singularity/)
- [nf-core Troubleshooting](https://nf-co.re/docs/usage/troubleshooting)

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_nfcore_variantcall" style="float: left"><b>Back to:</b>nf-core/variantcall</a>

<a href="/nextflow_varcal/nextflow/index" style="float: right"><b>Next:</b>Table of Contents</a></h5>
