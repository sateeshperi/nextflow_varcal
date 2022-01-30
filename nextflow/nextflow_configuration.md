---
layout: main
title: Nextflow Configuration
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_configuration
---
{% include _nextflow_nextflow_configuration_toc.html %}


<hr>
<center>This is part 10 of 14 of a <a href="/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>


## Nextflow configuration

*   A key Nextflow feature is the ability to decouple the workflow implementation, which describes the flow of data and operations to perform on that data, from the configuration settings required by the underlying execution platform.
*   This enables the workflow to be portable, allowing it to run on different computational platforms such as an institutional HPC or cloud infrastructure, without needing to modify the workflow implementation.
*   We have seen earlier that it is possible to provide a `process` with directives. These directives are process specific configuration settings.
*   Similarly, we have also provided parameters to our workflow which are parameter configuration settings. These configuration settings can be separated from the workflow implementation, into a configuration file.


We have been using an activated conda environment `varcal` to provide tools for our nextflow scripts. However we can create profiles which allows you to control the creation of a Conda environment by the Conda package manager and Singularity containers.

```bash
cd ~/nextflow_tutorial
touch nextflow.config
```

Paste the following

```groovy
//nextflow.config
profiles {
  conda {
    process.conda = "$HOME/nextflow_tutorial/environment.yml"
  }
  singularity {
    singularity.enabled = true
  }
}
```

Now make sure you have deactivated the conda environment using `conda deactivate` 

To run the variant_calling.nf using conda

```bash
nextflow run variant-calling.nf -profile conda 
```

>**Nextflow will proceed with creating a conda environment based on the `environment.yml` file in `nextflow.config`**
>```bash
>N E X T F L O W  ~  version 21.04.3
>Launching `variant-calling.nf` [backstabbing_yonath] - revision: 119a92be39
>V A R I A N T-C A L L I N G - N F   P I P E L I N E
>===================================
>genome       : nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta
>reads        : nextflow_tutorial/data/trimmed_fastq/SRR2584863_{1,2}.trim.fastq.gz
>outdir       : results
>
>[-        ] process > FASTQC           -
>[-        ] process > FASTQC           -
>[-        ] process > BWA_INDEX        -
>[-        ] process > BWA_ALIGN        -
>[-        ] process > SAMTOOLS_SORT    -
>[-        ] process > SAMTOOLS_INDEX   -
>[-        ] process > BCFTOOLS_MPILEUP -
>[-        ] process > BCFTOOLS_CALL    -
>[-        ] process > VCFUTILS         -
>Creating Conda env: nextflow_tutorial/environment.yml [cache nextflow_tutorial/work/conda/environment-bbe607cdc5e24ff213d1d83fc7d52ef1]
>```

## Configuration files

*   Settings in a configuration file are sets of name-value pairs (`name = value`).
*   The `name` is a specific property to set, while the `value` can be anything you can assign to a variable (see [nextflow scripting](/nextflow/nextflow_scripting)), for example, strings, booleans, or other variables.
*   It is also possible to access any variable defined in the host environment such as `$PATH`, `$HOME`, `$PWD`, etc.

>```groovy
>// nextflow.config
>my_home_dir = "$HOME"
>```

> Accessing variables in your configuration file - Generally, variables and functions defined in a configuration file are not accessible from the workflow script. Only variables defined using the `params` scope and the `env` scope (without `env` prefix) can be accessed from the workflow script.
>```groovy
>workflow {
>    MY_PROCESS( params.input )
>}
>```

*   Settings are also partitioned into scopes, which govern the behaviour of different elements of the workflow. For example, workflow parameters are governed from the params scope, while process directives are governed from the process scope. A full list of the available scopes can be found in the [documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes). **It is also possible to define your own scope.**

*   Configuration settings for a workflow are often stored in the file `nextflow.config` which is in the same directory as the workflow script. Configuration can be written in either of two ways. The first is using `dot` notation, and the second is using `brace` notation. Both forms of notation can be used in the same configuration file.

An example of dot notation:

```groovy
params.input = ''             // The workflow parameter "input" is assigned an empty string to use as a default value
params.outdir = './results'   // The workflow parameter "outdir" is assigned the value './results' to use by default.
```

An example of brace notation:

```groovy
params {
    input  = ''
    outdir = './results'
}
```

*   Configuration files can also be separated into multiple files and included into another using the includeConfig statement.

```groovy
// nextflow.config
params {
    input  = ''
    outdir = './results'
}

includeConfig 'system_resources.config'
```

```groovy
// system_resources.config
process {
    cpus = 1    // default cpu usage
    time = '1h' // default time limit
}
```

### How configuration files are combined

*   Configuration settings can be spread across several files. This also allows settings to be overridden by other configuration files. The priority of a setting is determined by the following order, ranked from highest to lowest.

1.  Parameters specified on the command line (`--param_name value`).
2.  Parameters provided using the `-params-file` option.
3.  Config file specified using the `-c` my_config option.
4.  The config file named `nextflow.config` in the current directory.
5.  The config file named `nextflow.config` in the workflow project directory (`$projectDir`: the directory where the script to be run is located).
6.  The config file `$HOME/.nextflow/config`.
7.  Values defined within the workflow script itself (e.g., `main.nf`).

> If configuration is provided by more than one of these methods, configuration is merged giving higher priority to configuration provided higher in the list.

*   Existing configuration can be completely ignored by using `-C <custom.config>` to use only configuration provided in the `custom.config` file.

> Configuring Nextflow vs Configuring a Nextflow workflow 
> Parameters starting with a single dash `-` (e.g., `-c my_config.config`) are configuration options for nextflow, while parameters starting with a double dash `--` (e.g., `--outdir`) are workflow parameters defined in the `params` scope.
> The majority of Nextflow configuration settings must be provided on the command-line, however a handful of settings can also be provided within a configuration file, such as `workdir = '/path/to/work/dir'` (`-w /path/to/work/dir`) or `resume = true` (`-resume`), and do not belong to a configuration scope.

## Configuring process behaviour

*   Earlier we saw that `process` directives allow the specification of settings for the task execution such as `cpus`, `memory`, `conda` and other resources in the pipeline script. This is useful when prototyping a small workflow script, however this ties the configuration to the workflow, making it less portable. A good practice is to separate the process configuration settings into another file.

The `process` configuration scope allows the setting of any process directives in the Nextflow configuration file.

For example:

```groovy
// nextflow.config
process {
    cpus = 2
    memory = 8.GB
    time = '1 hour'
    publishDir = [ path: params.outdir, mode: 'copy' ]
}
```

## Process selectors

*   The resources for a specific process can be defined using `withName:` followed by the process name ( either the simple name e.g., `'FASTQC'`, or the fully qualified name e.g., `'NFCORE_RNASEQ:RNA_SEQ:SAMTOOLS_SORT'`), and the directives within curly braces. For example, we can specify different `cpus` and `memory` resources for the processes `INDEX` and `FASTQC` as follows:

```groovy
// process_resources.config
process {
    withName: INDEX {
        cpus = 4
        memory = 8.GB
    }
    withName: FASTQC {
        cpus = 2
        memory = 4.GB
    }
}
```

*   When a workflow has many processes, it is inconvenient to specify directives for all processes individually, especially if directives are repeated for groups of processes. A helpful strategy is to annotate the processes using the `label` directive (processes can have multiple labels). The `withLabel` selector then allows the configuration of all processes annotated with a specific label, as shown below:

```groovy
// configuration_process_labels.nf
nextflow.enable.dsl=2

process P1 {

    label "big_mem"

    script:
    """
    echo P1: Using $task.cpus cpus and $task.memory memory.
    """
}

process P2 {

    label "big_mem"

    script:
    """
    echo P2: Using $task.cpus cpus and $task.memory memory.
    """
}

workflow {

    P1()
    P2()

}
```

```groovy
// configuration_process-labels.config
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
    }
}
```

Another strategy is to use process selector expressions. Both `withName:` and `withLabel:` allow the use of regular expressions to apply the same configuration to all processes matching a pattern. Regular expressions must be quoted, unlike simple process names or labels.

* The `|` matches either-or, e.g., `withName: 'INDEX|FASTQC'` applies the configuration to any process matching the name `INDEX` or `FASTQC`.
* The `!` inverts a selector, e.g.,` withLabel: '!small_mem'` applies the configuration to any process without the `small_mem label`.
* The `.*` matches any number of characters, e.g., `withName: 'NFCORE_RNASEQ:RNA_SEQ:BAM_SORT:.*'` matches all processes of the workflow `NFCORE_RNASEQ:RNA_SEQ:BAM_SORT`.
* 
A regular expression cheat-sheet can be found [here](https://www.jrebel.com/system/files/regular-expressions-cheat-sheet.pdf) if you would like to write more expressive expressions.

**Selector priority**

When mixing generic process configuration and selectors, the following priority rules are applied (from highest to lowest):

1.  `withName` selector definition.
2.  `withLabel` selector definition.
3.  Process specific directive defined in the workflow script.
4.  Process generic `process` configuration.

**Dynamic expressions**

*   A common scenario is that configuration settings may depend on the data being processed. Such settings can be dynamically expressed using a closure. For example, we can specify the `memory` required as a multiple of the number of `cpus`. Similarly, we can publish results to a subfolder based on the sample name.

```groovy
process FASTQC {

    input:
    tuple val(sample), path(reads)

    script:
    """
    fastqc -t $task.cpus $reads
    """
}
```

```groovy
// nextflow.config
process {
    withName: FASTQC {
        cpus = 2
        memory = { 2.GB * task.cpus }
        publishDir = { "fastqc/$sample" }
    }
}
```

## Configuring execution platforms

*   Nextflow supports a wide range of execution platforms, from running locally, to running on HPC clusters or cloud infrastructures. See https://www.nextflow.io/docs/latest/executor.html for the full list of supported executors.

![](images/nf-executors.PNG)

*   The default executor configuration is defined within the [`executor` scope](https://www.nextflow.io/docs/latest/config.html#scope-executor). For example, in the config below we specify the executor as `slurm` and the number of tasks the executor will handle in a parallel manner (`queueSize`) to 2.

```groovy
// nextflow.config
executor {
    name = 'slurm'
    queueSize = 2
}
```

The `process.executor` directive allows you to override the executor to be used by a specific process. This can be useful, for example, when there are short running tasks that can be run locally, and are unsuitable for submission to HPC executors (check for guidelines on best practice use of your execution system). Other process directives such as `process.clusterOptions`, `process.queue`, and `process.machineType` can be also be used to further configure processes depending on the executor used.

```groovy
//nextflow.config
executor {
    name = 'slurm'
    queueSize = 2
}
process {
    withLabel: 'short' {
        executor = 'local'
    }
}
```

## Configuring software requirements

*   An important feature of Nextflow is the ability to manage software using different technologies. It supports the Conda package management system, and container engines such as Docker, Singularity, Podman, Charliecloud, and Shifter.
*   These technologies allow one to package tools and their dependencies into a software environment such that the tools will always work as long as the environment can be loaded. This facilitates portable and reproducible workflows. Software environment specification is managed from the process scope, allowing the use of `process`selectors to manage which processes load which software environment. Each technology also has its own scope to provide further technology specific configuration settings.

### Software configuration using Conda

*   Conda is a software package and environment management system that runs on Linux, Windows, and Mac OS. Software packages are bundled into Conda environments along with their dependencies for a particular operating system (Not all software is supported on all operating systems). Software packages are tied to conda channels, for example, bioinformatic software packages are found and installed from the BioConda channel.

A Conda environment can be configured in several ways:

* Provide a path to an existing Conda environment.
* Provide a path to a Conda environment specification file (written in YAML).
* Specify the software package(s) using the `<channel>::<package_name>=<version>` syntax (separated by spaces), which then builds the Conda environment when the process is run.

```groovy
process {
    conda = "/home/user/miniconda3/envs/my_conda_env"
    withName: FASTQC {
        conda = "environment.yml"
    }
    withName: SALMON {
        conda = "bioconda::salmon=1.5.2"
    }
}
```

There is an optional `conda` scope which allows you to control the creation of a Conda environment by the Conda package manager. For example, `conda.cacheDir` specifies the path where the Conda environments are stored. By default this is in `conda` folder of the `work` directory.

## Software configuration using Docker

*   Docker is a container technology. Container images are lightweight, standalone, executable package of software that includes everything needed to run an application: code, runtime, system tools, system libraries and settings. Containerized software is intended to run the same regardless of the underlying infrastructure, unlike other package management technologies which are operating system dependant (See the published article on Nextflow). For each container image used, Nextflow uses Docker to spawn an independent and isolated container instance for each process task.

*   To use Docker, we must provide a container image path using the `process.container` directive, and also enable docker in the docker scope, `docker.enabled = true`. A container image path takes the form `(protocol://)registry/repository/image:version--build`. By default, Docker containers run software using a privileged user. This can cause issues, and so it is also a good idea to supply your user and group via the `docker.runOptions`.

```groovy
process.container = 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'
```

## Software configuration using Singularity

*   Singularity is another container technology, commonly used on HPC clusters. It is different to Docker in several ways. The primary differences are that processes are run as the user, and certain directories are automatically “mounted” (made available) in the container instance. Singularity also supports building Singularity images from Docker images, allowing Docker image paths to be used as values for `process.container`.

*   Singularity is enabled in a similar manner to Docker. A container image path must be provided using `process.container` and singularity enabled using `singularity.enabled = true`.

```groovy
process.container = 'https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0'
singularity.enabled = true
```

>Container protocols - The following protocols are supported:
>
>*  `docker://`: download the container image from the Docker Hub and convert it to the Singularity format (default).
>*  `library://`: download the container image from the Singularity Library service.
>*  `shub://`: download the container image from the Singularity Hub.
>*  `docker-daemon://`: pull the container image from a local Docker installation and convert it to a Singularity image file.
>*  `https://`: download the singularity image from the given URL.
>*  `file://`: use a singularity image on local computer storage.


### Configuration profiles

*   One of the most powerful features of Nextflow configuration is to predefine multiple configurations or `profiles` for different execution platforms. This allows a group of predefined settings to be called with a short invocation, `-profile <profile name>`.

*   Configuration profiles are defined in the `profiles` scope, which group the attributes that belong to the same profile using a common prefix.

```groovy
//configuration_profiles.config
profiles {

    standard {
        params.genome = '/local/path/ref.fasta'
        process.executor = 'local'
    }

    cluster {
        params.genome = '/data/stared/ref.fasta'
        process.executor = 'slurm'
        process.queue = 'long'
        process.memory = '10GB'
        process.conda = '/some/path/env.yml'
    }

    cloud {
        params.genome = '/data/stared/ref.fasta'
        process.executor = 'awsbatch'
        process.container = 'cbcrg/imagex'
        docker.enabled = true
    }

}
```

*   This configuration defines three different profiles: `standard`, `cluster` and `cloud` that set different process configuration strategies depending on the target runtime platform. By convention the standard profile is implicitly used when no other profile is specified by the user. To enable a specific profile use `-profile` option followed by the profile name:

```bash
nextflow run <your script> -profile cluster
```

> Configuration order - Settings from profiles will override general settings in the configuration file. However, it is also important to remember that configuration is evaluated in the order it is read in. For example, in the following example, the `publishDir` directive will always take the value ‘results’ even when the profile `hpc` is used. This is because the setting is evaluated before Nextflow knows about the hpc profile. If the `publishDir` directive is moved to after the profiles scope, then `publishDir` will use the correct value of `params.results`.

```groovy
params.results = 'results'
process.publishDir = params.results
profiles {
    hpc {
        params.results = '/long/term/storage/results'
    }
}
```

## Inspecting the Nextflow configuration

*   You can use the command nextflow config to print the resolved configuration of a workflow. This allows you to see what settings Nextflow will use to run a workflow.

```bash
nextflow config variant-calling.nf -profile singularity
```

>Output
>```bash
>singularity {
>   enabled = true
>}
>```


## Quick Recap

>*  Nextflow configuration can be managed using a Nextflow configuration file.
>*  Nextflow configuration files are plain text files containing a set of properties.
>*  You can define process specific settings, such as cpus and memory, within the `process` scope.
>*  You can assign different resources to different processes using the process selectors `withName` or `withLabel`.
You can define a profile for different configurations using the profiles scope. These profiles can be selected when launching a pipeline >* execution by using the `-profile` command-line option
>*  Nextflow configuration settings are evaluated in the order they are read-in.
>*  Workflow configuration settings can be inspected using `nextflow config <script> [options]`.

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_variant_calling" style="float: left"><b>Back to:</b>Nextflow Variant Calling</a>

<a href="/nextflow_varcal/nextflow/nextflow_modules" style="float: right"><b>Next:</b>NextFlow Modules</a></h5>
