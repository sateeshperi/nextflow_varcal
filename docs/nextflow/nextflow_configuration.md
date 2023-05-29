---
layout: main
title: Nextflow Configuration
categories: [nextflow]
tags: [cluster, nextflow, workflow, bioinformatics, tutorial]
permalink: /docs/nextflow/nextflow_configuration
---

{% include _docs_nextflow_nextflow_configuration_toc.html %}

<hr>
<center>This is part 10 of 14 of a <a href="/nextflow_varcal/docs/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

In this tutorial, we will discuss how to optimize Nextflow configuration and provide necessary information about Nextflow syntax.

## Nextflow configuration

A key Nextflow feature is the ability to decouple the workflow implementation, which describes the flow of data and operations to perform on that data, from the configuration settings required by the underlying execution platform. This enables the workflow to be portable, allowing it to run on different computational platforms such as an institutional HPC or cloud infrastructure, without needing to modify the workflow implementation.

We have seen earlier that it is possible to provide a `process` with directives. These directives are process-specific configuration settings. Similarly, we have also provided parameters to our workflow which are parameter configuration settings. These configuration settings can be separated from the workflow implementation, into a configuration file.

### Profiles for Conda, Docker, and Singularity

We have been using an activated conda environment `varcal` to provide tools for our Nextflow scripts. However, we can create profiles that allow you to control the creation of a Conda environment by the Conda package manager and Singularity containers.

Create a `nextflow.config` file with the following content:

```groovy
// nextflow.config
profiles {
  conda {
    process.conda = "/workspace/nextflow_tutorial/environment.yml"
    conda.useMamba = true
  }

  docker {
    docker.enabled = true
  }

  singularity {
    singularity.enabled = true
  }
}
```

To run the `variant_calling.nf` using the conda profile:

```bash
nextflow run variant-calling.nf -profile conda
```

Nextflow will proceed with creating a conda environment based on the `environment.yml` file specified in `nextflow.config`.

## Configuration files

Configuration files contain settings in the form of name-value pairs (`name = value`). The `name` is a specific property to set, while the `value` can be anything you can assign to a variable, such as strings, booleans, or other variables. It is also possible to access any variable defined in the host environment, such as `$PATH`, `$HOME`, `$PWD`, etc.

```groovy
// nextflow.config
my_home_dir = "$HOME"
```

Accessing variables in your configuration file: Generally, variables and functions defined in a configuration file are not accessible from the workflow script. Only variables defined using the `params` scope and the `env` scope (without `env` prefix) can be accessed from the workflow script.

```groovy
workflow {
    MY_PROCESS( params.input )
}
```

Settings are partitioned into scopes, which govern the behavior of different elements of the workflow. For example, workflow parameters are governed from the params scope, while process directives are governed from the process scope. A full list of the available scopes can be found in the [documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes). **It is also possible to define your own scope.**

Configuration settings for a workflow are often stored in the file `nextflow.config`, which is in the same directory as the workflow script. Configuration can be written in either of two ways: using `dot` notation or using `brace` notation. Both forms of notation can be used in the same configuration file.

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

Configuration files can also be separated into multiple files and included into another using the `includeConfig` statement.

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

This modular approach makes it easier to manage and organize your configuration settings across different files, enhancing the readability and maintainability of your workflow.

In summary, Nextflow provides a powerful and flexible configuration system that allows you to decouple the workflow implementation from its execution settings. By utilizing profiles, scopes, and different notations, you can create highly portable and modular workflows that can run on different computational platforms without the need to modify the workflow implementation.

### Combining configuration strategies

To make your workflow portable, modular, and adaptable, it is recommended to use a combination of the above-discussed strategies. Here's an example showing how to organize a Nextflow configuration:

1. Create a base configuration file (e.g., `base.config`) to define general settings for all processes:

```groovy
// base.config
process {
    cpus = 2
    memory = 8.GB
    time = '1 hour'
    publishDir = [ path: params.outdir, mode: 'copy' ]
}
```

2. Create a process-specific configuration file (e.g., `process_resources.config`) to define settings for individual processes:

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

3. Create a configuration file (e.g., `labels.config`) to define settings for processes using labels:

```groovy
// labels.config
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
    }
}
```

4. Include all these configuration files in the main `nextflow.config` file:

```groovy
// nextflow.config
includeConfig 'base.config'
includeConfig 'process_resources.config'
includeConfig 'labels.config'
```

5. Finally, if needed, use profiles for different execution environments (e.g., `conda`, `docker`, or `singularity`) as discussed earlier.

By combining these strategies, you can create a flexible, organized, and easily maintainable configuration system for your Nextflow workflows.

In conclusion, Nextflow's configuration system enables you to create modular and portable workflows that can run on various computational platforms. By utilizing profiles, scopes, different notations, and strategies to manage process behavior, you can efficiently separate the workflow implementation from its execution settings and adapt your workflows to various environments and requirements.

**Selector priority**

When mixing generic process configuration and selectors, the following priority rules are applied (from highest to lowest):

1. `withName` selector definition.
2. `withLabel` selector definition.
3. Process specific directive defined in the workflow script.
4. Process generic `process` configuration.

**Dynamic expressions**

Nextflow allows you to use dynamic expressions for certain configuration settings that depend on the data being processed. You can achieve this by using closures. For example, you can specify the required `memory` as a multiple of the number of `cpus`. Similarly, you can publish results to a subfolder based on the sample name.

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

Nextflow supports a wide range of execution platforms, from running locally to running on HPC clusters or cloud infrastructures. You can find the full list of supported executors in the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

The default executor configuration is defined within the [`executor` scope](https://www.nextflow.io/docs/latest/config.html#scope-executor). For example, in the config below, we specify the executor as `slurm` and the number of tasks the executor will handle in a parallel manner (`queueSize`) to 2.

```groovy
// nextflow.config
executor {
    name = 'slurm'
    queueSize = 2
}
```

The `process.executor` directive allows you to override the executor to be used by a specific process. This can be useful, for example, when there are short-running tasks that can be run locally and are unsuitable for submission to HPC executors (check for guidelines on the best practice use of your execution system). Other process directives such as `process.clusterOptions`, `process.queue`, and `process.machineType` can also be used to further configure processes depending on the executor used.

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

By understanding Nextflow's syntax and configuration options, you can create flexible, modular workflows that run efficiently on various computational platforms. Utilize selector priority, dynamic expressions, and platform configurations to optimize your workflow and adapt it to different environments and requirements.

## Configuring software requirements

Nextflow supports various technologies for managing software, including Conda, Docker, Singularity, Podman, Charliecloud, and Shifter. These technologies enable you to package tools and their dependencies into a software environment, ensuring that the tools work as long as the environment is loaded. This promotes portable and reproducible workflows. Software environment specification is managed from the process scope, allowing the use of `process` selectors to manage which processes load which software environment. Each technology also has its own scope for providing further technology-specific configuration settings.

### Software configuration using Conda

Conda is a software package and environment management system that runs on Linux, Windows, and Mac OS. Conda packages software into environments along with their dependencies for a specific operating system. Conda environments can be configured in several ways:

- Provide a path to an existing Conda environment.
- Provide a path to a Conda environment specification file (written in YAML).
- Specify the software package(s) using the `<channel>::<package_name>=<version>` syntax (separated by spaces), which builds the Conda environment when the process is run.

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

The optional `conda` scope allows you to control the creation of a Conda environment by the Conda package manager. For example, `conda.cacheDir` specifies the path where the Conda environments are stored, which is in the `conda` folder of the `work` directory by default.

## Software configuration using Docker

Docker is a container technology that allows you to package software into lightweight, standalone, executable containers. Docker containers run software consistently, regardless of the underlying infrastructure. To use Docker, you must provide a container image path using the `process.container` directive and enable Docker in the Docker scope with `docker.enabled = true`. You should also supply your user and group via `docker.runOptions`.

```groovy
process.container = 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'
```

## Software configuration using Singularity

Singularity is another container technology often used on HPC clusters. Singularity runs processes as the user and automatically mounts certain directories in the container instance. Singularity also supports building Singularity images from Docker images, allowing Docker image paths to be used as values for `process.container`.

To enable Singularity, provide a container image path using `process.container` and enable Singularity using `singularity.enabled = true`.

```groovy
process.container = 'https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0'
singularity.enabled = true
```

> Container protocols - The following protocols are supported:
>
> - `docker://`: download the container image from the Docker Hub and convert it to the Singularity format (default).
> - `library://`: download the container image from the Singularity Library service.
> - `shub://`: download the container image from the Singularity Hub.
> - `docker-daemon://`: pull the container image from a local Docker installation and convert it to a Singularity image file.
> - `https://`: download the Singularity image from the given URL.
> - `file://`: use a Singularity image on local computer storage.

By leveraging the power of Conda, Docker, and Singularity, you can create portable and reproducible workflows that run efficiently on various computational platforms.

## Software configuration using Podman, Charliecloud, and Shifter

These are other container technologies supported by Nextflow. Their usage is similar to Docker and Singularity, and they provide additional options for running containerized applications.

- **Podman** is a daemonless container engine for developing, managing, and running OCI Containers on your Linux System. It provides a Docker-compatible command line that eases the transition from other container engines.

- **Charliecloud** simplifies the deployment of containers in high-performance computing (HPC) environments. It uses Linux user namespaces to run containers with no privileged operations or daemons and minimal configuration changes on center resources.

- **Shifter** is a container solution for HPC centers. It is designed to be performant and integrates well with resource managers and MPI.

Configuration of these technologies in Nextflow follows a similar pattern to Docker and Singularity. For example, you would provide a container image path using `process.container` and enable the corresponding technology in its scope.

```groovy
process.container = 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
podman.enabled = true

process.container = 'docker://ubuntu:20.04'
charliecloud.enabled = true

process.container = 'docker://ubuntu:20.04'
shifter.enabled = true
```

Remember to consult the official documentation of these technologies and Nextflow for more detailed instructions and best practices. It's also important to know the specific requirements and limitations of your HPC center or cloud provider, as not all container technologies may be supported or recommended.

Nextflow's flexibility and compatibility with a wide range of technologies make it a powerful tool for creating and running complex computational pipelines, both locally and in distributed environments such as HPC clusters or cloud infrastructures.

---

<h5><a href="/nextflow_varcal/docs/nextflow/nextflow_variant_calling" style="float: left"><b>Back to:</b>Nextflow Variant Calling</a>

<a href="/nextflow_varcal/docs/nextflow/nextflow_modules" style="float: right"><b>Next:</b>NextFlow Modules</a></h5>
