---
layout: main
title: NF-Core @ HPC
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_slurm
---
{% include _nextflow_nextflow_slurm_toc.html %}


<hr>
<center>This is part 3 of 14 of a <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Executors

* A key feature of Nextflow is the ability to decouple the workflow implementation from the actual execution platform.
* This allows deployment of a workflow on any execution platform supported by Nextflow.

<img src="images/nf-executors.PNG" alt="drawing" width="500"/>

<br>

* The executor defines the underlying system where processes are executed. By default, a process uses the executor defined globally in the `nextflow.config` file.

* The `process.executor` directive allows you to configure which executor should be used by the process, overriding the default configuration.

> See the full list of executors [here](https://www.nextflow.io/docs/latest/executor.html).

<br>

## SLURM custom config

* Log in to the HPC cluster and request an interactive node.

* Navigate to your `nextflow_tutorial` directory and make sure you have downloaded the nf-core/viralintegration pipeline (`nf-core download nf-core/viralintegration -r 0.1.0`)

```bash
cd /workspace/nextflow_tutorial/nf-core-viralintegration-0.1.0/workflow/
```

* Let's create a `nfcore_custom_config` profile to submit jobs to an HPC cluster that uses the `slurm` scheduler.

```bash
mkdir conf/nf-core-slurm
nano conf/nf-core-slurm/nfcore_custom.config
```

* Copy and paste the following into the `nfcore_custom.config` file and save on exit.

```bash
/*
========================================================================================
    NF-CORE Custom Config File
========================================================================================
    Default config options for HPC compute environments
----------------------------------------------------------------------------------------
*/

//Profile config names for nf-core/configs

params {

  config_profile_description = ''
  config_profile_contact     = ''
  config_profile_url         = ''

  // Output options
  outdir                     = "results"
}

/*
========================================================================================
    Nextflow Metrics & Reports
========================================================================================
*/

timeline {
  enabled = true
  file    = "${params.outdir}/timeline.html"
}

report {
  enabled = true
  file    = "${params.outdir}/report.html"
}
trace {
  enabled = true
  fields  = 'task_id,name,status,exit,realtime,%cpu,%mem,rss,vmem,peak_rss,peak_vmem,rchar,wchar'
  file    = "${params.outdir}/trace.txt"
}

/*
========================================================================================
    Base Executor config
========================================================================================
*/

executor {
  queueSize = 2
}

/*
========================================================================================
    Profiles - slurm,singularity,conda,docker
========================================================================================
*/

profiles {
  slurm {
    process {
      executor     = 'slurm'
      queue        = ''
    }
    executor {
      queueSize    = 100
      pollInterval = '15 sec'
    }
  }

  singularity {
    singularity.enabled = true
  }
  
  conda {
    conda.enabled = true
  }
  
  docker {
    docker.enabled = true
  }

}
```

* Now run the pipeline with the following options to submit the Nextflow job to the cluster.

```bash
nextflow run main.nf --custom_config_base conf/nf-core-slurm -profile test,docker --outdir results
```

* You should see output similar to below in the directory `results-slurm`.
* Notice the change of `executor >  slurm`. Nextflow will now handle the submission and execution of the rest of the jobs on the cluster.

```bash

```

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_nfcore" style="float: left"><b>Back to:</b>NF-Core</a>

<a href="/nextflow_varcal/nextflow/nextflow_scripting" style="float: right"><b>Next:</b>Nextflow Scripting</a></h5>
