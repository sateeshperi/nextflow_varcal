---
layout: main
title: NF-Core @ HPC
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow_varcal/nextflow/nextflow_slurm
---
{% include _nextflow_nextflow_slurm_toc.html %}


<hr>
<center>This is part 3 of 14 of a <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Executors

* A key Nextflow feature is the ability to decouple the workflow implementation from the actual execution platform.
* This allows the deployment of a workflow on any executing platform support by the Nextflow.

<img src="images/nf-executors.PNG" alt="drawing" width="500"/>

<br>

* The executor defines the underlying system where processes are executed. By default a process uses the executor defined globally in the `nextflow.config` file.

* The process.executor directive allows you to configure what executor has to be used by the process, overriding the default configuration.

> See full list of executors [here](https://www.nextflow.io/docs/latest/executor.html)

<br>

## SLURM custom config

* Login to the HPC cluster and request an interactive node 


* Navigate to your `nextflow_tutorial` directory and make sure you have downloaded the nf-core/bactmap pipeline (`nf-core download nf-core/bactmap -r 1.0.0`)

```bash
cd nextflow_tutorial/nf-core-bactmap-1.0.0/workflow/
```

* Lets create a `nfcore_custom_config` profile to submit jobs to HPC cluster which uses `slurm` scheduler. 

```bash
mkdir conf/nf-core-slurm
nano conf/nf-core-slurm/nfcore_custom.config
```

* Copy paste the following into the `nfcore_custom.config` file and save on exit.

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
    Profiles - slurm,singularity,conda
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

}
```

* Now run the pipeline with following options to submit the nextflow job to the cluster.

```bash
nextflow run main.nf --custom_config_base conf/nf-core-slurm -profile test,singularity,slurm --outdir results-slurm
```

* You should output similar to below in the directory `results-slurm`.
* Notice the change of `executor >  slurm`. Nextflow will now further handle the submission and execution of the rest of the jobs to the cluster. 

```bash
N E X T F L O W  ~  version 21.04.3
Launching `main.nf` [focused_caravaggio] - revision: 2e456c2488


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/bactmap v1.0.0
------------------------------------------------------
Core Nextflow options
  runName                   : focused_caravaggio
  containerEngine           : singularity
  launchDir                 : nf-core-bactmap-1.0.0/workflow
  workDir                   : nf-core-bactmap-1.0.0/workflow/work
  projectDir                : nf-core-bactmap-1.0.0/workflow
  userName                  : sateesh
  profile                   : singularity,test,slurm
  configFiles               : nf-core-bactmap-1.0.0/workflow/nextflow.config

Input/output options
  input                     : https://raw.githubusercontent.com/nf-core/test-datasets/bactmap/samplesheet.csv
  outdir                    : results-slurm

Compulsory parameters
  reference                 : https://raw.githubusercontent.com/nf-core/test-datasets/bactmap/genome/NCTC13799.fna

Optional pipeline steps
  adapter_file              : nf-core-bactmap-1.0.0/workflow/assets/adapters.fas
  remove_recombination      : true
  rapidnj                   : true
  fasttree                  : true
  iqtree                    : true
  raxmlng                   : true

Max job request options
  max_cpus                  : 2
  max_memory                : 2.GB
  max_time                  : 1.h

Institutional config options
  custom_config_base        : conf/nf-core-slurm
  config_profile_description: 
  config_profile_contact    : 
  config_profile_url        : 

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
If you use nf-core/bactmap for your analysis please cite:

* The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
  https://github.com/nf-core/bactmap/blob/master/CITATIONS.md
------------------------------------------------------
executor >  slurm (47)
[4a/1c1cd0] process > NFCORE_BACTMAP:BACTMAP:INPUT_CHECK:SAMPLESHEET_CHECK (samplesheet.csv)                     [100%] 1 of 1 ✔
[72/05a811] process > NFCORE_BACTMAP:BACTMAP:BWA_INDEX (NCTC13799.fna)                                           [100%] 1 of 1 ✔
[18/2944ae] process > NFCORE_BACTMAP:BACTMAP:FASTP (ERR2172267)                                                  [100%] 3 of 3 ✔
[13/16d959] process > NFCORE_BACTMAP:BACTMAP:SUB_SAMPLING:MASH_SKETCH (ERR2172267)                               [100%] 3 of 3 ✔
[95/b7f3e1] process > NFCORE_BACTMAP:BACTMAP:SUB_SAMPLING:RASUSA (ERR2172267)                                    [100%] 3 of 3 ✔
[6c/3b0d9d] process > NFCORE_BACTMAP:BACTMAP:BWA_MEM (ERR2172267)                                                [100%] 3 of 3 ✔
[92/ba4225] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:SAMTOOLS_SORT (ERR2172267)                        [100%] 3 of 3 ✔
[52/da50f3] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:SAMTOOLS_INDEX (ERR2172267)                       [100%] 3 of 3 ✔
[31/173b95] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS (ERR2172267)    [100%] 3 of 3 ✔
[94/42ae0e] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT (ERR2172267) [100%] 3 of 3 ✔
[89/8ae1cf] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS (ERR2172267) [100%] 3 of 3 ✔
[5c/14c1a8] process > NFCORE_BACTMAP:BACTMAP:VARIANTS_BCFTOOLS:BCFTOOLS_MPILEUP (ERR2172267)                     [100%] 3 of 3 ✔
[4a/70d5df] process > NFCORE_BACTMAP:BACTMAP:VARIANTS_BCFTOOLS:BCFTOOLS_FILTER (ERR2172267)                      [100%] 3 of 3 ✔
[23/3a9c2b] process > NFCORE_BACTMAP:BACTMAP:VCF2PSEUDOGENOME (ERR2172267)                                       [100%] 3 of 3 ✔
[fe/27818e] process > NFCORE_BACTMAP:BACTMAP:ALIGNPSEUDOGENOMES                                                  [100%] 1 of 1 ✔
[c9/6920fa] process > NFCORE_BACTMAP:BACTMAP:SNPSITES                                                            [100%] 1 of 1 ✔
[c1/7a8d34] process > NFCORE_BACTMAP:BACTMAP:GUBBINS                                                             [100%] 1 of 1 ✔
[8c/2db914] process > NFCORE_BACTMAP:BACTMAP:CREATE_PHYLOGENY:RAPIDNJ                                            [100%] 1 of 1 ✔
[04/06194d] process > NFCORE_BACTMAP:BACTMAP:CREATE_PHYLOGENY:FASTTREE                                           [100%] 1 of 1 ✔
[73/b32064] process > NFCORE_BACTMAP:BACTMAP:CREATE_PHYLOGENY:IQTREE                                             [100%] 1 of 1 ✔
[1c/33b29c] process > NFCORE_BACTMAP:BACTMAP:CREATE_PHYLOGENY:RAXMLNG                                            [100%] 1 of 1 ✔
[41/c3698f] process > NFCORE_BACTMAP:BACTMAP:GET_SOFTWARE_VERSIONS                                               [100%] 1 of 1 ✔
[12/ee3d76] process > NFCORE_BACTMAP:BACTMAP:MULTIQC (1)                                                         [100%] 1 of 1 ✔
-[nf-core/bactmap] Pipeline completed successfully-
Completed at: 20-Oct-2021 16:43:22
Duration    : 9m 8s
CPU hours   : 0.2
Succeeded   : 47
```

* To use conda 

```bash
nextflow run main.nf --custom_config_base conf/nf-core-slurm -profile conda,test --outdir results-slurm
```

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_nfcore" style="float: left"><b>Back to:</b>NF-Core</a>

<a href="/nextflow_varcal/nextflow/nextflow_scripting" style="float: right"><b>Next:</b>Nextflow Scripting</a></h5>
