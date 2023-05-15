---
layout: main
title: NextFlow Modules
categories: [nextflow]
tags: [cluster, nextflow, workflow, bioinformatics, tutorial]
permalink: /nextflow/nextflow_modules
---

{% include _nextflow_nextflow_modules_toc.html %}

<hr>
<center>This is part 11 of 14 of a <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Modules in Nextflow DSL2

Modules in Nextflow DSL2 offer an extensive range of features that can be utilized according to the granularity of your pipeline creation needs.

> **MODULE** - A module in Nextflow is essentially a `process` that can be reused in various pipelines. It is designed to be as atomic as possible, meaning it cannot be further split into another module. For instance, a module file containing the process definition for a single tool, such as `FastQC`, is a prime example of this concept. Currently, this repository is intended to host only atomic module files. These should be added to the [modules/](https://github.com/nf-core/modules/tree/master/modules) directory of `nf-core/modules`, accompanied by the necessary documentation and tests.

### Key Points About Nextflow (DSL2) Modules

- Similar to most programming languages, Nextflow also supports the concept of creating reusable code blocks or modules.

- Nextflow (DSL2) introduces the provision to define module scripts, which can be included and shared across different workflow pipelines. This increases code reusability and enhances maintainability.

- A module file in Nextflow is simply a script containing one or more process definitions. These definitions can be imported into another Nextflow script, facilitating a modular and organized coding structure.

- A module can consist of various definitions including functions, processes, and workflows. This flexibility allows the creation of complex, yet manageable, pipeline structures.

### Syntax for Module Creation

```nextflow
process MODULE_NAME {
    // Declare inputs and outputs
    input:
    file input_data

    output:
    file 'output_data'

    // Declare the script to be executed
    script:
    """
    your_command ${input_data} > output_data
    """
}
```

In this example, `MODULE_NAME` is the name of the module, `input_data` is an input file, and `your_command` is the command to be run in this module.

This process block can be saved as a Nextflow script file and can be imported into other scripts as a module.

### Syntax for Module Import

```nextflow
include 'path/to/module' as MODULE_ALIAS
```

Here, `path/to/module` is the path to the Nextflow script file that contains the module, and `MODULE_ALIAS` is the alias name you want to use for this module in your script.

This `include` statement allows you to reuse the module across various workflow pipelines.

```bash
mkdir modules
cd modules
```

- Create a new file `fastqc.nf` in the current modules directory; paste the following and save.

> ```bash
> /*
> ========================================================================================
>    FASTQC module
> ========================================================================================
>    Website: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
> ========================================================================================
> */
>
> // Parameter definitions
> params.CONTAINER = "quay.io/biocontainers/fastqc:0.11.9--0"
> params.OUTPUT = "trim_fastqc_output"
>
> process FASTQC {
>
>    // where to store the results and in which way
>    publishDir( params.OUTPUT, mode: 'copy' )
>
>    // indicates to use as container the value indicated in the parameter
>    container( params.CONTAINER )
>
>    // show in the log which input file is analysed
>    tag( "${reads}" )
>
>    input:
>    tuple val( sample_id ), path( reads )
>
>    output:
>    path( "*_fastqc*" ), emit: fastqc_out
>
>    script:
>    """
>    fastqc ${reads}
>    """
> }
> ```

- Create a new file `bwa_index.nf` in the current modules directory; paste the following and save.

> ```bash
> /*
> ========================================================================================
>    BWA-INDEX module
> ========================================================================================
>    Website: http://www.htslib.org/doc/samtools.html
> ========================================================================================
> */
>
> // Parameter definitions
> params.CONTAINER = "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
> params.OUTPUT = "bwa_index"
>
> process BWA_INDEX {
>
>  // where to store the results and in which way
>  publishDir( params.OUTPUT, mode: 'copy' )
>
>  // indicates to use as container the value indicated in the parameter
>  container( params.CONTAINER )
>
>  // show in the log which input file is analysed
>  tag( "${genome}" )
>
>  input:
>  path( genome )
>
>  output:
>  tuple path( genome ), path( "*" ), emit: bwa_index
>
>  script:
>  """
>  bwa index ${genome}
>  """
> }
> ```

- Create a new file `bwa_align.nf` in the current modules directory; paste the following and save.

> ```bash
> /*
> ========================================================================================
>    BWA-ALIGN module
> ========================================================================================
>    Website: http://bio-bwa.sourceforge.net/
> ========================================================================================
> */
>
> // Parameter definitions
> params.CONTAINER = "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0"
> params.OUTPUT = "bwa_align"
>
> process BWA_ALIGN {
>
>  // where to store the results and in which way
>  publishDir( params.OUTPUT, mode: 'copy' )
>
>  // indicates to use as container the value indicated in the parameter
>  container( params.CONTAINER )
>
>  // show in the log which input file is analysed
>  tag( "${sample_id}" )
>
>  input:
>  tuple path( genome ), path( "*" ), val( sample_id ), path( reads )
>
>  output:
>  tuple val( sample_id ), path( "${sample_id}.aligned.bam" ), emit: aligned_bam
>
>  script:
>  """
>  INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
>  bwa mem \$INDEX ${reads} > ${sample_id}.aligned.sam
>  samtools view -S -b ${sample_id}.aligned.sam > ${sample_id}.aligned.bam
>  """
> }
> ```

- Create a new file `samtools.nf` in the current modules directory; paste the following and save.

> ```bash
> /*
> ========================================================================================
>    SAMTOOLS module
> ========================================================================================
>    Website: http://www.htslib.org/doc/samtools.html
> ========================================================================================
> */
>
> // Parameter definitions
> params.CONTAINER = "quay.io/biocontainers/samtools:1.14--hb421002_0"
> params.OUTPUT = "sorted_bam"
>
> process SAMTOOLS_SORT {
>
>  // where to store the results and in which way
>  publishDir( params.OUTPUT, mode: 'copy' )
>
>  // indicates to use as container the value indicated in the parameter
>  container( params.CONTAINER )
>
>  // show in the log which input file is analysed
>  tag( "${sample_id}" )
>
>  input:
>  tuple val( sample_id ), path( bam )
>
>  output:
>  tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam" ), emit: sorted_bam
>
>  script:
>  """
>  samtools sort -o "${sample_id}.aligned.sorted.bam" ${bam}
>  """
> }
>
> /*
> * Index the BAM file for visualization purpose
> */
> process SAMTOOLS_INDEX {
>
>  // where to store the results and in which way
>  publishDir( params.OUTPUT, mode: 'copy' )
>
>  // indicates to use as container the value indicated in the parameter
>  container( params.CONTAINER )
>
>  // show in the log which input file is analysed
>  tag( "${sample_id}" )
>
>  input:
>  tuple val( sample_id ), path( sortedbam )
>
>  output:
>  tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam" ), path("*"), emit: aligned_sorted_bam
>
>  script:
>  """
>  samtools index ${sortedbam}
>  """
> }
> ```

- Create a new file `bcftools.nf` in the current modules directory; paste the following and save.

> ```bash
> /*
> ========================================================================================
>    BCFTOOLS module
> ========================================================================================
>    Website: http://www.htslib.org/doc/bcftools.html
> ========================================================================================
> */
>
> // Parameter definitions
> params.CONTAINER = "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
> params.OUTPUT = "vcf"
>
> process BCFTOOLS_MPILEUP {
>
>  // where to store the results and in which way
>  publishDir( params.OUTPUT, mode: 'copy' )
>
>  // indicates to use as container the value indicated in the parameter
>  container( params.CONTAINER )
>
>  // show in the log which input file is analysed
>  tag( "${sample_id}" )
>
>  input:
>  tuple path( genome ), path( "*" ), val( sample_id ), path( sortedbam ), path("*")
>
>  output:
>  tuple val( sample_id ), path( "${sample_id}_raw.bcf" ), emit: raw_bcf
>
>  script:
>  """
>  bcftools mpileup -O b -o "${sample_id}_raw.bcf" -f ${genome} ${sortedbam}
>  """
> }
>
> /*
> * Detect the single nucleotide variants (SNVs).
> */
> process BCFTOOLS_CALL {
>
>  // where to store the results and in which way
>  publishDir( params.OUTPUT, mode: 'copy' )
>
>  // indicates to use as container the value indicated in the parameter
>  container( params.CONTAINER )
>
>  // show in the log which input file is analysed
>  tag( "${sample_id}" )
>
>  input:
>  tuple val( sample_id ), path( rawbcf )
>
>  output:
>  tuple val( sample_id ), path( "${sample_id}_variants.vcf" ), emit: variants_vcf
>
>  script:
>  """
>  bcftools call --ploidy 1 -m -v -o "${sample_id}_variants.vcf" ${rawbcf}
>  """
> }
>
> process VCFUTILS {
>  // where to store the results and in which way
>  publishDir( params.OUTPUT, mode: 'copy' )
>
>  // indicates to use as container the value indicated in the parameter
>  container( params.CONTAINER )
>
>  // show in the log which input file is analysed
>  tag( "${sample_id}" )
>
>  input:
>  tuple val( sample_id ), path( rawvcf )
>
>  output:
>  tuple val( sample_id ), path( "${sample_id}_final_variants.vcf" ), emit: final_variants_vcf
>
>  script:
>  """
>  vcfutils.pl varFilter ${rawvcf} > "${sample_id}_final_variants.vcf"
>  """
> }
> ```

## Importing Module Components in Nextflow

A component defined in a Nextflow module script can be readily imported into another Nextflow script using the `include` keyword. This provides flexibility and promotes the reusability of code.

Take for example, a process named `BWA_INDEX` that's defined in a module script `bwa_index.nf`. This process can be imported into the main execution context, and can therefore be invoked in the workflow scope.

By default, Nextflow looks for the script file `nextflow_tutorial/modules/bwa_index.nf`, resolving the path relative to the location of the including script.

> **Note**: Remember that relative paths must begin with the `./` prefix.

> **Remote Import Restriction**
>
> > Currently, it's not possible to include a script from a remote URL in the `from` statement.

### Multiple Inclusions

A Nextflow script can include any number of modules. When you need to include multiple components from the same module script, you can specify the component names in the same inclusion using curly brackets `{}`. Component names should be separated by a semi-colon `;` as shown below:

```nextflow
include { COMPONENT1; COMPONENT2; COMPONENT3 } from 'path/to/module'
```

A module script can also define one or more parameters using the same syntax as a Nextflow workflow script.

To illustrate this, let's create a new Nextflow script called `variant-calling.nf` in the `/workspace/nextflow_tutorial/` directory. This script will include multiple components from different module scripts:

> ```groovy
> /*
> ========================================================================================
>    Variant-Calling Nextflow Workflow
> ========================================================================================
>    Github   :
>    Contact  :
> ----------------------------------------------------------------------------------------
> */
>
> nextflow.enable.dsl=2
>
> // Display pipeline details
> println """\
>         V A R I A N T-C A L L I N G - N F   P I P E L I N E
>         ===================================
>         genome       : ${params.genome}
>         reads        : ${params.reads}
>         outdir       : ${params.outdir}
>         """
>         .stripIndent()
>
> /*
> ========================================================================================
>    Include Modules
> ========================================================================================
> */
>
> // Include modules with corresponding output directories
> include { FASTQC } from "./modules/fastqc" addParams(OUTPUT: "${params.outdir}/fastqc")
> include { BWA_INDEX } from "./modules/bwa_index" addParams(OUTPUT: "${params.outdir}/bwa_index")
> include { BWA_ALIGN } from "./modules/bwa_align" addParams(OUTPUT: "${params.outdir}/bwa_align")
> include { SAMTOOLS_SORT; SAMTOOLS_INDEX } from "./modules/samtools" addParams(OUTPUT: "${params.outdir}/sorted_bam")
> include { BCFTOOLS_MPILEUP; BCFTOOLS_CALL; VCFUTILS } from "./modules/bcftools" addParams(OUTPUT: "${params.outdir}/vcf")
>
> /*
> ========================================================================================
>    Create Channels
> ========================================================================================
> */
>
> // Create channels for reference genome and reads
> ref_ch = Channel.fromPath(params.genome, checkIfExists: true)
> reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
>
> /*
> ========================================================================================
>    WORKFLOW - Variant Calling
> ========================================================================================
> */
>
> workflow {
>    // Execute modules in the workflow
>    FASTQC(reads_ch)
>    BWA_INDEX(ref_ch)
>    BWA_ALIGN(BWA_INDEX.out.bwa_index.combine(reads_ch))
>    SAMTOOLS_SORT(BWA_ALIGN.out.aligned_bam)
>    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sorted_bam)
>    BCFTOOLS_MPILEUP(BWA_INDEX.out.bwa_index.combine(SAMTOOLS_INDEX.out.aligned_sorted_bam))
>    BCFTOOLS_CALL(BCFTOOLS_MPILEUP.out.raw_bcf)
>    VCFUTILS(BCFTOOLS_CALL.out.variants_vcf)
> }
>
> // Display pipeline execution summary upon completion
> workflow.onComplete {
>    println (workflow.success ? """
>        Pipeline execution summary
>        ---------------------------
>        Completed at: ${workflow.complete}
>        Duration    : ${workflow.duration}
>        Success     : ${workflow.success}
>        workDir     : ${workflow.workDir}
>        exit status : ${workflow.exitStatus}
>        """ : """
>        Failed: ${workflow.errorReport}
>        exit status : ${workflow.exitStatus}
>        """
>    )
> }
>
> /*
> ========================================================================================
>    THE END
> ========================================================================================
> */
> ```

Update the `nextflow.config` in the `/workspace/nextflow_tutorial/` directory to add the following specifications and save.

> ```bash
> /*
> ========================================================================================
>    NF-CORE Custom Config File
> ========================================================================================
>    Default config options for HPC compute environments
> ----------------------------------------------------------------------------------------
> */
>
> //Profile config names for nf-core/configs
>
> params {
>
>  config_profile_description = 'gitpod compliant config'
>  config_profile_contact     = ''
>  config_profile_url         = ''
>
>  // Input parameters
>
>  genome                     = "${launchDir}/data/ref_genome/ecoli_rel606.fasta"
>  reads                      = "${launchDir}/data/trimmed_fastq/*_{1,2}.trim.fastq.gz"
>
>  // Output options
>  outdir                     = "results"
> }
>
> /*
> ========================================================================================
>    Nextflow Metrics & Reports
> ========================================================================================
> */
>
> timeline {
>  enabled = true
>  file    = "${params.outdir}/timeline.html"
> }
>
> report {
>  enabled = true
>  file    = "${params.outdir}/report.html"
> }
> trace {
>  enabled = true
>  fields  = 'task_id,name,status,exit,realtime,%cpu,%mem,rss,vmem,peak_rss,peak_vmem,rchar,wchar'
>  file    = "${params.outdir}/trace.txt"
> }
>
> /*
> ========================================================================================
>    Base Executor config
> ========================================================================================
> */
>
> executor {
>  queueSize = 2
> }
>
> /*
> ========================================================================================
>    Profiles - slurm,singularity,conda
> ========================================================================================
> */
>
> profiles {
>  slurm {
>    process {
>      executor     = 'slurm'
>      queue        = 'my.q'
>    }
>    executor {
>      queueSize    = 100
>      pollInterval = '15 sec'
>    }
>  }
>
>  conda {
>    process.conda = "${launchDir}/environment.yml"
>  }
>
>  docker {
>    docker.enabled = true
>  }
>
>  singularity {
>    singularity.enabled = true
>  }
>
> }
> ```

To run the workflow using Docker profile on GitPod

```bash
nextflow run variant-calling.nf -profile docker
```

---

<details>
  <summary><b>CLICK HERE for SLURM submit script</b></summary>

<pre>
#!/bin/bash
#SBATCH --job-name=varcal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --partition=

cd nextflow_tutorial/
module load nextflow
nextflow run variant-calling.nf -profile singularity,slurm --outdir "results-slurm"
</pre>

## </details>

<br>

## Module Aliases in Nextflow

In Nextflow, a process component, like `BWA_INDEX`, can only be invoked once within the same workflow context. However, the inclusion of a module component allows you to specify a name `alias` using the `as` keyword in the `include` statement. This enables the inclusion and invocation of the same component multiple times in your script using different names.

For instance:

```groovy
nextflow.enable.dsl=2

// Including the same module twice with different aliases
include { BWA_INDEX } from './modules/bwa_index.nf'
include { BWA_INDEX as INDEX } from './modules/bwa_index.nf'

workflow {
  // Define channels
  ref_ch = Channel.fromPath( params.genome, checkIfExists: true  )
  reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )

  // Use the modules with different aliases
  BWA_INDEX( ref_ch )
  INDEX( ref_ch )
}
```

This allows the `BWA_INDEX` process to be invoked twice in the workflow under different aliases.

It's important to note that a module inherits the parameters defined before its `include` statement. Any parameters set after this statement are ignored. Therefore, it's advisable to define all pipeline parameters at the start of the script, prior to any `include` declarations.

The `addParams` option can be used to extend the module parameters without affecting the external scope. This ensures a clean and independent configuration for each module.

### Quick Recap

- A module file in Nextflow is a script containing one or more process definitions that can be imported into another Nextflow script.
- The `include` keyword is used to import a module into a workflow.
- Module aliases enable the reuse of the same module multiple times within a workflow.
- A module inherits parameters defined before its inclusion, and the `addParams` option allows for the extension of module parameters without affecting the external scope.

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_configuration" style="float: left"><b>Back to:</b>NextFlow Configuration</a>

<a href="/nextflow_varcal/nextflow/nextflow_sub_workflows" style="float: right"><b>Next:</b>NextFlow SubWorkflows</a></h5>
