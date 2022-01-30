---
layout: main
title: NextFlow Modules
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_modules
---
{% include _nextflow_nextflow_modules_toc.html %}


<hr>
<center>This is part 11 of 14 of a <a href="/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Modules

* In most programming languages there is the concept of creating code blocks/modules that can be reused.

* Nextflow (DSL2) allows the definition of module scripts that can be included and shared across workflow pipelines.

* A module file is nothing more than a Nextflow script containing one or more process definitions that can be imported from another Nextflow script.

* A module can contain the definition of a function, process and workflow definitions.

```bash
cd ~/nextflow_tutorial/
mkdir modules
cd modules
```


*   Create a new file `fastqc.nf` in the current `~/nextflow_tutorial/modules/` directory; paste the following and save.

>```bash
>/*
>========================================================================================
>    FASTQC module
>========================================================================================
>    Website: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
>========================================================================================
>*/
>
>// Parameter definitions
>params.CONTAINER = "quay.io/biocontainers/fastqc:0.11.9--0"
>params.OUTPUT = "trim_fastqc_output"
>
>process FASTQC {
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
>}
>```

*   Create a new file `bwa_index.nf` in the current `~/nextflow_tutorial/modules/` directory; paste the following and save.

>```bash
>/*
>========================================================================================
>    BWA-INDEX module
>========================================================================================
>    Website: http://www.htslib.org/doc/samtools.html
>========================================================================================
>*/
>
> // Parameter definitions
>params.CONTAINER = "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
>params.OUTPUT = "bwa_index"
>
>process BWA_INDEX {
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
>}
>```

*   Create a new file `bwa_align.nf` in the current `~/nextflow_tutorial/modules/` directory; paste the following and save.

>```bash
>/*
>========================================================================================
>    BWA-ALIGN module
>========================================================================================
>    Website: http://bio-bwa.sourceforge.net/
>========================================================================================
>*/
>
> // Parameter definitions
>params.CONTAINER = "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0"
>params.OUTPUT = "bwa_align"
>
>process BWA_ALIGN {
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
>}
>```

*   Create a new file `samtools.nf` in the current `~/nextflow_tutorial/modules/` directory; paste the following and save.

>```bash
>/*
>========================================================================================
>    SAMTOOLS module
>========================================================================================
>    Website: http://www.htslib.org/doc/samtools.html
>========================================================================================
>*/
>
>// Parameter definitions
>params.CONTAINER = "quay.io/biocontainers/samtools:1.14--hb421002_0"
>params.OUTPUT = "sorted_bam"
>
>process SAMTOOLS_SORT {
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
>}
>
>/*
> * Index the BAM file for visualization purpose
> */
>process SAMTOOLS_INDEX {
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
>}
>```

*   Create a new file `bcftools.nf` in the current `~/nextflow_tutorial/modules/` directory; paste the following and save.

>```bash
>/*
>========================================================================================
>    BCFTOOLS module
>========================================================================================
>    Website: http://www.htslib.org/doc/bcftools.html
>========================================================================================
>*/
>
>// Parameter definitions
>params.CONTAINER = "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
>params.OUTPUT = "vcf"
>
>process BCFTOOLS_MPILEUP {
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
>     bcftools mpileup -O b -o "${sample_id}_raw.bcf" -f ${genome} ${sortedbam}
>  """
>}
>
>/*
> * Detect the single nucleotide variants (SNVs).
> */
>process BCFTOOLS_CALL {
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
>}
>
>process VCFUTILS {
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
>}
>```

## Importing module components

* A component defined in a module script can be imported into another Nextflow script using the `include` keyword.

* The above snippets includes a process with name `BWA_INDEX` defined in the module script `bwa_index.nf` in the main execution context, as such it can be invoked in the workflow scope.

* Nextflow implicitly looks for the script file `~/nextflow_tutorial/modules/bwa_index.nf` resolving the path against the including script location.

> Note: Relative paths must begin with the `./` prefix.

> Remote
>> You can not include a script from a remote URL in the from statement.

---

### Multiple inclusions

* A Nextflow script allows the inclusion of any number of modules. When multiple components need to be included from the some module script, the component names can be specified in the same inclusion using the curly brackets `{}`. Component names are separated by a semi-colon `;` as shown below

* A module script can define one or more parameters using the same syntax of a Nextflow workflow script

* Create a new nextflow script called `variant-calling.nf` in the `~/nextflow_tutorial` as shown below with `include` statements.

>```groovy
>/*
>========================================================================================
>    Variant-Calling Nextflow Workflow
>========================================================================================
>    Github   : 
>    Contact  :     
>----------------------------------------------------------------------------------------
>*/
>
>nextflow.enable.dsl=2
>
>println """\
>         V A R I A N T-C A L L I N G - N F   P I P E L I N E
>         ===================================
>         genome       : ${params.genome}
>         reads        : ${params.reads}
>         outdir       : ${params.outdir}
>         """
>         .stripIndent()
>
>/*
>========================================================================================
>    Include Modules
>========================================================================================
>*/
>
>include { FASTQC }                                    from "./modules/fastqc" addParams(OUTPUT: "${params.output}/fastqc")
>include { BWA_INDEX  }                                from "./modules/bwa_index" addParams(OUTPUT: "${params.output}/bwa_index")
>include { BWA_ALIGN  }                                from "./modules/bwa_align" addParams(OUTPUT: "${params.output}/bwa_align")
>include { SAMTOOLS_SORT; SAMTOOLS_INDEX }             from "./modules/samtools" addParams(OUTPUT: "${params.output}/sorted_bam")
>include { BCFTOOLS_MPILEUP; BCFTOOLS_CALL; VCFUTILS } from "./modules/bcftools" addParams(OUTPUT: "${params.output}/vcf")
>
>/*
>========================================================================================
>    Create Channels
>========================================================================================
>*/
>
>ref_ch = Channel.fromPath( params.genome, checkIfExists: true  )  
>reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true ) 
>
>/*
>========================================================================================
>    WORKFLOW - Variant Calling
>========================================================================================
>*/
>
>workflow {
>
>    FASTQC( reads_ch )
>    BWA_INDEX( ref_ch )
>    BWA_ALIGN( BWA_INDEX.out.bwa_index.combine(reads_ch) )
>    SAMTOOLS_SORT( BWA_ALIGN.out.aligned_bam )
>    SAMTOOLS_INDEX( SAMTOOLS_SORT.out.sorted_bam )
>    BCFTOOLS_MPILEUP( BWA_INDEX.out.bwa_index.combine(SAMTOOLS_INDEX.out.aligned_sorted_bam) )
>    BCFTOOLS_CALL( BCFTOOLS_MPILEUP.out.raw_bcf )
>    VCFUTILS( BCFTOOLS_CALL.out.variants_vcf )
>
>}
>
>workflow.onComplete {
>
>    println ( workflow.success ? """
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
>}
>
>/*
>========================================================================================
>    THE END
>========================================================================================
>*/
>```

Let us submit this job the cluster. In order to do that we have to update our `nextflow.config` in `~/nextflow_tutorial` directory. Add the following cluster specifications to `nextflow.config` and save. 

>```bash
>/*
>========================================================================================
>    NF-CORE Custom Config File
>========================================================================================
>    Default config options for HPC compute environments
>----------------------------------------------------------------------------------------
>*/
>
>//Profile config names for nf-core/configs
>
>params {
>
>  config_profile_description = ''
>  config_profile_contact     = ''
>  config_profile_url         = ''
>
>  // Input parameters
>
>  genome                     = "$HOME/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta"
>  reads                      = "$HOME/nextflow_tutorial/data/trimmed_fastq/*_{1,2}.trim.fastq.gz"
>  
>  // Output options
>  outdir                     = "results"
>}
>
>/*
>========================================================================================
>    Nextflow Metrics & Reports
>========================================================================================
>*/
>
>timeline {
>  enabled = true
>  file    = "${params.outdir}/timeline.html"
>}
>
>report {
>  enabled = true
>  file    = "${params.outdir}/report.html"
>}
>trace {
>  enabled = true
>  fields  = 'task_id,name,status,exit,realtime,%cpu,%mem,rss,vmem,peak_rss,peak_vmem,rchar,wchar'
>  file    = "${params.outdir}/trace.txt"
>}
>
>/*
>========================================================================================
>    Base Executor config
>========================================================================================
>*/
>
>executor {
>  queueSize = 2
>}
>
>/*
>========================================================================================
>    Profiles - slurm,singularity,conda
>========================================================================================
>*/
>
>profiles {
>  sge {
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
>    process.conda = "$HOME/nextflow_tutorial/environment.yml"
>  }
>  
>  singularity {
>    singularity.enabled = true
>  }
>
>}
>```

Now, let us create a shell script `varcal.sh` to submit the `nextflow run` command as shown below:

> NOTE: change the email address to your own

>```bash
>#!/bin/bash -l
>
># Assign Job-Name instead of defauly which is the name of job-script
>#$ -N NF_Variant_Calling
># Submit job to the default queue "short"
>#$ -q my.q
>
>cd ~/nextflow_tutorial
>module load nextflow
>nextflow run variant-calling.nf -profile singularity,sge,short --output "results"
>```

Now submit the job to cluster for execution

```bash
sbatch varcal.sh
```

To check the status of job 

```bash
squeue
```

### Module aliases

* A process component, such as BWA_INDEX, can be invoked only once in the same workflow context.

* However, when including a module component it’s possible to specify a name `alias` using the keyword as in the include statement. This allows the inclusion and the invocation of the same component multiple times in your script using different names.

For example:

```groovy
nextflow.enable.dsl=2

include { BWA_INDEX } from './modules/bwa_index.nf'
include { BWA_INDEX as INDEX } from './modules/bwa_index.nf'
// OR can also be represented as
// include { BWA_INDEX; BWA_INDEX as INDEX } from './modules/bwa_index.nf'

workflow {
  // channels
  ref_ch = Channel.fromPath( params.genome, checkIfExists: true  )  
  reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true ) 
  // 
  BWA_INDEX( ref_ch )
  INDEX( ref_ch )
}
```

* The module inherits the parameters defined before the `include` statement, therefore any further parameter set later is ignored.

> Tip: Define all pipeline parameters at the beginning of the script before any include declaration.

* The option `addParams` can be used to extend the module parameters without affecting the external scope.

> Quick Recap
> A module file is a Nextflow script containing one or more process definitions that can be imported from another Nextflow script.
> To import a module into a workflow use the `include` keyword.

---

## Using nf-core modules

* **[Click Here for nf-core/modules github](https://github.com/nf-core/modules/tree/master/modules)**

```bash
cd ~/nextflow_tutorial
```

To initiaite a nextflow pipeline in nf-core style and pre-loaded templates and configuration, we can use `nf-core create` command:

```bash
nf-core create
```

You will prompted to enter `Workflow Name`, `Description`, `Author`. This will Initialise a pipeline git repository in your current directory. 

>Output
>```bash
>Workflow Name: variantcalling
>Description: simple variant calling pipeline
>Author: Sateesh Peri
>INFO     Creating new nf-core pipeline: 'nf-core/variantcalling'                                                                    create.py:67
>INFO     Initialising pipeline git repository                                                                                      create.py:168
>INFO     Done. Remember to add a remote and push to GitHub:                                                                        create.py:175
>          cd nextflow_tutorial/nf-core-variantcalling                                                                   
>          git remote add origin git@github.com:USERNAME/REPO_NAME.git                                                                           
>          git push --all origin                                                                                                                 
>INFO     This will also push your newly created dev branch and the TEMPLATE branch for syncing.                                    create.py:181
>INFO     !!!!!! IMPORTANT !!!!!!                                                                                                    create.py:58
>                                                                                                                                                
>         If you are interested in adding your pipeline to the nf-core community,                                                                
>         PLEASE COME AND TALK TO US IN THE NF-CORE SLACK BEFORE WRITING ANY CODE!                                                               
>                                                                                                                                                
>         Please read: https://nf-co.re/developers/adding_pipelines#join-the-community
>```


```bash
tree nf-core-variantcalling/
```

```bash
nf-core-variantcalling/
├── assets
│   ├── email_template.html
│   ├── email_template.txt
│   ├── multiqc_config.yaml
│   ├── nf-core-variantcalling_logo.png
│   ├── samplesheet.csv
│   ├── schema_input.json
│   └── sendmail_template.txt
├── bin
│   ├── check_samplesheet.py
│   └── scrape_software_versions.py
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
│   │   └── nf-core-variantcalling_logo.png
│   ├── output.md
│   ├── README.md
│   └── usage.md
├── lib
│   ├── nfcore_external_java_deps.jar
│   ├── NfcoreSchema.groovy
│   ├── NfcoreTemplate.groovy
│   ├── Utils.groovy
│   ├── WorkflowMain.groovy
│   └── WorkflowVariantcalling.groovy
├── LICENSE
├── main.nf
├── modules
│   ├── local
│   │   ├── functions.nf
│   │   ├── get_software_versions.nf
│   │   └── samplesheet_check.nf
│   └── nf-core
│       └── modules
│           ├── fastqc
│           │   ├── functions.nf
│           │   ├── main.nf
│           │   └── meta.yml
│           └── multiqc
│               ├── functions.nf
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
    └── variantcalling.nf

15 directories, 47 files
```

```bash
cd nf-core-variantcalling/
```


```bash
nf-core modules list remote
```

---
<details>
  <summary><b>CLICK HERE for currently available modules</b></summary>

<pre>

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.1


                                                                                 INFO     Modules available from nf-core/modules (master):               list.py:122

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Module Name                          ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ abacas                               │
│ adapterremoval                       │
│ agrvate                              │
│ allelecounter                        │
│ amps                                 │
│ arriba                               │
│ artic/guppyplex                      │
│ artic/minion                         │
│ assemblyscan                         │
│ ataqv/ataqv                          │
│ bakta                                │
│ bamaligncleaner                      │
│ bamtools/split                       │
│ bamutil/trimbam                      │
│ bandage/image                        │
│ bbmap/align                          │
│ bbmap/bbduk                          │
│ bbmap/bbsplit                        │
│ bbmap/index                          │
│ bcftools/concat                      │
│ bcftools/consensus                   │
│ bcftools/filter                      │
│ bcftools/index                       │
│ bcftools/isec                        │
│ bcftools/merge                       │
│ bcftools/mpileup                     │
│ bcftools/norm                        │
│ bcftools/query                       │
│ bcftools/reheader                    │
│ bcftools/stats                       │
│ bcftools/view                        │
│ bedtools/bamtobed                    │
│ bedtools/complement                  │
│ bedtools/genomecov                   │
│ bedtools/getfasta                    │
│ bedtools/intersect                   │
│ bedtools/makewindows                 │
│ bedtools/maskfasta                   │
│ bedtools/merge                       │
│ bedtools/slop                        │
│ bedtools/sort                        │
│ bedtools/subtract                    │
│ bismark/align                        │
│ bismark/deduplicate                  │
│ bismark/genomepreparation            │
│ bismark/methylationextractor         │
│ bismark/report                       │
│ bismark/summary                      │
│ blast/blastn                         │
│ blast/makeblastdb                    │
│ bowtie/align                         │
│ bowtie/build                         │
│ bowtie2/align                        │
│ bowtie2/build                        │
│ bwa/aln                              │
│ bwa/index                            │
│ bwa/mem                              │
│ bwa/sampe                            │
│ bwa/samse                            │
│ bwamem2/index                        │
│ bwamem2/mem                          │
│ bwameth/align                        │
│ bwameth/index                        │
│ cat/cat                              │
│ cat/fastq                            │
│ cellranger/count                     │
│ cellranger/mkfastq                   │
│ cellranger/mkgtf                     │
│ cellranger/mkref                     │
│ checkm/lineagewf                     │
│ chromap/chromap                      │
│ chromap/index                        │
│ clonalframeml                        │
│ cmseq/polymut                        │
│ cnvkit/batch                         │
│ cooler/cload                         │
│ cooler/digest                        │
│ cooler/dump                          │
│ cooler/merge                         │
│ cooler/zoomify                       │
│ csvtk/concat                         │
│ csvtk/split                          │
│ custom/dumpsoftwareversions          │
│ custom/getchromsizes                 │
│ cutadapt                             │
│ damageprofiler                       │
│ dastool/dastool                      │
│ dastool/scaffolds2bin                │
│ dedup                                │
│ deeptools/computematrix              │
│ deeptools/plotfingerprint            │
│ deeptools/plotheatmap                │
│ deeptools/plotprofile                │
│ delly/call                           │
│ diamond/blastp                       │
│ diamond/blastx                       │
│ diamond/makedb                       │
│ dragmap/align                        │
│ dragmap/hashtable                    │
│ dragonflye                           │
│ dshbio/exportsegments                │
│ dshbio/filterbed                     │
│ dshbio/filtergff3                    │
│ dshbio/splitbed                      │
│ dshbio/splitgff3                     │
│ ectyper                              │
│ emmtyper                             │
│ ensemblvep                           │
│ expansionhunter                      │
│ fargene                              │
│ fastani                              │
│ fastp                                │
│ fastqc                               │
│ fastqscan                            │
│ fasttree                             │
│ fgbio/callmolecularconsensusreads    │
│ fgbio/fastqtobam                     │
│ fgbio/groupreadsbyumi                │
│ fgbio/sortbam                        │
│ filtlong                             │
│ flash                                │
│ freebayes                            │
│ gatk4/applybqsr                      │
│ gatk4/baserecalibrator               │
│ gatk4/bedtointervallist              │
│ gatk4/calculatecontamination         │
│ gatk4/createsequencedictionary       │
│ gatk4/createsomaticpanelofnormals    │
│ gatk4/estimatelibrarycomplexity      │
│ gatk4/fastqtosam                     │
│ gatk4/filtermutectcalls              │
│ gatk4/gatherbqsrreports              │
│ gatk4/genomicsdbimport               │
│ gatk4/genotypegvcfs                  │
│ gatk4/getpileupsummaries             │
│ gatk4/haplotypecaller                │
│ gatk4/indexfeaturefile               │
│ gatk4/intervallisttools              │
│ gatk4/learnreadorientationmodel      │
│ gatk4/markduplicates                 │
│ gatk4/mergebamalignment              │
│ gatk4/mergevcfs                      │
│ gatk4/mutect2                        │
│ gatk4/revertsam                      │
│ gatk4/samtofastq                     │
│ gatk4/splitncigarreads               │
│ gatk4/variantfiltration              │
│ genmap/index                         │
│ genmap/mappability                   │
│ genrich                              │
│ gffread                              │
│ glnexus                              │
│ graphmap2/align                      │
│ graphmap2/index                      │
│ gstama/collapse                      │
│ gstama/merge                         │
│ gtdbtk/classifywf                    │
│ gubbins                              │
│ gunc/downloaddb                      │
│ gunc/run                             │
│ gunzip                               │
│ hicap                                │
│ hifiasm                              │
│ hisat2/align                         │
│ hisat2/build                         │
│ hisat2/extractsplicesites            │
│ hmmcopy/gccounter                    │
│ hmmcopy/readcounter                  │
│ hmmer/hmmalign                       │
│ homer/annotatepeaks                  │
│ homer/findpeaks                      │
│ homer/maketagdirectory               │
│ homer/makeucscfile                   │
│ idr                                  │
│ imputeme/vcftoprs                    │
│ iqtree                               │
│ ismapper                             │
│ isoseq3/cluster                      │
│ isoseq3/refine                       │
│ ivar/consensus                       │
│ ivar/trim                            │
│ ivar/variants                        │
│ jupyternotebook                      │
│ kallisto/index                       │
│ kallistobustools/count               │
│ kallistobustools/ref                 │
│ khmer/normalizebymedian              │
│ kleborate                            │
│ kraken2/kraken2                      │
│ krona/kronadb                        │
│ krona/ktimporttaxonomy               │
│ last/dotplot                         │
│ last/lastal                          │
│ last/lastdb                          │
│ last/mafconvert                      │
│ last/mafswap                         │
│ last/postmask                        │
│ last/split                           │
│ last/train                           │
│ leehom                               │
│ lima                                 │
│ lissero                              │
│ lofreq/call                          │
│ lofreq/callparallel                  │
│ lofreq/filter                        │
│ lofreq/indelqual                     │
│ macrel/contigs                       │
│ macs2/callpeak                       │
│ malt/build                           │
│ malt/run                             │
│ maltextract                          │
│ manta/germline                       │
│ manta/somatic                        │
│ manta/tumoronly                      │
│ mapdamage2                           │
│ mash/sketch                          │
│ mashtree                             │
│ maxbin2                              │
│ medaka                               │
│ megahit                              │
│ meningotype                          │
│ metabat2/jgisummarizebamcontigdepths │
│ metabat2/metabat2                    │
│ metaphlan3                           │
│ methyldackel/extract                 │
│ methyldackel/mbias                   │
│ minia                                │
│ miniasm                              │
│ minimap2/align                       │
│ minimap2/index                       │
│ mlst                                 │
│ mosdepth                             │
│ msisensor/msi                        │
│ msisensor/scan                       │
│ mtnucratio                           │
│ multiqc                              │
│ mummer                               │
│ muscle                               │
│ nanolyse                             │
│ nanoplot                             │
│ ncbigenomedownload                   │
│ nextclade                            │
│ ngmaster                             │
│ nucmer                               │
│ optitype                             │
│ pairix                               │
│ pairtools/dedup                      │
│ pairtools/flip                       │
│ pairtools/parse                      │
│ pairtools/restrict                   │
│ pairtools/select                     │
│ pairtools/sort                       │
│ pangolin                             │
│ paraclu                              │
│ pbbam/pbmerge                        │
│ pbccs                                │
│ peddy                                │
│ phantompeakqualtools                 │
│ phyloflash                           │
│ picard/collecthsmetrics              │
│ picard/collectmultiplemetrics        │
│ picard/collectwgsmetrics             │
│ picard/filtersamreads                │
│ picard/markduplicates                │
│ picard/mergesamfiles                 │
│ picard/sortsam                       │
│ pirate                               │
│ plasmidid                            │
│ plink/extract                        │
│ plink/vcf                            │
│ plink2/vcf                           │
│ pmdtools/filter                      │
│ porechop                             │
│ preseq/lcextrap                      │
│ prodigal                             │
│ prokka                               │
│ pycoqc                               │
│ pydamage/analyze                     │
│ pydamage/filter                      │
│ qcat                                 │
│ qualimap/bamqc                       │
│ qualimap/rnaseq                      │
│ quast                                │
│ racon                                │
│ rapidnj                              │
│ rasusa                               │
│ raxmlng                              │
│ rmarkdownnotebook                    │
│ roary                                │
│ rsem/calculateexpression             │
│ rsem/preparereference                │
│ rseqc/bamstat                        │
│ rseqc/inferexperiment                │
│ rseqc/innerdistance                  │
│ rseqc/junctionannotation             │
│ rseqc/junctionsaturation             │
│ rseqc/readdistribution               │
│ rseqc/readduplication                │
│ salmon/index                         │
│ salmon/quant                         │
│ samblaster                           │
│ samtools/ampliconclip                │
│ samtools/bam2fq                      │
│ samtools/depth                       │
│ samtools/faidx                       │
│ samtools/fastq                       │
│ samtools/fixmate                     │
│ samtools/flagstat                    │
│ samtools/idxstats                    │
│ samtools/index                       │
│ samtools/merge                       │
│ samtools/mpileup                     │
│ samtools/sort                        │
│ samtools/stats                       │
│ samtools/view                        │
│ scoary                               │
│ seacr/callpeak                       │
│ seqkit/split2                        │
│ seqsero2                             │
│ seqtk/mergepe                        │
│ seqtk/sample                         │
│ seqtk/subseq                         │
│ sequenzautils/bam2seqz               │
│ sequenzautils/gcwiggle               │
│ seqwish/induce                       │
│ shovill                              │
│ snpdists                             │
│ snpeff                               │
│ snpsites                             │
│ sortmerna                            │
│ spades                               │
│ spatyper                             │
│ sratools/fasterqdump                 │
│ sratools/prefetch                    │
│ staphopiasccmec                      │
│ star/align                           │
│ star/genomegenerate                  │
│ strelka/germline                     │
│ strelka/somatic                      │
│ stringtie/merge                      │
│ stringtie/stringtie                  │
│ subread/featurecounts                │
│ tabix/bgzip                          │
│ tabix/bgziptabix                     │
│ tabix/tabix                          │
│ tbprofiler/profile                   │
│ tiddit/cov                           │
│ tiddit/sv                            │
│ trimgalore                           │
│ ucsc/bed12tobigbed                   │
│ ucsc/bedclip                         │
│ ucsc/bedgraphtobigwig                │
│ ucsc/bigwigaverageoverbed            │
│ ucsc/liftover                        │
│ ucsc/wigtobigwig                     │
│ ultra/pipeline                       │
│ umitools/dedup                       │
│ umitools/extract                     │
│ unicycler                            │
│ untar                                │
│ unzip                                │
│ variantbam                           │
│ vcftools                             │
│ yara/index                           │
│ yara/mapper                          │
└──────────────────────────────────────┘
</pre>

</details>
---
<br>



```bash
nf-core modules install bwa/index
```

```bash

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.1



INFO     Installing 'bwa/index'                                                                                                   install.py:127
INFO     Downloaded 2 files to ./modules/nf-core/modules/bwa/index                                                        modules_command.py:268
```

```bash
nf-core modules install bwa/mem
```

```bash
INFO     Installing 'bwa/index'                                                                                                   install.py:127
INFO     Downloaded 2 files to ./modules/nf-core/modules/bwa/index                                                        modules_command.py:268
```

```bash
nf-core modules install samtools/sort
```

```bash
INFO     Installing 'samtools/sort'                                                                                               install.py:127
INFO     Downloaded 2 files to ./modules/nf-core/modules/samtools/sort                                                    modules_command.py:268
```

```bash
nf-core modules install samtools/index
```

```bash
INFO     Installing 'samtools/index'                                                                                              install.py:127
INFO     Downloaded 2 files to ./modules/nf-core/modules/samtools/index                                                   modules_command.py:268
```

```bash
nf-core modules install bcftools/mpileup
```

```bash
INFO     Installing 'bcftools/mpileup'                                                                                            install.py:127
INFO     Downloaded 2 files to ./modules/nf-core/modules/bcftools/mpileup                                                 modules_command.py:268
```

```bash
tree modules
```

```bash
modules
├── local
│   ├── functions.nf
│   ├── get_software_versions.nf
│   └── samplesheet_check.nf
└── nf-core
    └── modules
        ├── bwa
        │   ├── index
        │   │   ├── main.nf
        │   │   └── meta.yml
        │   └── mem
        │       ├── main.nf
        │       └── meta.yml
        ├── fastqc
        │   ├── functions.nf
        │   ├── main.nf
        │   └── meta.yml
        └── multiqc
            ├── functions.nf
            ├── main.nf
            └── meta.yml

<truncated>
```

* Every module folder will have
  1. `main.nf`  -
  2. `meta.yml` -











---

<h5><a href="/nextflow_varcal/nextflow/nextflow_configuration" style="float: left"><b>Back to:</b>NextFlow Configuration</a>

<a href="/nextflow_varcal/nextflow/nextflow_sub_workflows" style="float: right"><b>Next:</b>NextFlow SubWorkflows</a></h5>
