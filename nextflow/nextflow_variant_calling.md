---
layout: main
title: Variant Calling Workflow
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_variant_calling
---
{% include _nextflow_nextflow_variant_calling_toc.html %}


<hr>
<center>This is part 9 of 14 of a <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

```bash
mkdir workflow
cd workflow
```

```bash
conda activate varcal
```

## Variant Calling Workflow

Our variant calling workflow has the following steps:

1.  Index the reference genome for use by `bwa`.
2.  Align reads to reference genome using `bwa mem`.
3.  Convert the aligned SAM file to BAM format using `samtools`.
3.  Convert the format of the alignment to sorted BAM, with some intermediate steps.
4.  Calculate the read coverage of positions in the genome.
5.  Detect the single nucleotide variants (SNVs).
6.  Filter and report the SNVs in VCF (variant calling format).

![](images/variant_calling_workflow.png)

## Variant-Calling BASH script

```bash
set -e
mkdir /workspace/nextflow_tutorial/bash-results
cd /workspace/nextflow_tutorial/bash-results

genome=/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta

bwa index $genome

mkdir -p sam bam bcf vcf

for fq1 in /workspace/nextflow_tutorial/data/trimmed_fastq/*_1.trim.fastq.gz; do echo "working with file $fq1"; base=$(basename $fq1 _1.trim.fastq.gz);  echo "base name is $base" \

    fq1=/workspace/nextflow_tutorial/data/trimmed_fastq/${base}_1.trim.fastq.gz
   fq2=/workspace/nextflow_tutorial/data/trimmed_fastq/${base}_2.trim.fastq.gz
    sam=/workspace/nextflow_tutorial/results/sam/${base}.aligned.sam
    bam=/workspace/nextflow_tutorial/results/bam/${base}.aligned.bam
    sorted_bam=/workspace/nextflow_tutorial/results/bam/${base}.aligned.sorted.bam
    raw_bcf=/workspace/nextflow_tutorial/results/bcf/${base}_raw.bcf
    variants=/workspace/nextflow_tutorial/results/vcf/${base}_variants.vcf
    final_variants=/workspace/nextflow_tutorial/results/vcf/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done
```

## Variant Calling Nextflow Pipeline

```bash
conda activate varcal
```

---

## Finish the Variant-Calling Nextflow Workflow

```bash
code variant-calling.nf
```

```groovy
/*
========================================================================================
   Variant-Calling Nextflow Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Pipeline Input parameters

params.outdir = 'results'
params.genome = "/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta"
params.reads = "/workspace/nextflow_tutorial/data/trimmed_fastq/*_{1,2}.trim.fastq.gz"

println """\
         V A R I A N T-C A L L I N G - N F   P I P E L I N E
         ===================================
         genome       : ${params.genome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/*
========================================================================================
   Create Channels
========================================================================================
*/

ref_ch = Channel.fromPath( params.genome, checkIfExists: true )  
reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true ) 

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

    FASTQC( reads_ch )
    BWA_INDEX( ref_ch )
    BWA_ALIGN( BWA_INDEX.out.bwa_index.combine(reads_ch) ) // https://www.nextflow.io/docs/latest/process.html#understand-how-multiple-input-channels-work
    SAMTOOLS_SORT( BWA_ALIGN.out.aligned_bam )
    // Enter the rest of the processes for variant calling based on the bash script below

}

/*
========================================================================================
   Processes
========================================================================================
*/

/*
 * Align reads to reference genome & create BAM file.
 */
process FASTQC {
    tag{"FASTQC ${reads}"}
    label 'process_low'

    publishDir("${params.outdir}/fastqc_trim", mode: 'copy')
    
    input:
    tuple val( sample_id ), path( reads )

    output:
    path( "*_fastqc*" )

    script:
    """
    fastqc ${reads}
    """
}

/*
 * Index the reference genome for use by bwa and samtools.
 */
process BWA_INDEX {
  tag{"BWA_INDEX ${genome}"}
  label 'process_low'

  publishDir("${params.outdir}/bwa_index", mode: 'copy')
  
  input:
  path genome

  output:
  tuple path( genome ), path( "*" ), emit: bwa_index

  script:
  """
  bwa index ${genome} 
  """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BWA_ALIGN {
    tag{"BWA_ALIGN ${sample_id}"}
    label 'process_medium'

    publishDir("${params.outdir}/bwa_align", mode: 'copy')
    
    input:
    tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.bam" ), emit: aligned_bam

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem \$INDEX ${reads} > ${sample_id}.aligned.sam
    samtools view -S -b ${sample_id}.aligned.sam > ${sample_id}.aligned.bam
    """
}

/*
 * Convert the format of the alignment to sorted BAM.
 */
process SAMTOOLS_SORT {
  tag{"SAMTOOLS_SORT ${sample_id}"}
  label 'process_low'

  publishDir("${params.outdir}/bam_align", mode: 'copy')

  input:
  tuple val( sample_id ), path( bam )

  output:
  tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam" ), emit: sorted_bam

  script:
  """
  samtools sort -o "${sample_id}.aligned.sorted.bam" ${bam}
  """
}

/*
 * Index the BAM file for visualization purpose
 */
process SAMTOOLS_INDEX {

}

/*
 * Calculate the read coverage of positions in the genome.
 */
process BCFTOOLS_MPILEUP {

}

/*
 * Detect the single nucleotide variants (SNVs).
 */
process BCFTOOLS_CALL {

}

/*
 * Filter and report the SNVs in VCF (variant calling format).
 */
process VCFUTILS {

}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/
```

```bash
nextflow run variant-calling.nf 
```

> **[Understand how multiple input channels work](https://www.nextflow.io/docs/latest/process.html#understand-how-multiple-input-channels-work)**

## Handle completion event

*   This step shows how to execute an action when the pipeline completes the execution.

Note: that Nextflow processes define the execution of asynchronous tasks i.e. they are not executed one after another as they are written in the pipeline script as it would happen in a common imperative programming language.

The above script uses the `workflow.onComplete` event handler to print a confirmation message when the script completes.

```groovy
workflow.onComplete {

   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}
```
This code uses the ternary operator that is a shortcut expression that is equivalent to an if/else branch assigning some value to a variable.

>```groovy
>If expression is true ? "set value to a" : "else set value to b"
>```


## Metrics and reports

Nextflow is able to produce multiple reports and charts providing several runtime metrics and execution information.

*   The `-with-report` option enables the creation of the workflow execution report.

*   The `-with-trace` option enables the create of a tab separated file containing runtime information for each executed task, including: submission time, start time, completion time, cpu and memory used..

*   The `-with-timeline` option enables the creation of the workflow timeline report showing how processes where executed along time. This may be useful to identify most time consuming tasks and bottlenecks. See an example at this link.

*   The `-with-dag` option enables to rendering of the workflow execution direct acyclic graph representation. Note: this feature requires the installation of Graphviz, an open source graph visualization software, in your system.

More information can be found [here](https://www.nextflow.io/docs/latest/tracing.html).

```bash
nextflow run variant-calling.nf -resume -with-report -with-trace -with-timeline -with-dag dag.png
```

your final nextflow workflow DAG image should look something like this:

![](images/variant_calling_dag.png)


---

> Quick Recap
>*  Nextflow can combined tasks (processes) and manage data flows using channels into a single pipeline/workflow.
>*  A Workflow can be parameterise using params . These value of the parameters can be captured in a log file using `log.info`.
>*  Workflow steps are connected via their inputs and outputs using Channels.
>*   Nextflow can execute an action when the pipeline completes the execution using the `workflow.onComplete` event handler to print a confirmation message.
>*   Nextflow is able to produce multiple reports and charts providing several runtime metrics and execution information using the command line options `-with-report`, `-with-trace`, `-with-timeline` and produce a graph using `-with-dag`.

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_operators" style="float: left"><b>Back to:</b>Nextflow Operators</a>

<a href="/nextflow_varcal/nextflow/nextflow_configuration" style="float: right"><b>Next:</b>NextFlow Configuration</a></h5>
