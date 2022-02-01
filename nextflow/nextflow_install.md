---
layout: main
title: NextFlow Install
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_install
---
{% include _nextflow_nextflow_install_toc.html %}


<hr>
<center>This chapter is a part of <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Nextflow Install

### Requirements

* Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and [Java 8 (or later, up to 15)](http://www.oracle.com/technetwork/java/javase/downloads/index.html) to be installed.

* For the execution in a cluster of computers the use a shared file system is required to allow the sharing of tasks input/output files.

* Windows system is supported through WSL.

* Nextflow is distributed as a self-installing package, which means that it does not require any special installation procedure. Installation instructions can be found [here](https://www.nextflow.io/docs/latest/getstarted.html#installation)


```bash
mkdir -p ~/bin/
cd ~/bin/
curl -s https://get.nextflow.io | bash
chmod +x nextflow
echo 'PATH=~/bin:$PATH' >> ~/.bashrc
which nextflow      ## should point back to ~/bin
nextflow -v         ## check install by invoking Nextflow, getting version
```

---
## NF-CORE Install

```bash
conda install nf-core
```
or

```bash
pip install nf-core
```

<br>

## IDE

IDE choice will not be discussed in this tutorial. You can use what you are comfortable with, you will need to edit text files as well as run commands from the terminal, and will need to run those commands from a cluster node at one point. Some good options are VScode, with nextflow plugins and remote-ssh module setup and working. This is beyond the scope of the tutorial at this time.

*   **[Download Visual Studio Code](https://code.visualstudio.com/Download)**

*   [VSCode Nextflow Plugin](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow)

*   [VSCode Remote-SSH Plugin](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh)

## Environment Setup ()

*  Create `environment.yml`

```bash
channels:
  - bioconda
  - conda-forge
dependencies:
  - fastqc
  - trimmomatic
  - bwa
  - samtools
  - bcftools
  - multiqc
  - graphviz
```

```bash
mamba env create -n varcal -f environment.yml
```

## Data Download & Setup

*   **[Data Source - DataCarpentry Wrangling Genomics Lesson](https://datacarpentry.org/wrangling-genomics/02-quality-control/index.html)**

*  Download reference genome

```bash
cd ~/
mkdir -p nextflow_tutorial/data/ref_genome
cd nextflow_tutorial/data/ref_genome
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip data/ref_genome/ecoli_rel606.fasta.gz
ls
```


*  Download raw data

```bash
cd ~/
mkdir -p nextflow_tutorial/data/untrimmed_fastq/
cd nextflow_tutorial/data/untrimmed_fastq/

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz
```

```bash
conda activate varcal
```

*  Trim the reads

```bash
cd nextflow_tutorial/data/untrimmed_fastq/

for infile in *_1.fastq.gz
do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE ${infile} ${base}_2.fastq.gz \
                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
done
```

```bash
mkdir ../trimmed_fastq
mv *.trim* ../trimmed_fastq
cd ../trimmed_fastq
ls
```




---

<h5><a href="/nextflow_varcal/nextflow/index" style="float: left"><b>Back to:</b>Table of Contents</a>

<a href="/nextflow_varcal/nextflow/nextflow_intro" style="float: right"><b>Next:</b>Nextflow Introduction</a></h5>
