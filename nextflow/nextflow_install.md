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

## GITPOD 

Use a preconfigured Nextflow development environment using Gitpod.

>Requirements:
>* A GitHub, GitLab or BitBucket account (gives you 50 hours usage per month)
>* Web browser (Google Chrome, Firefox)

## Gitpod environment

To run Gitpod:

1) click the following URL:

gitpod.io/#https://github.com/sateeshperi/nextflow_tutorial.git

(which is our Github repository URL, prefixed with `gitpod.io/#`).

OR

2) click the button below to launch directly as well

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/sateeshperi/nextflow_tutorial.git)

OR

3) Install the Gitpod browser extension, following the instructions **[here](https://www.gitpod.io/docs/browser-extension)**. This will add a green Gitpod button to each Git repository for easy access to the environment. Then go to the git repository for the training (https://github.com/sateeshperi/nextflow_tutorial.git), and click the green button.

<img src="images/gitpodbutton.png" alt="drawing" width="300"/>


* Login to your Github account (and allow authorization).
* You may need to close the window at each stage and click the gitpod link again.
* For BitBucket or GitLab, you may have to sign into Gitpod and change your integration on the following page: [gitpod.io/integrations](https://gitpod.io/integrations).

Once you have signed in, gitpod will load (skip prebuild, if asked):

<img src="images/gitpod.png" alt="drawing" width="300"/>

* Gitpod gives you 50 hours per month to run the environment, with up to 4 parallel workspaces at a time. So feel free to come back at any time and try the course at your own pace.

* It is useful to save the URL that is now your Gitpod environment. Later, we can use this to return to the same closed environment (so please take note of it now).

> Explore your Gitpod IDE

You should see something similar to the following:

<img src="images/gitpod_startup.png" alt="drawing" width="300"/>

* The sidebar allows you to customize your Gitpod environment and perform basic tasks (Copy/Paste, Open files, search, git, etc.) Click the **Explorer** button to see which files are in this repository.

* The terminal allows you to run all the programs in the repository, for example **nextflow, conda and docker are installed**.

* The main window allows you to view and edit files. Clicking on a file in the explorer will open it within the main window.

> To save your files, choose your file of interest, then either use the left/top hand side bar (File/Save As…​) or use you’re keyboard shortcut to save as (On Mac: `CMD, shift, s`), then choose Show local. This will open up a explorer to choose where to save your file on your local machine.

* **The GitPod comes with Nextflow, Conda and Docker pre-installed**

### 1. Create `varcal` environment based on yml file

```bash
mamba env update -n base -f environment.yml
```

### 2. Download reference genome and raw reads

```bash
bash data/fetch_raw_data.sh
```

### 3. Trim raw reads using `trimmomatic`

```bash
bash data/trim.sh
```

## Reopening a Gitpod session

* Any running workspace will be automatically stopped after 30 minutes. You can open the environment again by going to `gitpod.io/workspaces` and finding your previous environment, then clicking the three dot button to the right, and selecting Open.

* If you save the URL from your previous Gitpod environment, you can just paste this into your browser to open the previous environment. Environments are saved for up to two weeks, but don’t rely on their existence, download any important files you want for posterity.

* Alternatively, you can start a new workspace by clicking the green gitpod button, or following the Gitpod URL: gitpod.io/#https://gitpod.io/#https://github.com/sateeshperi/nextflow_tutorial.git

---
<br>
---
<details>
  <summary><b>CLICK HERE for Manual Install Instructions</b></summary>

### Nextflow Requirements

* Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and [Java 8 (or later, up to 15)](http://www.oracle.com/technetwork/java/javase/downloads/index.html) to be installed.
* For the execution in a cluster of computers the use a shared file system is required to allow the sharing of tasks input/output files.
* Windows system is supported through WSL.
* Nextflow is distributed as a self-installing package, which means that it does not require any special installation procedure. Installation instructions can be found [here](https://www.nextflow.io/docs/latest/getstarted.html#installation)  
  
<pre>
mkdir -p ~/bin/
cd ~/bin/
curl -s https://get.nextflow.io | bash
chmod +x nextflow
echo 'PATH=~/bin:$PATH' >> ~/.bashrc
which nextflow      ## should point back to ~/bin
nextflow -v         ## check install by invoking Nextflow, getting version
</pre>
  
## NF-CORE Install  
<pre>
conda install nf-core
# or
pip install nf-core
</pre>
  
## Environment Setup

*  Create `environment.yml`
<pre>
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
</pre>
  
<pre>
mamba env create -n varcal -f environment.yml
</pre>
  
  
## Data Download & Setup

*   **[Data Source - DataCarpentry Wrangling Genomics Lesson](https://datacarpentry.org/wrangling-genomics/02-quality-control/index.html)**

*  Download reference genome  
<pre>
cd ~/
mkdir -p nextflow_tutorial/data/ref_genome
cd nextflow_tutorial/data/ref_genome
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip data/ref_genome/ecoli_rel606.fasta.gz
ls
</pre>
  
*  Download raw data
<pre>
cd ~/
mkdir -p nextflow_tutorial/data/untrimmed_fastq/
cd nextflow_tutorial/data/untrimmed_fastq/

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz
</pre>  
  
*  Trim the reads
<pre>
conda activate varcal
</pre>
  
  
<pre>
cd nextflow_tutorial/data/untrimmed_fastq/

for infile in *_1.fastq.gz
do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE ${infile} ${base}_2.fastq.gz \
                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
done
</pre>  
  
<pre>
mkdir ../trimmed_fastq
mv *.trim* ../trimmed_fastq
cd ../trimmed_fastq
ls
</pre>
    
## IDE

IDE choice will not be discussed in this tutorial. You can use what you are comfortable with, you will need to edit text files as well as run commands from the terminal, and will need to run those commands from a cluster node at one point. Some good options are VScode, with nextflow plugins and remote-ssh module setup and working. This is beyond the scope of the tutorial at this time.

*   **[Download Visual Studio Code](https://code.visualstudio.com/Download)**

*   [VSCode Nextflow Plugin](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow)

*   [VSCode Remote-SSH Plugin](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh)  

</details>
---
<br>

---

<h5><a href="/nextflow_varcal/nextflow/index" style="float: left"><b>Back to:</b>Table of Contents</a>

<a href="/nextflow_varcal/nextflow/nextflow_intro" style="float: right"><b>Next:</b>Nextflow Introduction</a></h5>
