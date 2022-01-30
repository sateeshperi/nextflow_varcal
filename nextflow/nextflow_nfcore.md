---
layout: main
title: NF-Core
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_nfcore
---
{% include _nextflow_nextflow_nfcore_toc.html %}


<hr>
<center>This is part 2 of 14 of a <a href="/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Learning Objectives
* Understand what nf-core is and how it relates to Nextflow, and how to use the nf-core helper tool to find nf-core pipelines.
* Run a nf-core pipeline (`nf-core/bactmap`) using a test dataset.
* Run a nf-core pipeline (`nf-core/bactmap`) using a test dataset on HPC

<br>

## What is nf-core?

**[nf-core](https://nf-co.re/)** is a community-led project to develop a set of best-practice pipelines built using Nextflow. Pipelines are governed by a **set of guidelines**, enforced by **community code reviews** and **automatic code testing**.

<img src="images/what_is_nf-core.png" alt="drawing" width="600"/>

<br>

## What is nf-core tools?

* nf-core provides a suite of helper tools aim to help users run and develop pipelines. These allow you to list all available pipelines and versions, with information about what versions you're running locally. There are also commands to help downloading pipelines for use offline.


### nf-core tools sub-commands

You can use the `--help` option to see the range of nf-core tools sub-commands.

```bash
nf-core --help
```

>```bash
>
>                                          ,--./,-.
>          ___     __   __   __   ___     /,-._.--~\
>    |\ | |__  __ /  ` /  \ |__) |__         }  {
>    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
>                                          `._,._,'
>
>    nf-core/tools version 2.1
>
>
>
>Usage: nf-core [OPTIONS] COMMAND [ARGS]...
>
>Options:
>  --version                  Show the version and exit.
>  -v, --verbose              Print verbose output to the console.
>  -l, --log-file <filename>  Save a verbose log to a file.
>  --help                     Show this message and exit.
>
>Commands:
>  list          List available nf-core pipelines with local info.
>  launch        Launch a pipeline using a web GUI or command line prompts.
>  download      Download a pipeline, nf-core/configs and pipeline...
>  licences      List software licences for a given workflow.
>  create        Create a new pipeline using the nf-core template.
>  lint          Check pipeline code against nf-core guidelines.
>  modules       Tools to manage Nextflow DSL2 modules as hosted on...
>  schema        Suite of tools for developers to manage pipeline schema.
>  bump-version  Update nf-core pipeline version number.
>  sync          Sync a pipeline TEMPLATE branch with the nf-core template.
>```

## Listing available nf-core pipelines

The simplest sub-command is `nf-core list`, which lists all available nf-core pipelines in the nf-core Github repository.

The output shows the latest version number and when that was released. If the pipeline has been pulled locally using Nextflow, it tells you when that was and whether you have the latest version.

Run the command below.

```bash
nf-core list
```

>```bash
>                                          ,--./,-.
>          ___     __   __   __   ___     /,-._.--~\
>    |\ | |__  __ /  ` /  \ |__) |__         }  {
>    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
>                                          `._,._,'
>
>    nf-core/tools version 2.1
>
>
>
>┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
>┃ Pipeline Name     ┃ Stars ┃ Latest Release ┃      Released ┃ Last Pulled ┃ Have latest release? ┃
>┡━━━━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
>│ rnaseq            │   384 │            3.4 │   1 weeks ago │           - │ -                    │
>│ fetchngs          │    32 │            1.3 │   4 weeks ago │           - │ -                    │
>│ eager             │    56 │          2.4.0 │   4 weeks ago │           - │ -                    │
>│ ampliseq          │    70 │          2.1.0 │   4 weeks ago │           - │ -                    │
>│ mhcquant          │    15 │          2.0.0 │  1 months ago │           - │ -                    │
>│ bacass            │    30 │          2.0.0 │  2 months ago │           - │ -                    │
>│ viralrecon        │    49 │            2.2 │  3 months ago │           - │ -                    │
>│ mag               │    63 │          2.1.0 │  3 months ago │           - │ -                    │
>│ bcellmagic        │    15 │          2.0.0 │  3 months ago │           - │ -                    │
>│ bactmap           │    24 │          1.0.0 │  4 months ago │           - │ -                    │
>│ smrnaseq          │    28 │          1.1.0 │  4 months ago │           - │ -                    │
>│ sarek             │   118 │          2.7.1 │  4 months ago │           - │ -                    │
>│ hic               │    28 │          1.3.0 │  5 months ago │           - │ -                    │
>│ metaboigniter     │     7 │          1.0.1 │  5 months ago │           - │ -                    │
>│ methylseq         │    72 │          1.6.1 │  5 months ago │           - │ -                    │
>│ diaproteomics     │     7 │          1.2.4 │  6 months ago │           - │ -                    │
>│ clipseq           │     8 │          1.0.0 │  6 months ago │           - │ -                    │
>│ pgdb              │     3 │          1.0.0 │  6 months ago │           - │ -                    │
>│ chipseq           │    88 │          1.2.2 │  6 months ago │           - │ -                    │
>│ scrnaseq          │    35 │          1.1.0 │  7 months ago │           - │ -                    │
>│ dualrnaseq        │     4 │          1.0.0 │  8 months ago │           - │ -                    │
>│ cageseq           │     7 │          1.0.2 │  9 months ago │           - │ -                    │
>│ nanoseq           │    46 │          1.1.0 │ 11 months ago │           - │ -                    │
>│ epitopeprediction │    12 │          1.1.0 │ 12 months ago │           - │ -                    │
>│ proteomicslfq     │    22 │          1.0.0 │ 12 months ago │           - │ -                    │
>│ hlatyping         │    25 │          1.2.0 │   1 years ago │           - │ -                    │
>│ atacseq           │    90 │          1.2.1 │   1 years ago │           - │ -                    │
>│ rnafusion         │    61 │          1.2.0 │   1 years ago │           - │ -                    │
>│ imcyto            │    11 │          1.0.0 │   1 years ago │           - │ -                    │
>│ slamseq           │     4 │          1.0.0 │   1 years ago │           - │ -                    │
>│ coproid           │     5 │            1.1 │   1 years ago │           - │ -                    │
>│ nascent           │     3 │            1.0 │   2 years ago │           - │ -                    │
>│ circrna           │     9 │            dev │             - │           - │ -                    │
>│ crisprvar         │     0 │            dev │             - │           - │ -                    │
>│ cutandrun         │    12 │            dev │             - │           - │ -                    │
>│ ddamsproteomics   │     4 │            dev │             - │           - │ -                    │
>│ demultiplex       │     9 │            dev │             - │           - │ -                    │
>│ denovohybrid      │     3 │            dev │             - │           - │ -                    │
>│ gwas              │     5 │            dev │             - │           - │ -                    │
>│ hicar             │     1 │            dev │             - │           - │ -                    │
>│ kmermaid          │    13 │            dev │             - │           - │ -                    │
>│ liverctanalysis   │     0 │            dev │             - │           - │ -                    │
>│ lncpipe           │    21 │            dev │             - │           - │ -                    │
>│ mnaseseq          │     7 │            dev │             - │           - │ -                    │
>│ pangenome         │    12 │            dev │             - │           - │ -                    │
>│ raredisease       │     8 │            dev │             - │           - │ -                    │
>│ rnavar            │     0 │            dev │             - │           - │ -                    │
>│ scflow            │     9 │            dev │             - │           - │ -                    │
>└───────────────────┴───────┴────────────────┴───────────────┴─────────────┴──────────────────────┘
>```

---

## Running nf-core pipelines

## Usage instructions and documentation

*   You can find general documentation and instructions for Nextflow and nf-core on the [nf-core website](https://nf-co.re/). Pipeline-specific documentation is bundled with each pipeline in the `/docs` folder. This can be read either locally, on GitHub, or on the nf-core website.

*   Each pipeline has its own webpage on nf-core website e.g. [nf-co.re/bactmap](https://nf-co.re/bactmap)

*   In addition to this documentation, each pipeline comes with basic command line reference. This can be seen by running the pipeline with the `--help` flag, for example:

```bash
nextflow run nf-core/bactmap -r 1.0.0 --help
```

>```bash
>------------------------------------------------------
>                                        ,--./,-.
>        ___     __   __   __   ___     /,-._.--~'
>  |\ | |__  __ /  ` /  \ |__) |__         }  {
>  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
>                                        `._,._,'
>  nf-core/bactmap v1.0.0
>------------------------------------------------------
>Typical pipeline command:
>
>  nextflow run nf-core/bactmap --input samplesheet.csv --reference ref.fasta -profile docker
>
>Input/output options
>  --input                             [string]  Path to a sample sheet describing paths to input fastq files
>  --outdir                            [string]  The output directory where the results will be saved. [default: ./results]
>  --email                             [string]  Email address for completion summary.
>
>Compulsory parameters
>  --reference                         [string]  Path to a fasta file of the reference sequence
>
>Optional pipeline steps
>  --trim                              [boolean] Trim reads [default: true]
>  --save_trimmed_fail                 [boolean] Saved failed read files after trimminng
>  --adapter_file                      [string]  path to file containing adapters in fasta format [default: ${baseDir}/assets/adapters.fas]
>  --subsampling_off                   [boolean] Turn off subsampling
>  --subsampling_depth_cutoff          [integer] Desired coverage depth when subsampling [default: 100]
>  --genome_size                       [string]  Specify genome size for subsampling rather than estimation using mash sketch
>  --remove_recombination              [boolean] Remove recombination using gubbins
>  --non_GATC_threshold                [number]  Maximum non GATC bases (i.e - and N) to allow in pseudogenome sequences [default: 0.5]
>  --rapidnj                           [boolean] Build a tree using the RapidNJ neighbour-joining algorithm
>  --fasttree                          [boolean] Build a tree using the FastTree approximate ML algorithm
>  --iqtree                            [boolean] Build a tree using the IQ-TREE ML algorithm
>  --raxmlng                           [boolean] Build a tree using the RAxML-NG ML algorithm
>
>Generic options
>  --enable_conda                      [boolean] enable conda rather than use containers
>  --validate_params                   [boolean] Boolean whether to validate parameters against the schema at runtime [default: true]
>  --show_hidden_params                [boolean] Show all params when using `--help`
>```

<br>

## Running pipelines with test profile

*   The `test` config profile specifies URLs for test data and all required parameters. With this, you can test any nf-core pipeline with the following command:

```bash
nextflow run nf-core/bactmap -r 1.0.0 -profile test,singularity
```

---

>**If you get an error saying `Failed to pull singularity image` with `status:255` follow the instructions below to create a `/tmp` in your >`$HOME ` folder and re-run the `nextflow run` command**
>
>```bash
>mkdir -p $HOME/tmp && export TMPDIR=$HOME/tmp
>```

---

sample output should looks like this:

```bash
executor >  local (47)[da/e5e356] process > NFCORE_BACTMAP:BACTMAP:INPUT_CHECK:SAMPLESHEET_CHECK (samplesheet.csv)                     [100%] 1 of 1 ✔
[c2/4ccc6f] process > NFCORE_BACTMAP:BACTMAP:BWA_INDEX (NCTC13799.fna)                                           [100%] 1 of 1 ✔[cb/62649d] process > NFCORE_BACTMAP:BACTMAP:FASTP (ERR2172267)                                                  [100%] 3 of 3 ✔
[ab/223a00] process > NFCORE_BACTMAP:BACTMAP:SUB_SAMPLING:MASH_SKETCH (ERR2172267)                               [100%] 3 of 3 ✔[18/bbd394] process > NFCORE_BACTMAP:BACTMAP:SUB_SAMPLING:RASUSA (ERR2172267)                                    [100%] 3 of 3 ✔[d3/13efe1] process > NFCORE_BACTMAP:BACTMAP:BWA_MEM (ERR2172267)                                                [100%] 3 of 3 ✔[3a/be837e] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:SAMTOOLS_SORT (ERR2172267)                        [100%] 3 of 3 ✔[a0/0dde98] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:SAMTOOLS_INDEX (ERR2172265)                       [100%] 3 of 3 ✔[15/f0dd04] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS (ERR2172267)    [100%] 3 of 3 ✔[ed/e12b16] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT (ERR2172265) [100%] 3 of 3 ✔[35/180be7] process > NFCORE_BACTMAP:BACTMAP:BAM_SORT_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS (ERR2172265) [100%] 3 of 3 ✔
[6f/ef2714] process > NFCORE_BACTMAP:BACTMAP:VARIANTS_BCFTOOLS:BCFTOOLS_MPILEUP (ERR2172267)                     [100%] 3 of 3 ✔[c9/bfd766] process > NFCORE_BACTMAP:BACTMAP:VARIANTS_BCFTOOLS:BCFTOOLS_FILTER (ERR2172267)                      [100%] 3 of 3 ✔
[b5/d502dd] process > NFCORE_BACTMAP:BACTMAP:VCF2PSEUDOGENOME (ERR2172266)                                       [100%] 3 of 3 ✔[b2/da2908] process > NFCORE_BACTMAP:BACTMAP:ALIGNPSEUDOGENOMES                                                  [100%] 1 of 1 ✔[5a/e60932] process > NFCORE_BACTMAP:BACTMAP:SNPSITES                                                            [100%] 1 of 1 ✔
[98/61b886] process > NFCORE_BACTMAP:BACTMAP:GUBBINS                                                             [100%] 1 of 1 ✔[e0/1b85b3] process > NFCORE_BACTMAP:BACTMAP:CREATE_PHYLOGENY:RAPIDNJ                                            [100%] 1 of 1 ✔
[14/9755ef] process > NFCORE_BACTMAP:BACTMAP:CREATE_PHYLOGENY:FASTTREE                                           [100%] 1 of 1 ✔[97/acabae] process > NFCORE_BACTMAP:BACTMAP:CREATE_PHYLOGENY:IQTREE                                             [100%] 1 of 1 ✔[36/e40b37] process > NFCORE_BACTMAP:BACTMAP:CREATE_PHYLOGENY:RAXMLNG                                            [100%] 1 of 1 ✔
[f0/8903f9] process > NFCORE_BACTMAP:BACTMAP:GET_SOFTWARE_VERSIONS                                               [100%] 1 of 1 ✔[2c/4fded7] process > NFCORE_BACTMAP:BACTMAP:MULTIQC (1)                                                         [100%] 1 of 1 ✔
-[nf-core/bactmap] Pipeline completed successfully-Completed at: 14-Oct-2021 19:37:51
Duration    : 4m 49s
CPU hours   : 0.2
Succeeded   : 47
```

 > Running with the test profile is a great way to confirm that you have Nextflow configured properly for your system before attempting to run with real data

The `results` folder will contain:

```bash
results/
├── bwa
├── fastp
├── fasttree
├── gubbins
├── iqtree
├── mash
├── multiqc
├── pipeline_info
├── pseudogenomes
├── rapidnj
├── rasusa
├── raxmlng
├── samtools
├── snpsites
└── variants
```

---

## Using nf-core pipelines offline

* Many of the techniques and resources described above require an active internet connection at run time - pipeline files, configuration profiles and software containers are all dynamically fetched when the pipeline is launched. This can be a problem for people using secure computing resources that do not have connections to the internet.

* To help with this, the `nf-core download` command automates the fetching of required files for running nf-core pipelines offline. The command can download a specific release of a pipeline with `-r/--release`. By default, the pipeline will download the pipeline code and the institutional nf-core/configs files.

```bash
nf-core download nf-core/bactmap -r 1.0.0
```

The command prompt will ask whether to download the singularity containers or none. Choose **`singularity`**

```bash
? Download software container images: (Use arrow keys)
   none
 » singularity
```

Next step choose **`n`** to define a shared singualrity cache dir

```bash
? Define $NXF_SINGULARITY_CACHEDIR for a shared Singularity image download folder? [y/n]: n
```

Next for the compression step, choose **`none`**

```bash
? Choose compression type: (Use arrow keys)
 » none
   tar.gz
   tar.bz2
   zip
```

```bash

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.1




In addition to the pipeline code, this tool can download software containers.
? Download software container images: singularity

Nextflow and nf-core can use an environment variable called $NXF_SINGULARITY_CACHEDIR that is a path to a directory where remote Singularity images are stored. This allows downloaded 
images to be cached in a central location.
? Define $NXF_SINGULARITY_CACHEDIR for a shared Singularity image download folder? [y/n]: n

If transferring the downloaded files to another system, it can be convenient to have everything compressed in a single file.
This is not recommended when downloading Singularity images, as it can take a long time and saves very little space.
? Choose compression type: tar.gz
INFO     Saving 'nf-core/bactmap'                                                                                                                                           download.py:160
          Pipeline revision: '1.0.0'                                                                                                                                                       
          Pull containers: 'singularity'                                                                                                                                                   
          Output file: 'nf-core-bactmap-1.0.0.tar.gz'                                                                                                                                      
INFO     Downloading workflow files from GitHub                                                                                                                             download.py:163
INFO     Downloading centralised configs from GitHub                                                                                                                        download.py:167
INFO     Found 18 containers                                                                                                                                                download.py:470
Pulling singularity images ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% • 18/18 completed
INFO     Compressing download..                                                                                                                                             download.py:187
INFO     Command to extract files: tar -xzf nf-core-bactmap-1.0.0.tar.gz                                                                                                    download.py:754
INFO     MD5 checksum for 'nf-core-bactmap-1.0.0.tar.gz': 39f6a27772c82bec184ceb7cd3f2697f                                                                                  download.py:795
```

Now you can run the nf-core/rnaseq pipeline with the provided test data using the test profile

```bash
cd nf-core-bactmap-1.0.0/workflow/
nextflow run main.nf -profile test,singularity
```

<br>

## Troubleshooting

If you run into issues running your pipeline you can you the nf-core website to [troubleshoot common mistakes and issues](https://nf-co.re/usage/troubleshooting).

## Extra resources and getting help

* If you still have an issue with running the pipeline then feel free to contact the nf-core community via the Slack channel . The nf-core Slack organisation has channels dedicated for each pipeline, as well as specific topics (eg. #help, #pipelines, #tools, #configs and much more). 

* The nf-core Slack can be found at `https://nfcore.slack.com` (NB: no hyphen in nfcore!). To join you will need an invite, which you can get at `https://nf-co.re/join/slack`.

* You can also get help by opening an issue in the respective pipeline repository on GitHub asking for help.

* If you have problems that are directly related to Nextflow and not our pipelines or the nf-core framework tools then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).


---
> Citation: If you use nf-core tools in your work, please cite the nf-core publication as follows:
>> The nf-core framework for community-curated bioinformatics pipelines. Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen. Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x. ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

## Quick Recap
>* nf-core is a community-led project to develop a set of best-practice pipelines built using the Nextflow workflow management system.
>* nf-core tools is a suite of helper tools that aims to help people run and develop pipelines.
>* nf-core pipelines can found using the nf-core helper tool `–list` option or from the nf-core website.
>* nf-core pipelines can be run using `nextflow run nf-core/<pipeline>` syntax, or launched and parameters configured using the nf-core helper tool launch option.

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_intro" style="float: left"><b>Back to:</b>Introduction</a>

<a href="/nextflow_varcal/nextflow/nextflow_slurm" style="float: right"><b>Next:</b>NF-Core @ HPC</a></h5>
