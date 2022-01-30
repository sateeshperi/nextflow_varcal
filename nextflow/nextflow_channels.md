---
layout: main
title: Nextflow Channels
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_channels
---
{% include _nextflow_nextflow_channels_toc.html %}


<hr>
<center>This is part 5 of 14 of a <a href="/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

```bash
cd ~/nextflow_tutorial
mkdir channels
cd channels
```

*   **Channels** are how Nextflow handles file management, allowing complex tasks to be split up, run in parallel, and reduces ‘admin’ required to get the right inputs to the right parts of the pipeline.

![](images/channel-files.png)

*   Channels connect processes via their `inputs` and `outputs`.
*   Channels can store multiple items, such as files (e.g., fastq files) or values.
*   The number of items a channel stores determines how many times a process will run using that channel as input. When the process runs using one item from the input channel, we will call that run a task.
*   **Channels are asynchronous**, which means that output data from a set of processes will not necessarily be output in the same order as they went in. 
*   However, the first element into a queue is the first out of the queue (**First in-First out**). This allows processes to run as soon as they receive input from a channel. **Channels only send data in one direction, from a producer (a process/operator), to a consumer (another process/operator)**.
*   Nextflow distinguishes between two different kinds of channels: 
    1. **`value`** channels
    2.  **`queue`** channels

### Value channels

*   **A value channel is bound to a single value**.
*   A value channel can be used an unlimited number times since its content is not consumed. This is also useful for processes that need to reuse input from a channel, for example, a reference genome sequence file that is required by multiple steps within a process, or by more than one process.
*   The `value` factory method is used to create a value channel. Values are put inside parentheses `()` to assign them to a channel.

**Create a nextflow file called `value_ch.nf`; add the following and `nextflow run value_ch.nf`:**

>```groovy
>//Creates an empty value channel.
>ch1 = Channel.value( 'GRCh38' )
>//Creates a value channel and binds a string to it.
>ch2 = Channel.value( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] )
>//Creates a value channel and binds a list object to it that will be emitted as a single item.
>ch3 = Channel.value( ['chr1' : 248956422, 'chr2' : 242193529, 'chr3' : 198295559] )
>//The value method can only take 1 argument, however, this can be a single list containing several elements.
>ch1.view()
>ch2.view()
>ch3.view()
>```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `value_ch.nf` [happy_bell] - revision: a80c27d5b2
>GRCh38
>[chr1, chr2, chr3, chr4, chr5]
>[chr1:248956422, chr2:242193529, chr3:198295559]
>```

>**Reminder:**
>* A [List object](https://www.tutorialspoint.com/groovy/groovy_lists.htm) can be defined by placing the values in square brackets [] separated by a comma.
>* A [Map object](https://www.tutorialspoint.com/groovy/groovy_maps.htm) is similar, but with `key:value pairs` separated by commas.
>
>Example:
>```groovy
>myList = [1776, -1, 33, 99, 0, 928734928763]
>
>myMap = [ p1 : "start", q2 : "end" ]
>```
---

### Queue channel

*   `Queue` channels are a type of channel in which data is consumed (used up) to make input for a process/operator.
*   Queue channels can be created in two ways:
    1.  As the outputs of a process.
    2.  Or a queue channels can be explicitly created using channel factory methods such as **[Channel.of](https://www.nextflow.io/docs/latest/channel.html#of)** or **[Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath)**.

> In Nextflow DSL1 queue channels can only be used once in a workflow, either connecting workflow input to process input, or process output to input for another process. **In DSL2 we can use a queue channel multiple times.**

> **What is the difference between a queue and value channel?**
> *    A queue channels is consumable and can store multiple values.
> *    A value channel a.k.a. singleton channel by definition is bound to a single value

## Queue channel factory

*   Channel factories are used to explicitly create channels. Channel factories are called using the `Channel.<method>` syntax, and return a specific instance of a Channel.

*   `Queue` (**consumable**) channels can be created using the following channel factory methods.
    1.  **[Channel.of](https://www.nextflow.io/docs/edge/channel.html#of)**
    2.  **[Channel.fromList](https://www.nextflow.io/docs/edge/channel.html#fromlist)**
    3.  **[Channel.fromPath](https://www.nextflow.io/docs/edge/channel.html#frompath)**
    4.  **[Channel.fromFilePairs](https://www.nextflow.io/docs/edge/channel.html#fromfilepairs)**
    5.  **[Channel.fromSRA](https://www.nextflow.io/docs/edge/channel.html#fromsra)**

---
<br>

### 1. **The `of` Channel factory**

*   When you want to create a channel containing multiple values you can use the channel factory `Channel.of`.
*   `Channel.of` allows the creation of a `queue` channel with the values specified as arguments, separated by a comma `,`.
*   You can specify a range of numbers as a single argument using the Groovy range operator `..`. This creates each value in the range (including the start and end values) as a value in the channel. The Groovy range operator can also produce ranges of dates, letters, or time. More information on the range operator can be found **[here](https://www.logicbig.com/tutorials/misc/groovy/range-operator.html)**.

**Create a nextflow file called `queue_of.nf`; add the following to create a `Channel.of` method and `nextflow run queue_of.nf`:**

```groovy
chromosome_ch = Channel.of( 'chr1','chr3','chr5','chr7' )
chromosome_ch = Channel.of(1..22, 'X', 'Y')
chromosome_ch.view()
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `queue_of.nf` [modest_miescher] - revision: fb9862a33f
>chr1
>chr3
>chr5
>chr7
>1
>2
>3
>4
>5
>6
>7
>8
>9
>10
>11
>12
>13
>14
>15
>16
>17
>18
>19
>20
>21
>22
>X
>Y
>```
>>*   The first line in this example creates a variable `chromosome_ch`. `chromosome_ch` is a queue channel containing the four values specified as arguments in the `of` method. The `.view` operator will print one line per item in a list. Therefore the view operator on the second line will print four lines, one for each element in the channel.
>>* Arguments passed to the `of` method can be of varying types e.g., combinations of numbers, strings, or objects. In the above examples we have examples of both string and number data types.

---
<br>

### 2. **The `fromList` Channel factory**

*   You can use the `Channel.fromList` method to create a queue channel from a list object.

**Create a new file `queue_list.nf`; add the following and `nextflow run queue_list.nf`:**

```groovy
aligner_list = ['salmon', 'kallisto']

aligner_ch = Channel.fromList(aligner_list)

aligner_ch.view()
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `queue_list.nf` [admiring_raman] - revision: b2dec22cea
>salmon
>kallisto
>```

>**`Channel.fromList` vs `Channel.of`**
>>In the above example, the channel has two elements. If you has used the `Channel.of(list)` it would have contained only 1 element [salmon, kallisto] and any operator or process using the channel would run once.

---
<details>
  <summary><b>CLICK HERE for an exercise on creating channels from a list</b></summary>

Create a Nextflow script `exercise_1.nf` that creates both a `queue` and `value` channel for the following list. Then print the contents of the channels using the view operator. How many lines does the queue and value channel print?

<pre>
ids = ['ERR908507', 'ERR908506', 'ERR908505']
</pre>

Solution
<pre>
ids = ['ERR908507', 'ERR908506', 'ERR908505']
queue_ch = Channel.fromList(ids)
value_ch = Channel.value(ids)
queue_ch.view()
value_ch.view()
</pre>

Output

<pre>
N E X T F L O W  ~  version 21.04.3
Launching `exercise_1.nf` [silly_hypatia] - revision: 5d8381b7d4
ERR908507
ERR908506
ERR908505
[ERR908507, ERR908506, ERR908505]
</pre>
<pre>
The queue channel `queue_ch` will print three lines.
The value channel `value_ch` will print one line.
</pre>
</details>
---
<br>

### 3. **The `fromPath` Channel factory**

*   The previous channel factory methods dealt with sending general values in a channel. A special channel factory method `fromPath` is used when wanting to pass files.

*   The `fromPath` factory method creates a queue channel emitting one or more files matching a file path.
*   The file path (written as a quoted string) can be the location of a single file or a **“glob pattern”** that matches multiple files or directories.
>   A glob pattern is specified as a string and is matched against directory or file names.
>    1.   An asterisk, `*`, matches any number of characters (including none).
>    2.   Two asterisks, `**`, works like `*` but will also search sub directories. This syntax is generally used for matching complete paths.
>    3.   Braces `{}` specify a collection of subpatterns. For example: `{bam,bai}` matches `bam` or `bai`
*   The file path can be a relative path (path to the file from the current directory), or an absolute path (path to the file from the system root directory - starts with `/`).

*   You can change the behaviour of `Channel.fromPath` method by changing its options. A list of `.fromPath` options is shown below:
*   In Nextflow, method parameters are separated by a `,` and parameter values specified with a colon `:`

| Name | Description |
| ---  | ---         |
| **glob** | When true, the characters `*`, `?`, `[]` and `{}` are interpreted as glob wildcards, otherwise they are treated as literal characters (default: true) |
| **type** | The type of file paths matched by the string, either file, dir or any (default: file) |
| **hidden** | When true, hidden files are included in the resulting paths (default: false) |
| **maxDepth** | Maximum number of directory levels to visit (default: no limit) |
| **followLinks** | When true, symbolic links are followed during directory tree traversal, otherwise they are managed as files (default: true) |
| **relative** | When true returned paths are relative to the top-most common directory (default: false) |
| **checkIfExists** | When true throws an exception if the specified path does not exist in the file system (default: false) |

**Create a new file `queue_path.nf`; add the following and `nextflow run queue_path.nf`:**

```groovy
println("Single File")
read_ch = Channel.fromPath("$HOME/nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz")
read_ch.view()

println("Glob Syntax")

// We can change the default options for the `fromPath` method to give an error if the file doesn’t exist using the `checkIfExists` parameter.
read_ch = Channel.fromPath( "$HOME/nextflow_tutorial/data/untrimmed_fastq/*.fastq.gz", checkIfExists: true )
read_ch.view()

// **Note The pattern must contain at least a star wildcard character.**
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `queue_path.nf` [awesome_yalow] - revision: 695fc25ebf
>Single File
>nextflow_tutorial/channels/$HOME/nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz
>Glob Syntax
>nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz
>nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz
>nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz
>nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz
>nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz
>nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz
>```

---

### 4. **The `fromFilePairs` Channel factory**

*   We have seen how to process files individually using `fromPath`. In Bioinformatics we often want to process files in pairs or larger groups, such as read pairs in sequencing.

For example in the `$HOME/nextflow_tutorial/data/untrimmed_fastq` directory we have three groups of read pairs:

| Sample_ID | read1 |	read2 |
| --- | --- | --- |
| SRR2589044 |	SRR2589044_1.fastq.gz |	SRR2589044_2.fastq.gz |
| SRR2584863 |	SRR2584863_1.fastq.gz |	SRR2584863_2.fastq.gz |
| SRR2584866 |	SRR2584866_1.fastq.gz |	SRR2584866_2.fastq.gz |

*   Nextflow provides a convenient factory method for this common bioinformatics use case. The `fromFilePairs` method creates a queue channel containing a **tuple** for every set of files matching a specific pattern (e.g., `/path/to/*_{1,2}.fastq.gz`).

*   **A tuple is a grouping of data, represented as a Groovy List**. 
    1. The first element of the tuple emitted from `fromFilePairs` is a string based on the shared part of the filenames (i.e., the `*` part of the glob pattern).
    2. The second element is the list of files matching the remaining part of the glob pattern (i.e., the `<string>_{1,2}.fastq.gz pattern`). This will include any files ending `_1.fastq.gz` or `_2.fastq.gz`.

**Create a new `queue_pairs.nf`; add the following and `nextflow run queue_pairs.nf`:**

```groovy
read_pair_ch = Channel.fromFilePairs("$HOME/nextflow_tutorial/data/untrimmed_fastq/*_{1,2}.fastq.gz", checkIfExists: true)
read_pair_ch.view()
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `queue_pairs.nf` [nauseous_austin] - revision: eb531c1c8a
>[SRR2589044, [nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz]]
>[SRR2584863, [nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz]]
>[SRR2584866, [nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz]]
>```

This will produce a queue channel, `read_pair_ch` , containing three elements:

*   Each element is a tuple that has:
    1. string value (the file prefix matched, e.g `SRR2584866`)
    2. and a list with the two files e,g. [`nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz`] .

*   The asterisk character `*`, matches any number of characters (including none), and the `{}` braces specify a collection of subpatterns. 
*   Therefore the `*_{1,2}.fastq.gz` pattern matches any file name ending in `_1.fastq.gz` or `_2.fastq.gz`.

> Read more information about the channel factory `fromFilePairs` [here](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)

---

### 5. **The `fromSRA` Channel factory**

*   The `fromSRA` method makes it possible to query the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) archive and returns a queue channel emitting the FASTQ files matching the specified selection criteria.

*   Read more about `fromSRA` Channel factory **[here](https://www.nextflow.io/docs/latest/channel.html#fromsra)**


FIXME: add few examples of complex multi map input for channels

> Quick Recap:
>> * Channels must be used to import data into Nextflow.
>> * Nextflow has two different kinds of channels, `queue` channels and `value` channels.
>> * Data in `value` channels can be used multiple times in workflow.
>> * Data in `queue` channels are consumed when they are used by a process or an operator.
>> * Channel factory methods, such as `Channel.of`, are used to create Channels.
>> * Channel factory methods have optional parameters e.g., `checkIfExists`, that can be used to alter the creation and behaviour of a channel.

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_scripting" style="float: left"><b>Back to:</b>Nextflow Scripting</a>

<a href="/nextflow_varcal/nextflow/nextflow_processes" style="float: right"><b>Next:</b>NextFlow Processes</a></h5>
