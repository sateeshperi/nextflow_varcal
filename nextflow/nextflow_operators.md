---
layout: main
title: Nextflow Operators
categories: [nextflow]
tags: [cluster,nextflow,workflow,bioinformatics,tutorial]
permalink: /nextflow/nextflow_operators
---
{% include _nextflow_nextflow_operators_toc.html %}


<hr>
<center>This is part 8 of 14 of a <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Operators

In the Channels episode we learnt how to create Nextflow channels to enable us to pass data and values around our workflow. If we want to modify the contents or behaviour of a channel, Nextflow provides methods called `operators`. We have previously used the `.view()` operator to view the contents of a channel. There are many more operator methods that can be applied to Nextflow channels that can be usefully separated into several groups:

*   **Filtering** operators: reduce the number of elements in a channel.
*   **Transforming** operators: transform the value/data in a channel.
*   **Splitting** operators: Split items in a channels into smaller chunks.
*   **Combining** operators: join channel together.
*   **Forking** operators: split a single channel into multiple channels.
*   **Maths** operators: apply simple math function on channels.
*   **Other**: Such as the view operator.

In this episode you will see examples, and get to use different types of operators.

## Using Operators

To use an operator, the syntax is the channel name, followed by a dot `.` , followed by the operator name and brackets `()`.

```groovy
channel_obj.<operator>()
```

### view

The `view` operator prints the items emitted by a channel to the console appending a new line character to each item in the channel.

```groovy
ch = channel.of('1', '2', '3')
ch.view()
```
We can also chain together the channel factory method `.of` and the operator `.view()` using the dot notation.

```groovy
ch = channel.of('1', '2', '3').view()
```

To make code more readable we can split the operators over several lines. The blank space between the operators is ignored and is solely for readability.
```groovy
ch = channel
      .of('1', '2', '3')
      .view()
```

>Output
>```groovy
>1
>2
>3
>```

### Closures

An optional closure `{}` parameter can be specified to customise how items are printed.

Briefly, a closure is a block of code that can be passed as an argument to a function. In this way you can define a chunk of code and then pass it around as if it were a string or an integer. By default the parameters for a closure are specified with the groovy keyword `$it` (it is for item).

For example here we use the the `view` operator and apply a closure to it, to add a `chr` prefix to each element of the channel using string interpolation.

```groovy
ch = channel
  .of('1', '2', '3')
  .view({ "chr$it" })
```

>Output
>```groovy
>chr1
>chr2
>chr3
>```
> Note: the `view()` operator doesn’t change the contents of the channel object.


```groovy
ch = channel
  .of('1', '2', '3')
  .view({ "chr$it" })

ch.view()  
```

>Output
>```groovy
>chr1
>chr2
>chr3
>1
>2
>3
>```

## Filtering operators

We can reduce the number of items in a channel by using filtering operators.

The `filter` operator allows you to get only the items emitted by a channel that satisfy a condition and discarding all the others. The filtering condition can be specified by using either:

*   regular expression,
*   a literal value,
*   a data type qualifier, e.g. Number (any integer,float …), String, Boolean
*   or any boolean statement.

**Data type qualifier**

Here we use the `filter` operator on the `chr_ch` channel specifying the data type qualifier `Number` so that only numeric items are returned. The Number data type includes both integers and floating point numbers. We will then use the `view` operator to print the contents.

```groovy
chr_ch = channel.of( 1..22, 'X', 'Y' )
autosomes_ch =chr_ch.filter( Number )
autosomes_ch.view()
```

To simplify the code we can chained together multiple operators, such as `filter` and `view` using a `.`.

The previous example could be rewritten like: The blank space between the operators is ignored and is used for readability.

```groovy
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter( Number )
  .view()
```

>Output
>```groovy
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
>```

**Regular expression**

To filter by a regular expression you have to do is to put ~ right in front of the string literal regular expression (e.g. `~"(^[Nn]extflow)`" or use slashy strings which replace the quotes with `/`. `~/^[Nn]extflow/`).

The following example shows how to filter a channel by using a regular expression `~/^1.*/` inside a slashy string, that returns only strings that begin with 1:

```groovy
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter(~/^1.*/)
  .view()
```

>Output
>```groovy
>1
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
>```

**Boolean statement**

A filtering condition can be defined by using a Boolean expression described by a closure `{}` and returning a boolean value. For example the following fragment shows how to combine a filter for a type qualifier `Number` with another filter operator using a Boolean expression to emit numbers less than 5:

```groovy
channel
  .of( 1..22, 'X', 'Y' )
  .filter(Number)
  .filter { it < 5 }
  .view()
```

>Output
>```groovy
>1
>2
>3
>4
>```

> Closures - In the above example we could remove the brackets around the filter condition e.g. `filter{ it<5}`, since it specifies a closure as the operator’s argument. This is language short for `filter({ it<5})`

**literal value**

Finally, if we only want to include elements of a specific value we can specify a literal value. In the example below we use the literal value `X` to filter the channel for only those elements containing the value `X`.

```groovy
channel
  .of( 1..22, 'X', 'Y' )
  .filter('X')
  .view()
```

>Output
>```groovy
>X
>```

## Modifying the contents of a channel

If we want to modify the items in a channel we use transforming operators.

### map - Applying a function to items in a channel

The `map` operator applies a function of your choosing to every item in a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a closure `{}` as shown in the example below:

```groovy
chr = channel
  .of( 'chr1', 'chr2' )
  .map ({ it.replaceAll("chr","") })

chr.view()
```

Here the map function uses the groovy string function `replaceAll` to remove the chr prefix from each element.

>Output
>```groovy
>1
>2
>```

We can also use the `map` operator to transform each element into a tuple.

In the example below we use the `map` operator to transform a channel containing fastq files to a new channel containing a tuple with the fastq file and the number of reads in the fastq file. We use the `countFastq` file method to count the number of records in a FASTQ formatted file.

We can change the default name of the closure parameter keyword from `it` to a more meaningful name `file` using `->`. When we have multiple parameters we can specify the keywords at the start of the closure, e.g. `file, numreads ->`.

```groovy
fq_ch = channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/*.fastq.gz")
               .map ({ file -> [file, file.countFastq()] })
               .view ({ file, numreads -> "file $file contains $numreads reads" })
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `oper.nf` [happy_plateau] - revision: 926b19594e
>file nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz contains 1553259 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz contains 1553259 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz contains 2768398 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz contains 2768398 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz contains 1107090 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz contains 1107090 reads
>```

We can then add a `filter` operator to only retain thoses fastq files with more than 25000 reads.

```groovy
channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/*.fastq.gz")
       .map ({ file -> [file, file.countFastq()] })
       .filter({ file, numreads -> numreads > 25000})
       .view ({ file, numreads -> "file $file contains $numreads reads" })
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `oper.nf` [cheeky_torvalds] - revision: fc32bf9f3f
>file nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz contains 1553259 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz contains 1553259 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz contains 2768398 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz contains 2768398 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz contains 1107090 reads
>file nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz contains 1107090 reads
>```

### Converting a list into multiple items

The `flatten` operator transforms a channel in such a way that every item in a `list` or `tuple` is flattened so that each single entry is emitted as a sole element by the resulting channel.

```groovy
list1 = [1,2,3]
ch = channel.of(list1)
            .view()
```

>Output
>```groovy
>[1, 2, 3]
>```

```groovy
ch =channel.of(list1)
           .flatten()
           .view()
```

>Output
>```groovy
>1
>2
>3
>```
This is similar to the channel factory `Channel.fromList`.

### Converting the contents of a channel to a single list item.

The reverse of the `flatten` operator is `collect`. The `collect` operator collects all the items emitted by a channel to a list and return the resulting object as a sole emission. This can be extremely useful when combing the results from the output of multiple processes, or a single process run multiple times.

```groovy
ch = channel.of(1, 2, 3, 4)
            .collect()
            .view()
```

It prints a single value:

>Output
>```groovy
>[1,2,3,4]
>```
The result of the collect operator is a `value` channel and can be used multiple times.

### Grouping contents of a channel by a key

The `groupTuple` operator collects `tuples` or `lists` of values by grouping together the channel elements that share the same key. Finally it emits a new tuple object for each distinct key collected.

For example:

```groovy
ch = channel.of( ['wt','wt_1.fq'], ['wt','wt_2.fq'], ["mut",'mut_1.fq'], ['mut', 'mut_2.fq'] )
            .groupTuple()
            .view()
```

>Output
>```groovy
>[wt, [wt_1.fq, wt_1.fq]]
>[mut, [mut_1.fq, mut_2.fq]]
>```

If we know the number of items to be grouped we can use the `groupTuple` `size` parameter. When the specified size is reached, the tuple is emitted. By default incomplete tuples (i.e. with less than size grouped items) are discarded (default).

For example:

```groovy
ch = channel.of( ['wt','wt_1.fq'], ['wt','wt_1.fq'], ["mut",'mut_1.fq'])
            .groupTuple(size:2)
            .view()
```

>Output
>```groovy
>[wt, [wt_1.fq, wt_1.fq]]
>```

This operator is useful to process altogether all elements for which there’s a common property or a grouping key.

### Merging Channels

Combining operators allow you to merge channels together. This can be useful when you want to combine the output channels from multiple processes to perform another task such as joint QC.

**mix**

The `mix` operator combines the items emitted by two (or more) channels into a single channel.

```groovy
ch1 = channel.of( 1,2,3 )
ch2 = channel.of( 'X','Y' )
ch3 = channel.of( 'mt' )

ch4 = ch1.mix(ch2,ch3).view()
```

>Output
>```groovy
>1
>2
>3
>X
>Y
>mt
>```

The items emitted by the resulting mixed channel may appear in any order, regardless of which source channel they came from. Thus, the following example it could be a possible result of the above example as well.

>Output
>```groovy
>1
>2
>X
>3
>mt
>Y
>```

**join**

The `join` operator creates a channel that joins together the items emitted by two channels for which exits a matching key. The key is defined, by default, as the first element in each item emitted.

```groovy
reads1_ch = channel.of(['wt', 'wt_1.fq'], ['mut','mut_1.fq'])
reads2_ch = channel.of(['wt', 'wt_2.fq'], ['mut','mut_2.fq'])
reads_ch = reads1_ch.join(reads2_ch)
                    .view()
```

>Output
>```groovy
>[wt, wt_1.fq, wt_2.fq]
>[mut, mut_1.fq, mut_2.fq]
>```

### Forking operators

Forking operators split a single channel into multiple channels.

**into**

The `into` operator connects a source channel to two or more target channels in such a way the values emitted by the source channel are copied to the target channels. Channel names are separated by a semi-colon. For example:

```groovy
channel.of( 'chr1', 'chr2', 'chr3' )
       .into({ ch1; ch2 })

ch1.view({"ch1 emits: $it"})
ch2.view({"ch2 emits: $it"})
```

>Output
>```groovy
>ch1 emits: chr1
>ch1 emits: chr2
>ch2 emits: chr1
>ch1 emits: chr3
>ch2 emits: chr2
>ch2 emits: chr3
>```

### Maths operators

The maths operators allows you to apply simple math function on channels. The maths operators are:

* count
* min
* max
* sum
* toInteger

**Counting items in a channel**

The `count` operator creates a channel that emits a single item: a number that represents the total number of items emitted by the source channel. For example:

```groovy
ch = channel.of(1..22,'X','Y')
            .count()
            .view()
```

>Output
>```groovy
>24
>```

### Splitting items in a channel

Sometimes you want to split the content of a individual item in a channel, like a file or string, into smaller chunks that can be processed by downstream operators or processes e.g. items stored in a CSV file.

Nextflow has a number of splitting operators that can achieve this:

|splitting operators | function |
| --- | --- |
| [splitCsv](https://www.nextflow.io/docs/latest/operator.html#splitcsv)	| The splitCsv operator allows you to parse text items emitted by a channel, that are formatted using the CSV format, and split them into records or group them into list of records with a specified length. |
| [splitFasta](https://www.nextflow.io/docs/latest/operator.html#splitfasta) |	The splitFasta operator allows you to split the entries emitted by a channel, that are formatted using the FASTA format. It returns a channel which emits text item for each sequence in the received FASTA content. |
| [splitFastq](https://www.nextflow.io/docs/latest/operator.html#splitfastq) |	The splitFastq operator allows you to split the entries emitted by a channel, that are formatted using the FASTQ format. It returns a channel which emits a text chunk for each sequence in the received item. |
| [splitText](https://www.nextflow.io/docs/latest/operator.html#splittext) |	The splitText operator allows you to split multi-line strings or text file items, emitted by a source channel into chunks containing n lines, which will be emitted by the resulting channel. |

**splitCsv**

**FIXME**

The `splitCsv` operator allows you to parse text items emitted by a channel, that are formatted using the CSV format, and split them into records or group them into list of records with a specified length. This is useful when you want to use a sample sheet.

In the simplest case just apply the `splitCsv` operator to a channel emitting a CSV formatted text files or text entries. For example:

Create a CSV file `samples.csv` of all files in `untrimmed_fastq` folder.

```bash
code /workspace/nextflow_tutorial/data/untrimmed_fastq/samples.csv
```

and paste the following and save:

```bash
sample_id,fastq_1,fastq_2
SRR2584863,nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz,nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz
SRR2584866,nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz,nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz
SRR2589044,nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz,nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz
```

We can use the `splitCsv()` operator to split the channel contaning a CSV file into three elements.

```groovy
csv_ch = channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/samples.csv")
                .splitCsv()

csv_ch.view()
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `oper.nf` [peaceful_kalman] - revision: 8c2443bc38
>[sample_id, fastq_1, fastq_2]
>[SRR2584863, nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz]
>[SRR2584866, nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz]
>[SRR2589044, nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz]
>```
The above example shows hows the CSV file `samples.csv` is parsed and is split into three elements.

**Accessing values**

Values can be accessed by its positional index using the square brackets syntax `[index]`. So to access the first column you would use `[0]` as shown in the following example:

```groovy
csv_ch=channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/samples.csv")
              .splitCsv()

csv_ch.view({it[0]})
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `oper.nf` [goofy_cuvier] - revision: bfdd56ff2f
>sample_id
>SRR2584863
>SRR2584866
>SRR2589044
>```

**Column headers**

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its name, as shown in the following example:

```groovy
csv_ch=channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/samples.csv")
              .splitCsv(header:true)

csv_ch.view({it.fastq_1})
```

>Output
>```groovy
>N E X T F L O W  ~  version 21.04.3
>Launching `oper.nf` [big_galileo] - revision: c2d5776695
>nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz
>nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz
>nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz
>```

### Tab delimited files

If you want to split a `tab` delimited file or file separated by another character use the `sep` parameter of the split `splitCsv` operator.

For example:

```groovy
Channel.of("val1\tval2\tval3\nval4\tval5\tval6\n")
       .splitCsv(sep: "\t")
       .view()
```

>Output
>```groovy
>[val1, val2, val3]
>[val4, val5, val6]
>```

## More resources

> See the [operators documentation](https://www.nextflow.io/docs/latest/operator.html) on Nextflow web site for more details

>Quick Recap
>*  Nextflow operators are methods that allow you to modify, set or view channels.
>*  Operators can be separated in to several groups; filtering , transforming , splitting , combining , forking and Maths operators
>*  To use an operator use the dot notation after the Channel object e.g. `my_ch.view()`.
>*  You can parse text items emitted by a channel, that are formatted using the CSV format, using the `splitCsv` operator.

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_workflow" style="float: left"><b>Back to:</b>Nextflow Workflow</a>

<a href="/nextflow_varcal/nextflow/nextflow_variant_calling" style="float: right"><b>Next:</b>NextFlow Variant-Calling Pipeline</a></h5>
