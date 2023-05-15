---
layout: main
title: Nextflow Operators
categories: [nextflow]
tags: [cluster, nextflow, workflow, bioinformatics, tutorial]
permalink: /nextflow/nextflow_operators
---

{% include _nextflow_nextflow_operators_toc.html %}

<hr>
<center>This is part 8 of 14 of a <a href="/nextflow_varcal/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

## Operators in Nextflow

Operators in Nextflow are methods that allow us to modify the contents or behavior of channels. Channels, as we learned in the previous section, are used to pass data and values around our workflow.

Operators can be categorized into several groups:

1. **Filtering Operators**: These reduce the number of elements in a channel.
2. **Transforming Operators**: These transform the value/data in a channel.
3. **Splitting Operators**: These split items in a channel into smaller chunks.
4. **Combining Operators**: These join channels together.
5. **Forking Operators**: These split a single channel into multiple channels.
6. **Math Operators**: These apply simple math functions on channels.
7. **Other Operators**: Such as the view operator.

Here is a diagram to illustrate how operators work in Nextflow:

### Using Operators

The syntax to use an operator is as follows: `channel_obj.<operator>()`.

For example, the `view` operator prints the items emitted by a channel to the console, appending a new line character to each item in the channel. Here's how you can use it:

```groovy
ch = channel.of('1', '2', '3')
ch.view()
```

You can also chain together the channel factory method `.of` and the operator `.view()` using the dot notation:

```groovy
ch = channel.of('1', '2', '3').view()
```

For readability, you can split the operators over several lines. The blank space between the operators is ignored:

```groovy
ch = channel
      .of('1', '2', '3')
      .view()
```

This will output:

```groovy
1
2
3
```

In the following sections, we will explore more examples and use cases of different types of operators.

### Closures

Closures in Groovy are blocks of code that can be defined and passed around as if they were an object. They can be used as parameters in methods, which is particularly useful when using operators in Nextflow. The default parameter for a closure is `$it`, which refers to the current item.

For instance, you can use a closure with the `view` operator to customize how items are printed. Here's an example where we add a `chr` prefix to each element of the channel:

```groovy
ch = channel
  .of('1', '2', '3')
  .view({ "chr$it" })
```

This will output:

```groovy
chr1
chr2
chr3
```

Note that the `view()` operator doesn’t change the contents of the channel object. If you call `view()` again on the same channel, it will print the original items:

```groovy
ch = channel
  .of('1', '2', '3')
  .view({ "chr$it" })

ch.view()
```

This will output:

```groovy
chr1
chr2
chr3
1
2
3
```

## Filtering Operators

Filtering operators are used to reduce the number of items in a channel. The `filter` operator, for example, allows you to get only the items emitted by a channel that satisfy a certain condition.

The filtering condition can be specified using a regular expression, a literal value, a data type qualifier (like Number, String, Boolean), or any boolean statement.

Here's an example where we use the `filter` operator to get only the numeric items from a channel:

```groovy
chr_ch = channel.of( 1..22, 'X', 'Y' )
autosomes_ch = chr_ch.filter( Number )
autosomes_ch.view()
```

This will output:

```groovy
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
```

For readability, you can chain together multiple operators using the dot notation:

```groovy
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter( Number )
  .view()
```

This will give the same output as the previous example.

### Regular Expression

Regular expressions can be used to filter items in a channel. By placing a `~` in front of the string literal regular expression, you can filter items that match the pattern. For example, the regular expression `~/^1.*/` will return only strings that begin with 1:

```groovy
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter(~/^1.*/)
  .view()
```

This will output:

```groovy
1
10
11
12
13
14
15
16
17
18
19
```

### Boolean Statement

A filtering condition can also be defined using a Boolean expression within a closure `{}` that returns a boolean value. For instance, you can combine a filter for a type qualifier `Number` with another filter operator using a Boolean expression to emit numbers less than 5:

```groovy
channel
  .of( 1..22, 'X', 'Y' )
  .filter(Number)
  .filter { it < 5 }
  .view()
```

This will output:

```groovy
1
2
3
4
```

Note: In the above example, we could remove the brackets around the filter condition e.g. `filter{ it<5}`, since it specifies a closure as the operator’s argument. This is a shorthand for `filter({ it<5})`.

### Literal Value

If you want to include elements of a specific value, you can specify a literal value. In the example below, we use the literal value `X` to filter the channel for only those elements containing the value `X`.

```groovy
channel
  .of( 1..22, 'X', 'Y' )
  .filter('X')
  .view()
```

This will output:

```groovy
X
```

These examples illustrate how versatile and powerful Nextflow's filtering operators can be. They allow you to manipulate and control the data flowing through your channels with precision.

## Modifying the Contents of a Channel

Transforming operators are used to modify the items in a channel.

### map - Applying a Function to Items in a Channel

The `map` operator applies a function to every item in a channel and returns a new channel with the modified items. The function is expressed within a closure `{}`. For example:

```groovy
chr = channel
  .of( 'chr1', 'chr2' )
  .map ({ it.replaceAll("chr","") })

chr.view()
```

In this example, the `map` function uses the Groovy string function `replaceAll` to remove the "chr" prefix from each element. The output will be:

```groovy
1
2
```

The `map` operator can also transform each element into a tuple. In the following example, we transform a channel containing fastq files into a new channel containing a tuple with the fastq file and the number of reads in the fastq file:

```groovy
fq_ch = channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/*.fastq.gz")
               .map ({ file -> [file, file.countFastq()] })
               .view ({ file, numreads -> "file $file contains $numreads reads" })
```

We can add a `filter` operator to retain only those fastq files with more than 25000 reads:

```groovy
channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/*.fastq.gz")
       .map ({ file -> [file, file.countFastq()] })
       .filter({ file, numreads -> numreads > 25000})
       .view ({ file, numreads -> "file $file contains $numreads reads" })
```

### Converting a List into Multiple Items

The `flatten` operator transforms a channel so that every item in a list or tuple is emitted as a separate element by the resulting channel:

```groovy
list1 = [1,2,3]
ch = channel.of(list1)
            .flatten()
            .view()
```

This will output:

```groovy
1
2
3
```

This is similar to the channel factory `Channel.fromList`.

### Converting the Contents of a Channel to a Single List Item

The `collect` operator collects all the items emitted by a channel into a list and returns the resulting object as a single emission. This can be useful when combining the results from multiple processes, or a single process run multiple times:

```groovy
ch = channel.of(1, 2, 3, 4)
            .collect()
            .view()
```

This will output:

```groovy
[1,2,3,4]
```

The result of the `collect` operator is a `value` channel and can be used multiple times.

### Grouping Channel Contents by a Key in Nextflow

In Nextflow, the `groupTuple` operator is utilized to group tuples or lists of values in a channel based on a shared key. This operator then emits a new tuple for each distinct key collected.

For example, consider the following Nextflow syntax:

```groovy
ch = channel.of( ['wt','wt_1.fq'], ['wt','wt_2.fq'], ['mut','mut_1.fq'], ['mut', 'mut_2.fq'] )
            .groupTuple()
            .view()
```

The output produced would be:

```groovy
[wt, [wt_1.fq, wt_2.fq]]
[mut, [mut_1.fq, mut_2.fq]]
```

If the number of items to be grouped is known, the `groupTuple` `size` parameter can be used. Once the specified size is reached, the tuple is emitted. By default, incomplete tuples (i.e., with fewer than the specified size of grouped items) are discarded.

For instance, consider the following example:

```groovy
ch = channel.of( ['wt','wt_1.fq'], ['wt','wt_2.fq'], ['mut','mut_1.fq'])
            .groupTuple(size:2)
            .view()
```

The output in this case would be:

```groovy
[wt, [wt_1.fq, wt_2.fq]]
```

The `groupTuple` operator is particularly helpful when processing elements with a common property or grouping key.

### Merging Channels in Nextflow

Nextflow allows you to merge channels using combining operators, which can be useful when you need to combine output channels from multiple processes to perform tasks like joint QC.

**mix**

The `mix` operator combines items emitted by two or more channels into a single channel. Consider the following example:

```groovy
ch1 = channel.of( 1, 2, 3 )
ch2 = channel.of( 'X', 'Y' )
ch3 = channel.of( 'mt' )

ch4 = ch1.mix(ch2, ch3).view()
```

The output could be:

```groovy
1
2
3
X
Y
mt
```

The items emitted by the resulting mixed channel may appear in any order, regardless of which source channel they originated from. Therefore, the following output could also be a possible result of the above example:

```groovy
1
2
X
3
mt
Y
```

In summary, Nextflow provides the `groupTuple` operator for grouping channel elements by a key, and the `mix` operator for merging channels. These operators can be helpful when processing data with shared properties or combining outputs from multiple channels.

### Joining Channels in Nextflow

The `join` operator in Nextflow allows you to create a channel that joins items emitted by two channels based on a matching key. By default, the key is defined as the first element in each item emitted.

For example:

```groovy
reads1_ch = channel.of(['wt', 'wt_1.fq'], ['mut','mut_1.fq'])
reads2_ch = channel.of(['wt', 'wt_2.fq'], ['mut','mut_2.fq'])
reads_ch = reads1_ch.join(reads2_ch)
                    .view()
```

Output:

```groovy
[wt, wt_1.fq, wt_2.fq]
[mut, mut_1.fq, mut_2.fq]
```

### Forking Operators in Nextflow

Forking operators enable you to split a single channel into multiple channels.

**into**

The `into` operator connects a source channel to two or more target channels, copying values emitted by the source channel to the target channels. Channel names are separated by a semi-colon. For example:

```groovy
channel.of( 'chr1', 'chr2', 'chr3' )
       .into({ ch1; ch2 })

ch1.view({"ch1 emits: $it"})
ch2.view({"ch2 emits: $it"})
```

Output:

```groovy
ch1 emits: chr1
ch1 emits: chr2
ch2 emits: chr1
ch1 emits: chr3
ch2 emits: chr2
ch2 emits: chr3
```

### Maths Operators in Nextflow

Maths operators allow you to apply simple math functions on channels. These operators include:

- count
- min
- max
- sum
- toInteger

**Counting items in a channel**

The `count` operator creates a channel that emits a single item: a number representing the total number of items emitted by the source channel. For example:

```groovy
ch = channel.of(1..22,'X','Y')
            .count()
            .view()
```

Output:

```groovy
24
```

### Splitting Items in a Channel in Nextflow

In some cases, you may want to split the content of an individual item in a channel, such as a file or string, into smaller chunks that can be processed by downstream operators or processes (e.g., items stored in a CSV file).

Nextflow offers several splitting operators for this purpose:

| Splitting Operators                                                        | Function                                                                                                                                                                    |
| -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [splitCsv](https://www.nextflow.io/docs/latest/operator.html#splitcsv)     | The splitCsv operator parses text items emitted by a channel in CSV format and splits them into records or groups them into lists of records with a specified length.       |
| [splitFasta](https://www.nextflow.io/docs/latest/operator.html#splitfasta) | The splitFasta operator splits entries emitted by a channel in FASTA format. It returns a channel that emits a text item for each sequence in the received FASTA content.   |
| [splitFastq](https://www.nextflow.io/docs/latest/operator.html#splitfastq) | The splitFastq operator splits entries emitted by a channel in FASTQ format. It returns a channel that emits a text chunk for each sequence in the received item.           |
| [splitText](https://www.nextflow.io/docs/latest/operator.html#splittext)   | The splitText operator splits multi-line strings or text file items emitted by a source channel into chunks containing n lines, which are emitted by the resulting channel. |

In summary, Nextflow provides a variety of operators for joining, forking, applying.

### splitCsv

The `splitCsv` operator in Nextflow allows you to parse text items emitted by a channel that are formatted using the CSV format and split them into records or group them into lists of records with a specified length. This is useful when working with sample sheets, for example.

In the simplest case, just apply the `splitCsv` operator to a channel emitting CSV-formatted text files or text entries. Consider the following example:

Create a CSV file `samples.csv` containing all files in the `untrimmed_fastq` folder:

```bash
code /workspace/nextflow_tutorial/data/untrimmed_fastq/samples.csv
```

Paste the following and save:

```bash
sample_id,fastq_1,fastq_2
SRR2584863,nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz,nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz
SRR2584866,nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz,nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz
SRR2589044,nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz,nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz
```

You can use the `splitCsv()` operator to split the channel containing a CSV file into three elements:

```groovy
csv_ch = channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/samples.csv")
                .splitCsv()

csv_ch.view()
```

Output:

```groovy
[sample_id, fastq_1, fastq_2]
[SRR2584863, nextflow_tutorial/data/untrimmed_fastq/SRR2584863_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2584863_2.fastq.gz]
[SRR2584866, nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2584866_2.fastq.gz]
[SRR2589044, nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz, nextflow_tutorial/data/untrimmed_fastq/SRR2589044_2.fastq.gz]
```

This example demonstrates how the CSV file `samples.csv` is parsed and split into three elements.

**Accessing values**

Values can be accessed by their positional index using the square brackets syntax `[index]`. To access the first column, you would use `[0]` as shown in the following example:

```groovy
csv_ch = channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/samples.csv")
                .splitCsv()

csv_ch.view({it[0]})
```

Output:

```groovy
sample_id
SRR2584863
SRR2584866
SRR2589044
```

**Column headers**

When the CSV file begins with a header line defining the column names, you can specify the parameter `header: true`, which allows you to reference each value by its name, as shown in the following example:

```groovy
csv_ch = channel.fromPath("/workspace/nextflow_tutorial/data/untrimmed_fastq/samples.csv")
                .splitCsv(header:true)

csv_ch.view({it.fastq_1})
```

Output:

```groovy
nextflow_tutorial/data/untrimmed_fastq/SRR2584866_1.fastq.gz
nextflow_tutorial/data/untrimmed_fastq/SRR2589044_1.fastq.gz
```

### Tab delimited files

If you want to split a `tab` delimited file or file separated by another character, use the `sep` parameter of the `splitCsv` operator.

For example:

```groovy
Channel.of("val1\tval2\tval3\nval4\tval5\tval6\n")
       .splitCsv(sep: "\t")
       .view()
```

Output:

```groovy
[val1, val2, val3]
[val4, val5, val6]
```

## More resources

- See the [operators documentation](https://www.nextflow.io/docs/latest/operator.html) on the Nextflow website for more details.
- Quick Recap:
  - Nextflow operators are methods that allow you to modify, set, or view channels.
  - Operators can be separated into several groups: filtering, transforming, splitting, combining, forking, and math operators.
  - To use an operator, use the dot notation after the Channel object, e.g., `my_ch.view()`.
  - You can parse text items emitted by a channel that are formatted using the CSV format using the `splitCsv` operator.

---

<h5><a href="/nextflow_varcal/nextflow/nextflow_workflow" style="float: left"><b>Back to:</b>Nextflow Workflow</a>

<a href="/nextflow_varcal/nextflow/nextflow_variant_calling" style="float: right"><b>Next:</b>NextFlow Variant-Calling Pipeline</a></h5>
