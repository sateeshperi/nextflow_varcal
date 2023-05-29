---
layout: main
title: Nextflow Scripting
categories: [nextflow]
tags: [cluster, nextflow, workflow, bioinformatics, tutorial]
permalink: /docs/nextflow/nextflow_scripting
---

{% include _docs_nextflow_nextflow_scripting_toc.html %}

<hr>
<center>This is part 4 of 14 of a <a href="/nextflow_varcal/docs/nextflow/" target="_blank">Introduction to NextFlow</a>.</center>
<hr>

<br>

```bash
mkdir scripting
cd scripting/
```

## Language Basics

- A **Domain Specific Language (DSL)** is a computer language that's targeted to a particular kind of problem, rather than a general-purpose language that's aimed at any kind of software problem.
- Nextflow is a **Domain Specific Language (DSL)** implemented on top of the Groovy programming language, which in turn is a superset of the Java programming language. This means that Nextflow can run any Groovy and Java code.
- **However, it is not necessary to learn Groovy to use Nextflow DSL, but it can be useful in edge cases where you need more functionality than the DSL provides.**

Below are a few basics of the Groovy syntax.

### Printing values

- To print something is as easy as using the `println` method and passing the text to print in quotes.

**Create a new file `groovy.nf` and add the following:**

```groovy
//groovy.nf
println("Hello, World!")
```

```bash
nextflow run groovy.nf
```

> Output
>
> ```bash
> N E X T F L O W  ~  version 21.04.3
> Launching `groovy.nf` [dreamy_jepsen] - revision: 98447b37f7
> Hello, World!
> ```

### Comments

- When we write any code, it is useful to document it using comments. In Nextflow, comments use the same syntax as in the C-family programming languages. This can be confusing for people familiar with the `#` syntax for commenting in other languages.

```groovy
// This is a single line comment. Everything after the // is ignored.

/*
   Comments can also
   span multiple
   lines.
 */
```

### Variables

- In any programming language, you need to use variables to store different types of information. A variable is a pointer to a space in the computer’s memory that stores the value associated with it.

- Variables are assigned using `=` and can have any value. **Groovy is dynamically-typed which means the variable’s data types is based on its value.**

Groovy knows various types of data. Four common ones are:

- `String` − These are text literals which are represented in the form of a chain of characters enclosed in quotes. For example, "Hello World".
- `int` − This is used to represent whole numbers. An example is 1234.
- `Boolean` − This represents a Boolean value which can either be true or false.
- `float` - This is used to represent floating-point numbers, like 12.34.

A more complete list can be found **[here](https://www.tutorialspoint.com/groovy/groovy_data_types.htm)**.

---

**Add the following example data types to `groovy.nf` and run `nextflow run groovy.nf`:**

```groovy
//int − This is used to represent whole numbers.
my_int = 1

//float − This is used to represent floating-point numbers.
my_float = 3.1499392

//Boolean − This represents a Boolean value which can either be true or false.
my_bool = false

//String - These are text literals which are represented in the form of a chain of characters.
my_string = "chr1"

// A block of text that spans multiple lines can be defined by delimiting it with triple single `'''` or double quotes `"""`:
my_text = """
          This is a multi-line
          text using triple quotes.
          """

// To display the value of a variable to the screen in Groovy, we can use the `println` method passing the variable name as a parameter.
println(my_int)
println(my_float)
println(my_bool)
println(my_string)
println(my_text)

// String Interpolation
// To use a variable inside a single or multi-line double-quoted string "", prefix the variable name with a $ to show it should be interpolated.
println("processing chromosome $my_int")
println("value of pi is $my_float")
```

> Output
>
> ```bash
> N E X T F L O W  ~  version 21.04.3
> Launching `groovy.nf` [angry_cajal] - revision: 6ebb8cff41
> Hello, World!
> 1
> 3.1499392
> false
> chr1
>
>          This is a multi-line
>          using triple quotes.
>
> processing chromosome 1
> value of pi is 3.1499392
> ```

> **Variable names inside single quoted strings do not support String interpolation.**

---

### def

- Local variables are defined using the `def` keyword. It should always be used when defining variables local to a function or a closure.

**Add the following line to `groovy.nf` and run `nextflow run groovy.nf`:**

```groovy
def x = 'local_variable_def'
println(x)
```

> Output
>
> ```groovy
> N E X T F L O W  ~  version 21.04.3
> Launching `groovy.nf` [compassionate_marconi] - revision: 932bc80461
> Hello, World!
> 1
> 3.1499392
> false
> chr1
>
>          This is a multi-line
>          text using triple quotes.
>
> processing chromosome 1
> value of pi is 3.1499392
> local_variable_def
> ```

---

### Lists

- To store multiple values in a variable, we can use a List. A List (also known as an array) object can be defined by placing the list items in square brackets and separating items by commas `,`.
- Negative numbers can be used as indices in Groovy. When we do so, the index `-1` gives us the last element in the list, `-2` the second to last, and so on. Because of this, kmers[3] and kmers[-1] point to the same element here.

**Create a new `lists.nf` file; add the following and run `nextflow run lists.nf`:**

```groovy
kmers = [11,21,27,31]
// You can access a given item in the list with square-bracket notation []. These positions are numbered starting at 0, so the first element has an index of 0.
println(kmers[0])
// Lists can also be indexed with negative indexes
println(kmers[-1])
// The first three elements Lists elements using a range.
println(kmers[0..2])
// String interpolation for Lists - To use an expression like `kmer[0..2]` inside a double quoted String `""`, we use the `${expression}` syntax, similar to Bash/shell scripts
println("The first three elements in the Lists are: ${kmers[0..2]}")
```

> Output
>
> ```groovy
> N E X T F L O W  ~  version 21.04.3
> Launching `lists.nf` [intergalactic_minsky] - revision: 068f0e0a7c
> 11
> 31
> [11, 21, 27]
> The first three elements in the Lists are: [11, 21, 27]
> ```

---

### List Methods

- Lists implement a number of useful methods that can perform operations on their contents. See more **[here](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html)**. When using a method on a type of object, you need to prefix the method with the variable name.

**Add the following to `lists.nf` and run `nextflow run lists.nf`:**

```groovy
// To get the length of the list
println(kmers.size())

// Inside a string, we need to use the ${} syntax
println("list size is:  ${kmers.size()}")

// To retrieve items in a list
println(kmers.get(1))
```

> Output
>
> ```groovy
> N E X T F L O W  ~  version 21.04.3
> Launching `lists.nf` [fervent_meninsky] - revision: bf70a3d08a
> 11
> 31
> [11, 21, 27]
> The first three elements in the Lists are: [11, 21, 27]
> 4
> list size is:  4
> 21
> ```

---

<details>
  <summary><b>CLICK HERE for more common list methods</b></summary>

<pre>
```
mylist = [1,2,3]
println mylist
println mylist + [1]
println mylist - [1]
println mylist * 2
println mylist.reverse()
println mylist.collect{ it+3 }
println mylist.unique().size()
println mylist.count(1)
println mylist.min()
println mylist.max()
println mylist.sum()
println mylist.sort()
println mylist.find{it%2 == 0}
println mylist.findAll{it%2 == 0}
```
</pre>

Output

<pre>
```
[1, 2, 3]
[1, 2, 3, 1]
[2, 3]
[1, 2, 3, 1, 2, 3]
[3, 2, 1]
[4, 5, 6]
3
1
1
3
6
[1, 2, 3]
2
[2]
```
</pre>

## </details>

### Maps

- It can be difficult to remember the index of a value in a list, so we can use Groovy Maps (also known as associative arrays) that have an arbitrary type of key instead of an integer value.
- The syntax is very much aligned. To specify the key, use a colon before the value [key:value]. Multiple values are separated by a comma. Note: the key value is not enclosed in quotes.

**Create a new `maps.nf` file; add the following and run `nextflow run maps.nf`:**

```groovy
roi = [ chromosome : "chr17", start: 7640755, end: 7718054, genes: ['ATP1B2','TP53','WRAP53']]

// Maps can be accessed in a conventional square-bracket syntax or as if the key was a property of the map or using the dot notation. **Note: When retrieving a value, the key value is enclosed in quotes.**
println(roi['chromosome'])

// Use a dot notation
println(roi.start)

// Use the get method
println(roi.get('genes'))
```

> Output
>
> ```groovy
> N E X T F L O W  ~  version 21.04.3
> Launching `maps.nf` [pedantic_hamilton] - revision: c289105bc9
> chr17
> 7640755
> [ATP1B2, TP53, WRAP53]
> ```

---

<details>
  <summary><b>CLICK HERE for syntax on adding/modifying a map</b></summary>

<pre>
Add the following to `lists.nf` and run `nextflow run maps.nf`:
```
// Use the square brackets
roi['chromosome'] = 'chr19'

// or

// Use a dot notation        
roi.chromosome = 'chr19'  

// Use the put method              
roi.put('genome', 'hg38')

println("Adding or modifying list results:")
println(roi['chromosome'])
println(roi['genome'])
```
</pre>

Output

<pre>
```groovy
N E X T F L O W  ~  version 21.04.3
Launching `maps.nf` [elegant_lumiere] - revision: c3c0dc5ad4
chr17
7640755
[ATP1B2, TP53, WRAP53]
Adding or modifying list results:
chr19
hg38
```
</pre>

## </details>

### Closures

- **Closures are the swiss army knife of Nextflow/Groovy programming**. In a nutshell a closure is is a block of code that can be passed as an argument to a function. This can be useful to create a re-usable function.

- We can assign a closure to a variable in same way as a value using the `=`.

> ```groovy
> square = { it * it }
> ```

- **The curly brackets `{}` around the expression `it * it` tells the script interpreter to treat this expression as code**. it is an implicit variable that is provided in closures. It’s available when the closure doesn’t have an explicitly declared parameter and represents the value that is passed to the function when it is invoked.

- We can pass the function `square` as an argument to other functions or methods. Some built-in functions take a function like this as an argument. One example is the `collect` method on lists that iterates through each element of the list transforming it into a new value using the closure.

**Create a new `closures.nf` file; add the following as an example and `nextflow run closures.nf`:**

```groovy
square = { it * it }
x = [ 1, 2, 3, 4 ]
println(x)
println(x.collect(square))
```

> Output
>
> ```groovy
> N E X T F L O W  ~  version 21.04.3
> Launching `closures.nf` [furious_mandelbrot] - revision: 73a27ee028
> [1, 2, 3, 4]
> [1, 4, 9, 16]
> ```

#### Closure parameters

- By default, closures take a single parameter called **`it`**. To define a different name, use the `variable ->` syntax in Nextflow.
- In the example below, the variable `num` is assigned as the closure input parameter instead of `it`:

> ```groovy
> square = { num -> num * num }
> ```

Let's define a closure to add the prefix `chr` to each element of the list in a Nextflow script.

**Add the following to `closures.nf` and execute it using `nextflow run closures.nf`:**

```groovy
prefix = {"chr${it}"}
x = x.collect(prefix)
println x
```

> Output:
>
> ```groovy
> N E X T F L O W  ~  version 21.04.3
> Launching `closures.nf` [exotic_swirles] - revision: 016a18e635
> [1, 2, 3, 4]
> [1, 4, 9, 16]
> [chr1, chr2, chr3, chr4]
> ```

---

### Multiple map parameters

- In Nextflow, it’s also possible to define closures with multiple, custom-named parameters using the `->` syntax. Separate the custom-named parameters by a comma before the `->` operator.

**Add the following to `closures.nf` and execute it using `nextflow run closures.nf`:**

```groovy
tp53 = [chromosome: "chr17", start: 7661779, end: 7687538, genome: 'GRCh38', gene: "TP53"]
// Perform subtraction of end and start coordinates
region_length = {start, end -> end - start }
// Add the region length to the map tp53
tp53.length = region_length(tp53.start, tp53.end)
println(tp53)
```

> Output
>
> ```groovy
> N E X T F L O W  ~  version 21.04.3
> Launching `closures.nf` [desperate_kalam] - revision: 9b078bd0b9
> [1, 2, 3, 4]
> [1, 4, 9, 16]
> [chr1, chr2, chr3, chr4]
> [chromosome:chr17, start:7661779, end:7687538, genome:GRCh38, gene:TP53, length:25759]
> ```

This demonstrates how to work with closures and their parameters in Nextflow, allowing for more flexible and dynamic code.

---

<details>
  <summary><b>CLICK HERE for another example</b></summary>

In Nextflow, the method `each()` when applied to a map can take a closure with two arguments, to which it passes the key-value pair for each entry in the map object.

Add the following to `closures.nf` and execute it using `nextflow run closures.nf`:

```groovy
// Closure with two parameters
printMap = { a, b -> println "$a with value $b" }
// Map object
my_map = [ chromosome : "chr17", start : 1, end : 83257441 ]
// Each iterates through each element
my_map.each(printMap)
```

Output

```groovy
N E X T F L O W  ~  version 21.04.3
Launching `closures.nf` [cheesy_khorana] - revision: 68b4c2dbbb
[1, 2, 3, 4]
[1, 4, 9, 16]
[chr1, chr2, chr3, chr4]
[chromosome:chr17, start:7661779, end:7687538, genome:GRCh38, gene:TP53, length:25759]
chromosome with value chr17
start with value 1
end with value 83257441
```

</details>

<br>

> Learn more about closures in the [Groovy documentation](http://groovy-lang.org/closures.html).

---

### Conditional Execution

#### If statement

- One of the most important features of any programming language is the ability to execute different code under different conditions. The simplest way to do this in Nextflow is to use the `if` construct. The `if` statement uses the same syntax common to other programming languages such as Java, C, JavaScript, etc.

> ```groovy
> if( < boolean expression > ) {
>    // true branch
> }
> else {
>    // false branch
> }
> ```
>
> > The else branch is optional. Also, curly brackets are optional when the branch defines just a single statement.

**Create a new `conditional.nf`; add the following and execute it using `nextflow run conditional.nf`:**

```groovy
x = 12
if( x > 10 )
    println "$x is greater than 10"

// Null, empty strings, and empty collections are evaluated to false in Groovy

list = [1, 2, 3]

if( list )
  println list
else
  println 'The list is empty'
```

> Output
>
> ```groovy
> N E X T F L O W  ~  version 21.04.3
> Launching `cond.nf` [amazing_jennings] - revision: 24d52ebd9f
> 12 is greater than 10
> [1, 2, 3]
> ```

---

## More resources

- **[Nextflow Scripting Documentation](https://www.nextflow.io/docs/latest/script.html#nextflow-scripting)**
- The complete Groovy language documentation is available at **[this link](http://groovy-lang.org/documentation.html#languagespecification)**.

> **Quick Recap**
>
> > - Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming language.
> > - To define a variable, simply assign a value to it e.g `a = 1`.
> > - Comments use the same syntax as in the C-family programming languages: `//` or multiline `/* */`.
> > - Multiple values can be stored in lists [value1, value2, value3, …] or maps [chromosome: 1, start :1].
> > - Lists are indexed and sliced with square brackets (e.g., list[0] and list[2..9])
> > - A closure is an expression (block of code) encased in {} e.g. `{it * it}`.

---

<h5><a href="/nextflow_varcal/docs/nextflow/nextflow_slurm" style="float: left"><b>Back to:</b>NF-Core @ HPC</a>

<a href="/nextflow_varcal/docs/nextflow/nextflow_channels" style="float: right"><b>Next:</b>NextFlow Channels</a></h5>
