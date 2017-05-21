Fast Core-Genome Identity (ANI) Estimation
========================================================================

### Download and compile

Follow [`INSTALL.txt`](INSTALL.txt) to compile the code

### Usage

* **One to One.** Single query and reference genome:

```sh
mashcgi -q [QUERY_GENOME] -s [REFERENCE_GENOME] -o [OUTPUT_FILE] 
```

Here QUERY\_GENOME and REFERENCE_GENOME are the query genome assemblies (in fasta format) for which ANI is required.

* **One to Many.** Single query genome and multiple reference genomes:

```sh
mashcgi -q [QUERY_GENOME] -sl [REFERENCE_LIST] -o [OUTPUT_FILE]
```

REFERENCE\_LIST is supposed to be a file containing paths to reference genomes, 1 per line.

* **Many to Many.** Multiple query genomes and multiple reference genomes:

```sh
mashcgi --ql [QUERY_LIST] -sl [REFERENCE_LIST] -o [OUTPUT_FILE]
```
Again, QUERY\_LIST and REFERENCE\_LIST are files containing paths to genomes, 1 per line.

In all above use cases, OUTPUT\_FILE will contain space delimited row(s) with query genome, reference genome, ANI value, count of bidirectional fragment mappings, and total query fragments. Additional log is printed to stderr. We suggest users to do minimal quality check of the input genome assemblies (both reference and query), especially the N50 be ≥10 Kbp.
