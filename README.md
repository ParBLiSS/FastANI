Fast Core-Genome Identity (ANI) Estimation
========================================================================

### Download and compile

Follow [`INSTALL.txt`](INSTALL.txt) to compile the code

### Usage

* Multiple reference genomes and single query genome:

```sh
mashcgi -sl [REFERENCE_LIST] -q [QUERY_GENOME] -o [OUTPUT_FILE]
```

Here, REFERENCE\_LIST is a file containing paths to reference genomes, 1 per line. QUERY\_GENOME is the query genome for which ANI is computed against the reference.

* Multiple reference genomes and multiple query genomes:

```sh
mashcgi -sl [REFERENCE_LIST] --ql [QUERY_LIST] -o [OUTPUT_FILE]
```
Again, REFERENCE\_LIST and QUERY\_LIST are files containing paths to genomes, 1 per line.

* Single reference and query genome:

```sh
mashcgi -s [REFERENCE_GENOME] -q [QUERY_GENOME] -o [OUTPUT_FILE] 
```

In both types of use cases, OUTPUT\_FILE will contain space delimited row(s) with reference genome, ANI value, count of bidirectional fragment mappings, total query fragments. Additional log is printed to stderr. We strongly suggest users to do minimal quality check of the input genome assemblies (both reference and query), especially the N50 be ≥10 Kbp.
