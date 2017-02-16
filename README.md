Fast Core-Genome Identity (ANI) Estimation
========================================================================

### Download and compile

Follow [`INSTALL.txt`](INSTALL.txt) to compile the code

### Usage

* Multiple reference genomes and single query genome:

```sh
mashcgi -sl [REFERENCE_LIST] -q [QUERY_GENOME] -o [OUTPUT_FILE] --pi 80
```

Here, REFERENCE_LIST is a file containing paths to reference genomes, 1 per line. QUERY_GENOME is the query genome for which ANI is computed against the reference.

* Single reference and query genome:

```sh
mashcgi -s [REFERENCE_GENOME] -q [QUERY_GENOME] -o [OUTPUT_FILE] --pi 80
```

In both types of use cases, OUTPUT_FILE will be a space delimited file with reference genome, ANI value, count of bidirectional fragment mappings, reference genome N50 and reference genome total length.
