Fast Core-Genome Identity (ANI) Estimation
========================================================================

### Download and compile

Follow [`INSTALL.txt`](INSTALL.txt) to compile the code

### Run

Inside the build directory, 

```sh
mashcgi --sl [REFERENCE_LIST] -q [QUERY_GENOME] -o [OUTPUT_FILE] -m [FRAGMENT_LENGTH] --pi 80
```

Here, REFERENCE_LIST is a file containing paths to reference genomes, 1 per line. QUERY_GENOME is the query genome for which ANI is computed against the reference. FRAGMENT_LENGTH is the size of each fragment cut from query genome, current default is 5,000 bp. Each of these fragments is mapped onto the reference.

OUTPUT_FILE will contain list of reference genomes and the corresponding ANI values. Third column will denote the count of query fragments mapped onto the reference. ANI values with >=50 mapped fragments should be trusted.
