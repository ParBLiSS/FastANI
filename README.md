Fast Core-Genome Identity (ANI) Estimation
========================================================================

### Download and compile

Follow [`INSTALL.txt`](INSTALL.txt) to compile the code

### Run

Inside the build directory, 

```sh
mashcgi --sl [REFERENCE_LIST] -q [QUERY_GENOME] -o [OUTPUT_FILE]  --pi 80
```

Here, REFERENCE_LIST is a file containing paths to reference genomes, 1 per line. QUERY_GENOME is the query genome for which ANI is computed against the reference. QUERY_GENOME should be fragmented into pieces of non-overlapping 5K sequences before execution. 

OUTPUT_FILE will contain list of reference genomes and the corresponding ANI values. Third column will denote the count of fragments mapped onto the reference. ANI values will >=50 mapped fragments should be trusted.
