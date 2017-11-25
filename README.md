FastANI
========================================================================
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)

FastANI is a fast alignment-free estimator of Average Nucleotide Identity (ANI) between two genomes. ANI is defined as mean nucleotide identity of orthologous gene pairs between two microbial genomes. This software supports comparison of both complete and draft genomes. FastANI follows a similar workflow as described by [Goris et al. 2007](http://www.ncbi.nlm.nih.gov/pubmed/17220447). However, it avoids expensive sequence alignments and uses [Mashmap](https://github.com/marbl/MashMap) as its sequence mapping engine. Based on our experiments with complete and draft genomes, FastANI's accuracy is on par with [BLAST-based ANI solver](http://enve-omics.ce.gatech.edu/ani/) and yields two to three orders of magnitude speedup. Therefore, it is useful for ANI analysis of large number of genome pairs. Detailed results and comparisons with existing methods are described in the paper. 

### Download and Compile

Clone the software from Github and follow [`INSTALL.txt`](INSTALL.txt) to compile the code. 

### Usage Summary

* **Produce help page.** Quickly check the software usage and available command line options.

```sh
$ ./fastANI -h
```

* **One to One.** Compute ANI between single query and single reference genome:

```sh
$ ./fastANI -q [QUERY_GENOME] -s [REFERENCE_GENOME] -o [OUTPUT_FILE] 
```

Here QUERY\_GENOME and REFERENCE\_GENOME are the query genome assemblies in fasta or multi-fasta format.

* **One to Many.** Compute ANI between single query genome and multiple reference genomes:

```sh
$ ./fastANI -q [QUERY_GENOME] --sl [REFERENCE_LIST] -o [OUTPUT_FILE]
```

For above use case, REFERENCE\_LIST should be a file containing directory paths to reference genomes, one per line.

* **Many to Many.** When there are multiple query genomes and multiple reference genomes:

```sh
$ ./fastANI --ql [QUERY_LIST] --sl [REFERENCE_LIST] -o [OUTPUT_FILE]
```
Again, QUERY\_LIST and REFERENCE\_LIST are files containing paths to genomes, one per line.

In all above use cases, OUTPUT\_FILE will contain space delimited row(s) with query genome, reference genome, ANI value, count of bidirectional fragment mappings, and total query fragments. Alignment fraction (wrt. the query genome) is simply the ratio of mappings and total fragments. Two genome assemblies are provided in [data](data) folder to do a quick test run. 

We suggest users to do minimal quality check of their input genome assemblies (both reference and query), especially the N50 be ≥10 Kbp.

### An Example Run

* **One to One.** Compare Escherichia coli and Shigella flexneri genomes provided in the [data](data) folder.

```sh
$ ./fastANI -q data/Shigella_flexneri_2a_01.fna -s data/Escherichia_coli_str_K12_MG1655.fna -o fastani.out 
```

Expect output log in the following format in the console:

```sh
$ ./fastANI -q data/Shigella_flexneri_2a_01.fna -s data/Escherichia_coli_str_K12_MG1655.fna -o fastani.out 
>>>>>>>>>>>>>>>>>>
Reference = [data/Escherichia_coli_str_K12_MG1655.fna]
Query = [data/Shigella_flexneri_2a_01.fna]
Kmer size = 16
Fragment length = 3000
ANI output file = fastani.out
>>>>>>>>>>>>>>>>>>
....
INFO, skch::main, Time spent post mapping : 0.00310319 sec
```

The output is saved in file fastani.out . It should contain the ANI estimate between E. coli and S. flexneri genomes.

```sh
$ cat fastani.out
data/Shigella_flexneri_2a_01.fna data/Escherichia_coli_str_K12_MG1655.fna 97.7443 1305 1608
```

### Visualize Conserved Regions b/w Two Genomes

FastANI supports visualization of the reciprocal mappings computed between two genomes. 
Getting this visualization requires a one to one comparison using FastANI as discussed above, except an additional flag --visualize should be provided. 
This flag forces FastANI to output mapping file (with `.visual` extension) that contains information of all the reciprocal mappings. Finally, an [R script](scripts) is provided in the repository which uses [genoPlotR](https://cran.r-project.org/web/packages/genoPlotR/index.html) package to plot these mappings. Here we show an example run using two genomes: Bartonella quintana ([GenBank: CP003784.1](https://www.ncbi.nlm.nih.gov/nuccore/CP003784.1)) and Bartonella henselae ([NCBI Reference Sequence: NC_005956.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_005956.1)).

```sh
$ ./fastANI -q B_quintana.fna -s B_henselae.fna --visualize -o fastani.out
$ Rscript scripts/visualize.R B_quintana.fna B_henselae.fna fastani.out.visual
```

Using above commands, we get a plot file fastani.out.visual.pdf displayed below. Each red line segment denotes a reciprocal mapping between two genomes, indicating their evolutionary conserved regions.

![FastANI Visualization](https://alurulab.cc.gatech.edu/sites/all/images/fastani/fastani.out.visual.jpg)

### Parallelization

As of now, FastANI doesn't support parallelization internally. However, for one-to-many or many-to-many genome comparisons, users can simply divide their reference database into multiple chunks, and execute them as parallel processes. We provide a [helper script](scripts) to do this splitting.

### Troubleshooting

Users are welcome to report any issue or feedback related to FastANI by posting a [Github issue](https://github.com/ParBLiSS/FastANI/issues).

### Publication

To be updated soon...
