## LICENSE

See LICENSE file

## CITATION AND CONTACT

For updates, latest version, please check: [http://tda.gatech.edu](http://tda.gatech.edu)


### Contact
    `tdalab@cc.gatech.edu`


Citation for the partitioner (BibTeX):

```
    @article { Herrmann19-SISC,
        author =  {Julien Herrmann and M. Yusuf {\"O}zkaya and Bora U{\c{c}}ar and Kamer Kaya and {\"{U}}mit V. {\c{C}}ataly{\"{u}}rek},
        title = {Multilevel Algorithms for Acyclic Partitioning of Directed Acyclic Graphs},
        journal = {{SIAM} Journal on Scientific Computing (SISC)},
        KEYWORDS = {directed graph ;  acyclic partitioning ;  multilevel partitioning},
        year = {2019},
        volume = {41},
        number = {4},
        pages = {A2117-A2145},
        doi = {10.1137/18M1176865},
        URL = {https://doi.org/10.1137/18M1176865},
    }
```

Citation for the Scheduler (BibTeX):

```
    @inproceedings { Ozkaya19-IPDPS,
        title = {A scalable clustering-based task scheduler for homogeneous processors using DAG partitioning},
        booktitle = {33rd IEEE International Parallel and Distributed Processing Symposium},
        author = {M. Yusuf \"Ozkaya and Anne Benoit and Bora U{\c{c}}ar and Julien Herrmann and {\"{U}}mit V. {\c{C}}ataly{\"{u}}rek},
        month = {May},
        year = {2019},
    }
```

## HOW TO BUILD

1.  install `scons`

2.  Prepare `config.py` file.

    if undirected initial partitioning will be used:

    - set `metis` and/or `scotch` to 1.

    - set up the paths for `metis`/`scotch` external libraries

    set up the compilers of choice

3.  Run `scons` to build

## HOW TO EXECUTE

1.  Build the project

2.  For the partitioner:
        run `./exe/rMLGP`

    if trying to run scheduler:
        run `./exe/rMLGPschedule`

    the applications will show the command line parameter usage

    The file name and number of parts are required. All other paramters are optional.

    Example usage:

    ./exe/rMLGP 2mm_10_20_30_40.dot 4 --print 1 --ratio 1.1


## GRAPH FORMAT

The application accepts a number of different formats:

  - `.dgraph` (`.dg`),
  - `.dot` (graphviz dot language), and
  - `.mtx` (Matrix Market File Format)

You can check out sample graphs or the websites for the respective graph formats.

The application has an option to generate binary versions of the input files.
This is, by default, enabled. If there is a binary version of the input,
the application prefers reading it instead of actual file.
Thus, if there is a change in the input file, where there is also a binary
version of a file with the same name, it will not be visible to the program.
The user should remove any `.bin` files for the inputs in this case.


## API USAGE

The dagP API provides five functions.

```
int dagP_init_parameters(MLGP_option *opt, const int nbPart);
int dagP_init_filename(MLGP_option* opt, char* file_name);
int dagP_opt_reallocUBLB(MLGP_option *opt, const int nbPart);
int dagP_read_graph(char* file_name, dgraph *G, const MLGP_option *opt);
ecType dagP_partition_from_dgraph(dgraph *G, const MLGP_option *opt, idxType* parts);
int dagP_free_graph(dgraph* G);
int dagP_free_option(MLGP_option *opt);s
```

First function initializes dagP options.

 - Options can be updated after initializing it using this command.

Second function initializes a file name for output files: (`file_name + "_dagP.out"`)

Third function can be used in case the number of parts to partition to is changed after parameter initialization.


Fourth function reads a directed graph from a file

 - The input graph format is inferred from the file extension.

 - If there is a binary version of the input file and the `use_binary_input` option `true`, then the partitioner will give priority to binary version of the input file.

Fifth function runs the partitioning algorithm on the given graph with the given options, then writes the number of partition assignment of the nodes to the `parts` array.

Last two functions free the respective variables.


An example code snipped using the API is available, named `useapi.cpp`.