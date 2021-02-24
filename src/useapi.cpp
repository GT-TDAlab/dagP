#include <cstdio>
#include "dagP.h"

int main(int argc, char *argv[])
{
    int nbParts = 2;
    // get input file name (and optionally the number of parts) as a command line argument
    if (argc < 2) {
        printf("Program requires the input file as command line argument!\n");
        printf ("Usage: %s <input-file> [<number-of-parts>]", argv[0]);
        return 0;
    }
    else if (argc > 3) {
        printf("Too many arguments given!\n");
        printf ("Usage: %s <input-file> [<number-of-parts>]", argv[0]);
        return 0;
    }
    if (argc == 3) {
        nbParts = atoi(argv[2]);
    }

    MLGP_option opt;
    dgraph G;

    // initialize the ``opt'' variable with number of parts goal = 2
    dagP_init_parameters (&opt, 2); // initializes default parameters

    dagP_init_filename(&opt, argv[1]); // initialize the input file name and then the  output file name

    printf ("input: %s\nnbPart: %d\n", opt.file_name, opt.nbPart);

    // You can reallocate the related variables for a different number of parts later on,
    // instead of providing it on the init_parameters function:
    // reallocate the part upper bound and lower bound arrays for a different number of parts
    // ub and lb arrays by default filled with -1. To set arbitrary values, simply run a for loop
    dagP_opt_reallocUBLB (&opt, nbParts);

    // opt.print = 10; // parameters can then be updated by direct access
    // opt.co_match = 1;
    // opt.conpar = 1;
    // opt.inipart = 6; // opt.inipart = IP_UNDIR;
    opt.runs = 5; // number of runs before selecting best edge cut

    // read the graph from file with name from command-line arguments
    dagP_read_graph (argv[1], &G, &opt);

    // allocate `parts` array with number of nodes limit = G.nVrtx + 1
    idxType *parts = (idxType*) calloc((G.nVrtx+1), sizeof(idxType));
    if (parts == NULL)
        printf("Could not allocate `parts` array.\n");

    // to the partitioning return value is edge cut, part assignments are written over parts array
    ecType x = dagP_partition_from_dgraph(&G, &opt, parts);

    // // node id's are 1-indexed
    // for(idxType i=1; i<= G.nVrtx; ++i){
    //     printf("part[node:%d] = %d\n", i, parts[i]);
    // }

    printf ("edge cut: %d\n", (int) x);

    free(parts);
    dagP_free_option(&opt); // frees the arrays in opt
    dagP_free_graph(&G); // frees the arrays in the graph and resets the variables
    return 0;
}


/*
Compile and run with:
```
scons
g++ -Wall -o useapi src/useapi.cpp -I src/recBisection/ -I src/common lib/libdagp.a -L/Users/myozka/Documents/metis-5.1.0/build/Darwin-x86_64/libmetis -lmetis -lm
./useapi 2mm_10_20_30_40.dot 2
```
*/