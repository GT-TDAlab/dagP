#ifndef VCYCLE2WAY_H
#define VCYCLE2WAY_H

#include <stdio.h>
#include <memory.h>
#include <string.h>

#include "dgraph.h"
#include "info.h"
#include "option.h"

typedef struct coarsen {
    /*Coarsening process*/
    dgraph *graph; /*Current graph*/
    struct coarsen *previous_coarsen; /*Previous graph in the coarsening*/
    struct coarsen *next_coarsen; /*Graph obtained by applying the matching below*/

    /*Matching description*/
    idxType  *leader; /*If leader[i] == leader[j], nodes i and j are in the same node in next_graph (next coarsen)*/
    idxType  *new_index; /*new_index[i] is the new index of node i in next_coarsen*/

    /*Partition description*/
    idxType  *part; /*part[i] is the index of the partition set node i belongs to (starting from 0)*/
    dgraph   *quotient; /*quotient graph for the given partition*/
    idxType  *flag; /*nodes that doesn't have the same flag can't be matched in the coarsening phase*/
    int      nbpart;
} coarsen;



coarsen* initializeCoarsen(dgraph *graph);
void freeCoarsen(coarsen *coars);
void freeCoarsenHeader(coarsen *coars);
dgraph *initializeNextGraph(coarsen *coars, MLGP_option opt);
coarsen *projectBack(coarsen* coars);

coarsen *VCycle2way(dgraph* graph, idxType* flag, const MLGP_option opt, MLGP_info* info);


#endif
