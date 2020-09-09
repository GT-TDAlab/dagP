#ifndef RVCYCLE_H_
#define RVCYCLE_H_

#include <stdio.h>
#include <memory.h>
#include <string.h>

#include "vcycle2way.h"


typedef struct rcoarsen {
    /*Coarsening process*/
    coarsen  *coars; /*Current coarsen*/
    struct rcoarsen *previous_rcoarsen; /*Previous rcoarsen in the recursion*/
    struct rcoarsen *next_rcoarsen1; /*One of the next rcoarsen in the recursion*/
    struct rcoarsen *next_rcoarsen2; /*One of the next rcoarsen in the recursion*/

    /*Link between the recursion*/
    idxType* previous_index; /*previous_index[i] is the index of i in previous_rcoarsen*/
    idxType* next_index; /*next_index[i] is the new index of node i in next_rcoarsen1 or next_rcoarsen2*/

} rcoarsen;

#include "info.h"

rcoarsen* initializeRCoarsen(coarsen* icoars);
void freeRCoarsenHeader(rcoarsen *rcoars);
rcoarsen* rVCycle(dgraph *igraph, const MLGP_option opt, rMLGP_info* info);

#endif
