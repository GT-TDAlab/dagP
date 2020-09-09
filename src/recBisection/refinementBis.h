#ifndef REFINEMENTBIS_H
#define REFINEMENTBIS_H

#include "option.h"
// #include "vcycle.h"
#include "dgraph.h"
#include "utils.h"
#include "info.h"
// #include "dgraphPartition.h"
// #include "refinement.h"
#include "debug.h"
// #include "vcycle.h"
#include "dgraphTraversal.h"


typedef struct {
    ecType inFrom[2];
    ecType outTo[2];
} inFromOutToCnts;

void recBisRefinementStep(MLGP_option opt, int level, dgraph* graph, idxType* part, idxType nbpart, MLGP_info* info);
void refinementPostOrder_Max(dgraph *G, idxType *toporderpart,
                             double* lb_pw, double* ub_pw,
                             idxType *part, ecType* edgecut, MLGP_option opt, idxType *partsize);


#endif
