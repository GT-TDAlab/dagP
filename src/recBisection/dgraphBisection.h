#ifndef DGRAPHBISECTION_H
#define DGRAPHBISECTION_H

#include "dgraph.h"
#include "utils.h"
#include "initialBisection.h"
#include "refinementBis.h"

void greedyGraphGrowingTwoPhases(dgraph *G, double* ub_pw, idxType *part, MLGP_option opt, int fromPart, int otherPart, ecType *edgeCut);
void greedyGraphGrowingIn(dgraph *G, double* ub_pw, idxType *part, MLGP_option opt, int fromPart, int otherPart, ecType *edgeCut);
void greedyGraphGrowingEqualTopOrder(dgraph *G, double* ub_pw, idxType *part, MLGP_option opt, int fromPart, int otherPart, ecType *edgeCut);


#endif