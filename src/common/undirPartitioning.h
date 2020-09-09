#ifndef UNDIRPARTIONING_H_
#define UNDIRPARTIONING_H_

#ifdef dagP_METIS
#include "metis.h"

void dgraph_to_metis(dgraph* G, idx_t* in, idx_t* adj, idx_t* vw, idx_t* adjw);
#else
typedef int idx_t;
#endif

#include "option.h"
#include "rvcycle.h"

void conpar(dgraph* G, idxType* constraint, const MLGP_option opt, MLGP_info* info);
void undirPartitioningLibrary(dgraph* G, idxType* part, const MLGP_option opt, int seed);
void undirBisectionRandom(dgraph* graph, idxType* part);
void undirBisectionGGG(dgraph* G, idxType* part, MLGP_option opt);
#endif
