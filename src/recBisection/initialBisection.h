#ifndef INITIALBISECTION_H_
#define INITIALBISECTION_H_

#include "option.h"
#include "vcycle2way.h"

/* returns #parts found */
void fixUndirBisection(ecType* edgecut, coarsen* coars, idxType* part, MLGP_option opt, MLGP_info* info);
void undirBisectionLibrary(ecType* edgecut, coarsen* coars, idxType* part, MLGP_option opt, int seed, MLGP_info* info);

void initialBisection(MLGP_option opt, coarsen*  coars, MLGP_info* info);
#endif
