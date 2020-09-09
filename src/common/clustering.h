#ifndef MATCHING_H_
#define MATCHING_H_

#include "dgraph.h"
#include "option.h"
#include "dgraphTraversal.h"
#include "vcycle2way.h"


void computeClustering(const MLGP_option opt, int co_level, coarsen *C, idxType* part);

#endif
