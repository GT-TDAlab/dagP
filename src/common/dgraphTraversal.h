#ifndef DGRAPHTRAVERSAL_H
#define DGRAPHTRAVERSAL_H

#include "dgraph.h"
#include "utils.h"

void natural_with_part(dgraph *G, idxType* toporder, idxType *part, idxType *toporderpart, idxType nbpart);
void topologic_natural_with_part(dgraph *G, idxType* toporder, idxType *part, idxType *toporderpart, idxType nbpart);

void DFStopsort_with_part(dgraph *G, idxType* toporder, idxType *part, idxType *toporderpart, idxType nbpart);
void BFStopsort_with_part(dgraph *G, idxType* toporder, idxType *part, idxType *toporderpart, idxType nbpart);
void randDFStopsort_with_part(dgraph* G, idxType *part, int nbpart, idxType *toporder);

void topSortOnParts(dgraph *G, idxType *part, idxType *toporder, idxType nbpart);
void topSortOnParts_BFS(dgraph *G, idxType *part, idxType *toporder, idxType nbpart);
void randTopSortOnParts(dgraph *G, idxType *part, idxType *toporder, idxType nbpart);
void topsortPriorities(dgraph *G, idxType *toporder, ecType *priorities);
void checkTopsort(dgraph *G, idxType *toporder);

void DFStopsort(dgraph *G, idxType *toporder);
void BFStopsort(dgraph *G, idxType *toporder);
void randDFStopsort(dgraph *G, idxType *toporder);
void randBFStopsort(dgraph *G, idxType *toporder);
void DFSsort(dgraph *G, idxType *toporder);
void BFSsort(dgraph *G, idxType *toporder);
void randDFSsort(dgraph *G, idxType *toporder);
void randBFSsort(dgraph *G, idxType *toporder);
void randTopsort(dgraph *G, idxType *toporder);
void mixtopsort(dgraph *G, idxType *toporder,int priority,int first);
void randTopsort(dgraph *G, idxType *toporder);


#endif
