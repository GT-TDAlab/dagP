#ifndef DGRAPH_H_
#define DGRAPH_H_

/* vertex ids start from 1 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

typedef int idxType; /* indices */
typedef long long ecType;  /* edge costs */
typedef int vwType;  /* vertex weights */

#define ecType_MAX  LLONG_MAX
#define vwType_MAX  INT_MAX
#define idxType_MAX INT_MAX

#define VWORKMULT 2 /* vertex work array's size */


#define DG_FRMT_UN   0  /* no weight ; no cost */
#define DG_FRMT_VW   1  /* vertex weight */
#define DG_FRMT_EC   2  /* edge cost */
#define DG_FRMT_VWEC 3  /* vertex weight edge cost */


typedef  struct{
    idxType nVrtx; /* num vertices */
    idxType nEdge; /* num edges */

    idxType *inStart;
    idxType *inEnd;
    idxType *in; /* for each vertex v, we will sore the edges
                  u->v in in[inStart[v]...inEnd[v]] */

    idxType *outStart;
    idxType *outEnd;
    idxType *out;/* for each vertex v, we will sore the edges
                  v-> in out[outStart[v]...outEnd[v]] */

    vwType *vw;
    ecType *ecIn; /* we store edge weights twice; this is at j for the edge (i->j)  */
    ecType *ecOut;/* we store edge weights twice; this is at i for the edge (i->j)  */
    int *hollow; /* == 1 if this node is from another part of the graph (not to be executed) */

    idxType *sources;
    idxType nbsources;

    idxType *targets;
    idxType nbtargets;

    /* information fields */
    int frmt;
    double totvw;
    double totec;
    idxType maxindegree;
    idxType maxoutdegree;
    vwType maxVW;
    ecType maxEC;

} dgraph;

void buildPartmatrix(dgraph *G, idxType * part, idxType nbpart, ecType** partmatrix);
int computePartmatrix(dgraph *G, idxType *part, ecType** partmatrix, idxType node, idxType movetopart);
int thereIsCycle(idxType nbpart, ecType** partmatrix);
int partSizeChecker(idxType* partsize, double ub_pw[], int nbpart, vwType maxvw);

void allocateDGraphData(dgraph *G, idxType nVrtx, idxType nEdge, int frmt);
void fillOutFromIn(dgraph *G);
void sortNeighborLists(dgraph *G);
void setVertexWeights(dgraph *G, vwType *vw);
void transitivelyReduce(dgraph *G);
void freeDGraphData(dgraph *G);
void dgraph_info(dgraph* G, int* maxindegree, int* minindegree, double* aveindegree, int* maxoutdegree, int* minoutdegree, double* aveoutdegree);
void set_dgraph_info(dgraph* G);

idxType sourcesList(dgraph *G, idxType* sources);
idxType sourcesListPart(dgraph* G, idxType* sources, idxType *part, int part_idx);
idxType outputList(dgraph *G, idxType* outputs);
int checkAcyclicity(dgraph *G, idxType *part, idxType nbpart);
void oneDegreeFirst(dgraph* G, idxType* order);
idxType farthestNode(dgraph* G, idxType startnode);
void computeToplevels(dgraph *G, idxType* toplevels);
void computeWeightedToplevels(dgraph* G, ecType* toplevels);
void computeBottomlevels(dgraph *G, idxType* bottomlevels);
void computeWeightedBottomlevels(dgraph* G, idxType* bottomlevels);
void computeToplevelsWithTopOrder(dgraph* G, idxType* toplevels, idxType* toporder);
void analyzeDGraph(dgraph *G);
void computeDistances(dgraph* G, idxType sourceNode, idxType* dist);
void connectedComponents(dgraph* G, idxType** components, idxType* sizes, idxType* nbcomp);


double computeLatency(dgraph* G, idxType* part, double l1, double l2);
void checkGraphAcyclicity(dgraph* G);
void reverseTab(idxType* tab, idxType size);

void printNodeInfo(dgraph *G,idxType* part,idxType node);

int printPartWeights(dgraph* G, idxType* part);

void components(dgraph *G, idxType *cmpts, idxType *nCmpts);
void transitiveReduction(dgraph* G);
void addSingleSourceTarget(dgraph* G, idxType* flag);
void checkDgraph(dgraph* G);
void copyDgraph(dgraph* Gcopy, dgraph* G);
ecType edgeCut(dgraph* G, idxType* part);
void reverseGraph(dgraph *G);
void randomizeWeights(dgraph *G, vwType vmin, vwType vmax, ecType emin, ecType emax);
void applyCCR(dgraph *G, double CCR);

idxType nbPart(dgraph* G, idxType* part, vwType* partsize);
ecType volume(dgraph* G, idxType* part, idxType nbPart);

#endif

