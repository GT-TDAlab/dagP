#include <stdio.h>
#include <stdlib.h>

#include <limits.h>
#include <math.h>
#include <float.h>

#include "clustering.h"
#include "utils.h"


void computeMatchNOrder(dgraph *G, const MLGP_option opt, idxType *norder)
{
    switch (opt.co_norder) {
        case CO_NORD_RANDDFSTOP:
            randDFStopsort(G, norder);
            break;
        case CO_NORD_RANDBFSTOP:
            randBFStopsort(G, norder);
            break;
        case CO_NORD_RANDTOP:
            randTopsort(G, norder);
            break;
        case CO_NORD_DFSTOP:
            DFStopsort(G, norder);
            break;
        case CO_NORD_BFSTOP:
            BFStopsort(G, norder);
            break;
        case CO_NORD_RANDDFS:
            randDFSsort(G, norder);
            break;
        case CO_NORD_RANDBFS:
            randBFSsort(G, norder);
            break;
        case CO_NORD_DFS:
            DFSsort(G, norder);
            break;
        case CO_NORD_BFS:
            BFSsort(G, norder);
            break;
        default:
            // if nothing just do Random
            // case CO_NORD_RAND:
            shuffleTab(1, G->nVrtx, norder+1);
            break;
    }

    if (opt.co_norder_reverse == 1)
        reverseTab(norder, G->nVrtx);

}



int lookNeighbors(coarsen* C, idxType i, idxType* neighidx, int* isout, idxType* toplevels, int* forbidenu, int* forbidenv, int* matched, const MLGP_option opt) {
    /* This function looks at the neighbors of i in G according to the options
     * It modifies the list of idexes neighidx and keep in memory if the neighbor is an in-edge or an out-edge with isout table
     * It returns the number of valid neighbors
     */
    dgraph* G = C->graph;
    double weight_lim = (0.1 * G->totvw) / opt.nbPart;
    idxType k;
    int nbneigh = 0;
    if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT)) {
        int only_input_i = (G->inEnd[i] - G->inStart[i] + 1 == 1);
        for (k = G->inStart[i]; k <= G->inEnd[i]; k++) {
            if ((opt.co_eorder == CO_EORD_NONE) && (nbneigh > 0))
                break;
            idxType node = G->in[k];
            if (node > G->nVrtx)
                u_errexit("In matching node = %d for i = %d, j=%d\n", (int) node, i, k);
            if (matched[node] == 1)
                continue;
            if (forbidenu[node] == 1)
                continue;
            int only_output_node = (G->outEnd[node] - G->outStart[node] + 1 == 1);
            if ((abs(toplevels[node] - toplevels[i]) > 1) && (!only_input_i) && (!only_output_node))
                continue;
            if (G->vw[i] + G->vw[node] > weight_lim)
                continue;
            if (C->flag[i] != C->flag[node])
                continue;
            neighidx[nbneigh] = k;
            isout[nbneigh] = 0;
            nbneigh++;
        }
    }
    if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT)) {
        int only_output_i = (G->outEnd[i] - G->outStart[i] + 1) == 1;
        for (k = G->outStart[i]; k <= G->outEnd[i]; k++) {
            if ((opt.co_eorder == CO_EORD_NONE) && (nbneigh > 0))
                break;
            idxType node = G->out[k];
#if _dGP_DEBUG >= PD_HIGH
            if (node > G->nVrtx)
                u_errexit("In matching node = %d for i = %d, k=%d\n", (int) node, i, k);
#endif
            if (matched[node] == 1)
                continue;
            if (forbidenv[node] == 1)
                continue;
            int only_input_node = (G->inEnd[node] - G->inStart[node] + 1) == 1;
            if ((abs(toplevels[i] - toplevels[node]) > 1) && (!only_input_node) && (!only_output_i))
                continue;
            if (G->vw[i] + G->vw[node] > weight_lim)
                continue;
            if (C->flag[i] != C->flag[node])
                continue;
            neighidx[nbneigh] = k;
            isout[nbneigh] = 1;
            nbneigh++;
        }
    }
    return nbneigh;
}

void selectBestNeighbors(dgraph* G, idxType i, idxType* innode, idxType* outnode, idxType* neighidx, idxType nbneigh, int* isout, const MLGP_option opt) {
    /* This function select the best neighbor of i in the neighidx table
     * It modifies innode and outnode to the correspond values for the best neighbor
     */
    idxType k;
    if (isout[0] == 1) {
        *innode = i;
        *outnode = G->out[neighidx[0]];
    } else { // if (isout[0] == 0) {
        *innode = G->in[neighidx[0]];
        *outnode = i;
    }
    double best_cost = 0.0, newcost;
    switch (opt.co_eorder) {
        case CO_EORD_WEIGHT :
            for (k = 0; k < nbneigh; k++) {
                if (isout[k] == 1) {
                    newcost = G->ecOut[neighidx[k]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *innode = i;
                        *outnode = G->out[neighidx[k]];
                    }
                }
                if (isout[k] == 0) {
                    newcost = G->ecIn[neighidx[k]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *innode = G->in[neighidx[k]];
                        *outnode = i;
                    }
                }
            }
            break;

        case CO_EORD_RATIO :
            for (k = 0; k < nbneigh; k++) {
                if (isout[k] == 1) {
                    newcost = (double) G->ecOut[neighidx[k]] / (double) G->vw[G->out[neighidx[k]]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *innode = i;
                        *outnode = G->out[neighidx[k]];
                    }
                }
                if (isout[k] == 0) {
                    newcost = (double) G->ecIn[neighidx[k]] / (double) G->vw[G->in[neighidx[k]]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *innode = G->in[neighidx[k]];
                        *outnode = i;
                    }
                }
            }
            break;
    }
}

void matchingCoTop(coarsen* C, idxType *toplevels, idxType *norder, idxType *leader, const MLGP_option opt)
{
    dgraph* G = C->graph;
    idxType i, n, k;
    int* forbidenu, *forbidenv, *isout, *matched;
    idxType* neighidx;
    neighidx = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "neighidx");
    isout = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "isout");
    matched = (int*) calloc(G->nVrtx+1, sizeof(int));
    if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT))
        forbidenu = (int *) calloc(G->nVrtx + 1, sizeof(int));
    if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT))
        forbidenv = (int *) calloc(G->nVrtx + 1, sizeof(int));

    for (i=1; i<=G->nVrtx; i++) {
        leader[i] = i;
    }
    for (n=1; n<=G->nVrtx; n++) {
		i = norder[n];
		if (matched[i] == 1)
	    	continue;
        int nbneigh = lookNeighbors(C, i, neighidx, isout, toplevels, forbidenu, forbidenv, matched, opt);

        if (nbneigh == 0)
            continue;
        idxType innode, outnode;
        selectBestNeighbors(G, i, &innode, &outnode, neighidx, nbneigh, isout, opt);

        leader[outnode] = innode;
        matched[innode] = 1;
        matched[outnode] = 1;
        if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT))
            for (k=G->inStart[outnode]; k<= G->inEnd[outnode]; k++)
                if (abs(toplevels[outnode] - toplevels[G->in[k]]) <= 1)
                    forbidenu[G->in[k]] = 1;
        if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT))
            for (k=G->outStart[innode]; k<= G->outEnd[innode]; k++)
                if (abs(toplevels[innode] - toplevels[G->out[k]]) <= 1)
                    forbidenv[G->out[k]] = 1;
    }
    free(neighidx);
    free(isout);
    free(matched);
    if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT))
        free(forbidenu);
    if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT))
        free(forbidenv);
}

int lookNeighborsAggForbidden(coarsen* C, idxType i, idxType* neighidx, int* isout, idxType* toplevels, idxType* leader, int* matchedu, int* matchedv, idxType* nbbadinput, idxType* badinputnode, idxType* nbbadoutput, idxType* badoutputnode, vwType* weight, const MLGP_option opt) {
    /* This function looks at the neighbors of i in G according to the options for agglomerative coarsening
     * It modifies the list of idexes neighidx and keep in memory if the neighbor is an in-edge or an out-edge with isout table
     * It returns the number of valid neighbors
     */
    dgraph* G = C->graph;
    int nbneighbors = 0;
    double weight_lim = (0.1 * G->totvw) / opt.nbPart;
    idxType j,k;
    if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT)) {
        if (nbbadinput[i] == 1) {
            idxType node = badinputnode[i];
            if ((weight[leader[node]] + G->vw[i] <= weight_lim) && (abs(toplevels[i] - toplevels[node]) <= 1) &&
                (matchedv[node] == 0) && ((matchedu[node] == 1) || (nbbadoutput[node] == 0))) {
                for (k = G->inStart[i]; k <= G->inEnd[i]; k++)
                    if (G->in[k] == badinputnode[i]) {
                        neighidx[nbneighbors] = k;
                        break;
                    }
                isout[nbneighbors] = 0;
                if (C->flag[i] == C->flag[G->in[k]])
                    nbneighbors++;
            }
        }
        if (nbbadinput[i] == 0) {
            for (j = G->inStart[i]; j <= G->inEnd[i]; j++) {
                if ((opt.co_eorder == CO_EORD_NONE) && (nbneighbors > 0))
                    break;
                idxType node = G->in[j];
                if ((matchedv[node] == 1) ||
                    (abs(toplevels[i] - toplevels[node]) > 1) ||
                    ((matchedu[node] == 0) && (nbbadoutput[node] > 0)) ||
                    (weight[leader[G->in[j]]] + G->vw[i] > weight_lim))
                    continue;
                if (C->flag[i] != C->flag[node])
                    continue;
                neighidx[nbneighbors] = j;
                isout[nbneighbors] = 0;
                nbneighbors++;
            }
        }
    }
    if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT)) {
        if (nbbadoutput[i] == 1) {
            idxType node = badoutputnode[i];
            if ((weight[leader[node]] + G->vw[i] <= weight_lim) && (abs(toplevels[node] - toplevels[i]) <= 1) &&
                (matchedu[node] == 0) && ((matchedv[node] == 1) || (nbbadinput[node] == 0))) {
                for (k = G->outStart[i]; k <= G->outEnd[i]; k++)
                    if (G->out[k] == node) {
                        neighidx[nbneighbors] = k;
                        break;
                    }
                isout[nbneighbors] = 1;
                if (C->flag[i] == C->flag[G->out[k]])
                    nbneighbors++;

            }
        }
        if (nbbadoutput[i] == 0) {
            for (j = G->outStart[i]; j <= G->outEnd[i]; j++) {
                if ((opt.co_eorder == CO_EORD_NONE) && (nbneighbors > 0))
                    break;
                idxType node = G->out[j];
                if ((matchedu[node] == 1) ||
                    (abs(toplevels[node] - toplevels[i]) > 1) ||
                    ((matchedv[node] == 0) && (nbbadinput[node] > 0)) ||
                    (weight[leader[node]] + G->vw[i] > weight_lim))
                    continue;
                if (C->flag[i] != C->flag[node])
                    continue;
                neighidx[nbneighbors] = j;
                isout[nbneighbors] = 1;
                nbneighbors++;
            }
        }
    }
    return nbneighbors;
}

int lookNeighborsAggCheck(coarsen* C, idxType i, idxType* neighidx, int* isout, idxType* toplevels, idxType* leader, int* matchedu, int* matchedv, vwType* weight, const MLGP_option opt) {
    /* This function looks at the neighbors of i in G according to the options for agglomerative coarsening check
     * It modifies the list of idexes neighidx and keep in memory if the neighbor is an in-edge or an out-edge with isout table
     * It returns the number of valid neighbors
     */
    dgraph* G = C->graph;
    int nbneighbors = 0;
    double weight_lim = (0.1 * G->totvw) / opt.nbPart;
    idxType j,k;
    if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT)) {
        for (j = G->inStart[i]; j <= G->inEnd[i]; j++) {
            if ((opt.co_eorder == CO_EORD_NONE) && (nbneighbors > 0))
                break;
            idxType node = G->in[j];
            if ((matchedv[node] == 1) ||
                (abs(toplevels[i] - toplevels[node]) > 1) ||
                (weight[leader[G->in[j]]] + G->vw[i] > weight_lim))
                continue;
            if (C->flag[i] != C->flag[node])
                continue;
            neighidx[nbneighbors] = j;
            isout[nbneighbors] = 0;
            nbneighbors++;
        }
    }
    if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT)) {
        for (j = G->outStart[i]; j <= G->outEnd[i]; j++) {
            if ((opt.co_eorder == CO_EORD_NONE) && (nbneighbors > 0))
                break;
            idxType node = G->out[j];
            if ((matchedu[node] == 1) ||
                (abs(toplevels[node] - toplevels[i]) > 1) ||
                (weight[leader[node]] + G->vw[i] > weight_lim))
                continue;
            if (C->flag[i] != C->flag[node])
                continue;
            neighidx[nbneighbors] = j;
            isout[nbneighbors] = 1;
            nbneighbors++;
        }
    }
    return nbneighbors;
}

void selectBestNeighborsAgg(dgraph* G, idxType* node, int* nodeisout, idxType* neighidx, idxType nbneighbors, int* isout, idxType* leader, vwType* weight, ecType* edgesweight, const MLGP_option opt) {
    /* This function select the best neighbor of i in the neighidx table for agglomerative coarsening
     * It modifies node and nodeisout to the corresponding values for the best neighbor
     */
    idxType k;
    if (isout[0] == 1)
        *node = G->out[neighidx[0]];
    if (isout[0] == 0)
        *node = G->in[neighidx[0]];
    double best_cost = 0.0, newcost;
    *nodeisout = isout[0];
    switch (opt.co_eorder) {
        case CO_EORD_WEIGHT :
            for (k = 0; k < nbneighbors; k++) {
                if (isout[k] == 1) {
                    newcost = (double) G->ecOut[neighidx[k]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *node = G->out[neighidx[k]];
                        *nodeisout = isout[k];
                    }
                } else {/*if (isout[k] == 0){*/
                    newcost = (double) G->ecIn[neighidx[k]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *node = G->in[neighidx[k]];
                        *nodeisout = isout[k];
                    }
                }
            }
            break;

        case CO_EORD_RATIO :
            for (k = 0; k < nbneighbors; k++) {
                if (isout[k] == 1) {
                    newcost = (double) G->ecOut[neighidx[k]] / (double) G->vw[G->out[neighidx[k]]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *node = G->out[neighidx[k]];
                        *nodeisout = isout[k];
                    }
                } else {/*if (isout[k] == 0){*/
                    newcost = (double) G->ecIn[neighidx[k]] / (double) G->vw[G->in[neighidx[k]]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *node = G->in[neighidx[k]];
                        *nodeisout = isout[k];
                    }
                }
            }
            break;

        case CO_EORD_AGGWEIGHT :
            for (k = 0; k < nbneighbors; k++) {
                if (isout[k] == 1)
                    edgesweight[leader[G->out[neighidx[k]]]] += (double) G->ecOut[neighidx[k]];
                else
                    edgesweight[leader[G->in[neighidx[k]]]] += (double) G->ecIn[neighidx[k]];
            }
            for (k = 0; k < nbneighbors; k++) {
                if (isout[k] == 1) {
                    newcost = edgesweight[leader[G->out[neighidx[k]]]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *node = G->out[neighidx[k]];
                        *nodeisout = isout[k];
                    }
                    edgesweight[leader[G->out[neighidx[k]]]] = 0;
                } else {/*if (isout[k] == 0){*/
                    newcost = edgesweight[leader[G->in[neighidx[k]]]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *node = G->in[neighidx[k]];
                        *nodeisout = isout[k];
                    }
                    edgesweight[leader[G->in[neighidx[k]]]] = 0;
                }
            }
            break;

        case CO_EORD_AGGRATIO :
            for (k = 0; k < nbneighbors; k++) {
                if (isout[k] == 1)
                    edgesweight[leader[G->out[neighidx[k]]]] += (double) G->ecOut[neighidx[k]];
                else
                    edgesweight[leader[G->in[neighidx[k]]]] += (double) G->ecIn[neighidx[k]];
            }
            for (k = 0; k < nbneighbors; k++) {
                if (isout[k] == 1) {
                    newcost = edgesweight[leader[G->out[neighidx[k]]]] / weight[leader[G->out[neighidx[k]]]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *node = G->out[neighidx[k]];
                        *nodeisout = isout[k];
                    }
                    edgesweight[leader[G->out[neighidx[k]]]] = 0;
                } else {/*if (isout[k] == 0){*/
                    newcost = edgesweight[leader[G->in[neighidx[k]]]] / weight[leader[G->in[neighidx[k]]]];
                    if (newcost > best_cost) {
                        best_cost = newcost;
                        *node = G->in[neighidx[k]];
                        *nodeisout = isout[k];
                    }
                    edgesweight[leader[G->in[neighidx[k]]]] = 0;
                }
            }
            break;
    }
}


void clusteringCoTop(coarsen* C, idxType *toplevels, idxType *norder, idxType* leader, const MLGP_option opt)
{
    dgraph* G = C->graph;
    idxType i, n, k;

    double t1 = -u_wseconds();

    idxType *neighidx = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "neighidx");
    int *isout = (idxType*) umalloc(sizeof(int)*(G->nVrtx+1), "isout");
    vwType *weight = (vwType*) calloc(G->nVrtx+1, sizeof(vwType));
    int *matchedu = (int*) calloc(G->nVrtx+1, sizeof(int));
    int *matchedv = (int*) calloc(G->nVrtx+1, sizeof(int));
    idxType *nbbadoutput = (idxType *) calloc(G->nVrtx + 1, sizeof(idxType));
    idxType *badoutputnode = (idxType *) calloc(G->nVrtx + 1, sizeof(idxType));
    idxType *nbbadinput = (idxType *) calloc(G->nVrtx + 1, sizeof(idxType));
    idxType *badinputnode = (idxType *) calloc(G->nVrtx + 1, sizeof(idxType));

    if (weight == NULL || matchedv == NULL || matchedu == NULL || nbbadinput == NULL || nbbadoutput == NULL || badoutputnode == NULL || badinputnode == NULL)
        u_errexit("Failed allocation in clusteringCoTop\n");
    ecType* edgesweight;
    switch (opt.co_eorder) {
        case CO_EORD_AGGWEIGHT :
        case CO_EORD_AGGRATIO :
            edgesweight = (ecType*) calloc(G->nVrtx+1, sizeof(ecType));
            break;
    }

    for (i=1; i<=G->nVrtx; i++) {
        leader[i] = i;
        weight[i] = G->vw[i];
    }
    t1 += u_wseconds();

    double t2  = -u_wseconds();
    for (n=1; n<=G->nVrtx; n++) {
        i = norder[n];
        if ((matchedu[i] == 1)||(matchedv[i] == 1))
            continue;

        int nbneighbors = lookNeighborsAggForbidden(C, i, neighidx, isout, toplevels, leader, matchedu, matchedv, nbbadinput, badinputnode, nbbadoutput, badoutputnode, weight, opt);

        if (nbneighbors == 0)
            continue;

        idxType node;
        int nodeisout;
        selectBestNeighborsAgg(G, &node, &nodeisout, neighidx, nbneighbors, isout, leader, weight, edgesweight, opt);

        leader[i] = leader[node];
        weight[leader[node]] += G->vw[i];
        if (nodeisout == 1) {
            matchedu[i] = 1;
            if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT))
                for (k = G->outStart[i]; k <= G->outEnd[i]; k++) {
                    idxType child = G->out[k];
                    if (abs(toplevels[child] - toplevels[i]) <= 1) {
                        if ((nbbadinput[child] == 0) || (leader[badinputnode[child]] != leader[i]))
                            nbbadinput[child]++;
                        badinputnode[child] = i;
                    }
                }
            if (matchedv[node] == 0) {
                matchedv[node] = 1;
                if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT))
                    for (k = G->inStart[node]; k <= G->inEnd[node]; k++) {
                        idxType dad = G->in[k];
                        if (abs(toplevels[dad] - toplevels[node]) <= 1) {
                            if ((nbbadoutput[dad] == 0) || (leader[badoutputnode[dad]] != leader[node]))
                                nbbadoutput[dad]++;
                            badoutputnode[dad] = node;
                        }
                    }
            }
        }
        if (nodeisout == 0) {
            matchedv[i] = 1;
            if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT))
                for (k = G->inStart[i]; k <= G->inEnd[i]; k++) {
                    idxType dad = G->in[k];
                    if (abs(toplevels[dad] - toplevels[i]) <= 1) {
                        if ((nbbadoutput[dad] == 0) || (leader[badoutputnode[dad]] != leader[i]))
                            nbbadoutput[dad]++;
                        badoutputnode[dad] = i;
                    }
                }
            if (matchedu[node] == 0) {
                matchedu[node] = 1;
                if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT))
                    for (k = G->outStart[node]; k <= G->outEnd[node]; k++) {
                        idxType child = G->out[k];
                        if (abs(toplevels[child] - toplevels[node]) <= 1) {
                            if ((nbbadinput[child] == 0) || (leader[badinputnode[child]] != leader[node]))
                                nbbadinput[child]++;
                            badinputnode[child] = node;
                        }
                    }
            }
        }
    }
    t2 += u_wseconds();

    free(matchedu);
    free(matchedv);
    free(weight);
    free(isout);
    free(neighidx);
    free(nbbadoutput);
    free(badoutputnode);
    free(nbbadinput);
    free(badinputnode);

    switch (opt.co_eorder) {
        case CO_EORD_AGGWEIGHT :
        case CO_EORD_AGGRATIO :
            free(edgesweight);
            break;
    }
}


int createCycleDown(dgraph* G, idxType currentnode, idxType* leader, idxType targetedleader, int* hasbeentraversed, idxType* traversed, idxType* nbtraversed, idxType* toplevels, int weleft);

int createCycleUp(dgraph* G, idxType currentnode, idxType* leader, idxType targetedleader, int* hasbeentraversed, idxType* traversed, idxType* nbtraversed, idxType* toplevels, int weleft){
    idxType j,k;
    for (j = G->inStart[currentnode]; j <= G->inEnd[currentnode]; j++) {
        idxType innode = G->in[j];
        if (leader[innode] != leader[currentnode])
            continue;
        if (abs(toplevels[innode] - toplevels[currentnode]) > 1)
            continue;
        if ((leader[innode] == targetedleader)&&(weleft)){
            return 1;
        }
        if (hasbeentraversed[innode] == 1)
            continue;
        hasbeentraversed[innode] = 1;
        traversed[*nbtraversed] = innode;
	    (*nbtraversed)++;
        int cycle = createCycleDown(G, innode, leader, targetedleader, hasbeentraversed, traversed, nbtraversed, toplevels, leader[innode] != targetedleader);
        if (cycle == 1){
            return 1;
        }
    }
    return 0;
}

int createCycleDown(dgraph* G, idxType currentnode, idxType* leader, idxType targetedleader, int* hasbeentraversed, idxType* traversed, idxType* nbtraversed, idxType* toplevels, int weleft){
    idxType j,k;
    for (j = G->outStart[currentnode]; j <= G->outEnd[currentnode]; j++) {
        idxType outnode = G->out[j];
        if (abs(toplevels[outnode] - toplevels[currentnode]) > 1)
            continue;
        if ((leader[outnode] == targetedleader)&&(weleft)) {
            return 1;
        }
        if (hasbeentraversed[outnode] == 1)
            continue;
        hasbeentraversed[outnode] = 1;
        traversed[*nbtraversed] = outnode;
        (*nbtraversed)++;
        int cycle = createCycleUp(G, outnode, leader, targetedleader, hasbeentraversed, traversed, nbtraversed, toplevels, leader[outnode] != targetedleader);
        if (cycle == 1){
            return 1;
        }
    }
    return 0;
}


void clusteringCoCyc(coarsen* C, idxType *toplevels, idxType *norder, idxType* leader, const MLGP_option opt)
{
    dgraph* G = C->graph;
    idxType i, n, k, j, nbtraversed = 0;

    double t1 = -u_wseconds();

    idxType *neighidx = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "neighidx");
    int *isout = (idxType*) umalloc(sizeof(int)*(G->nVrtx+1), "isout");
    vwType *weight = (vwType*) calloc(G->nVrtx+1, sizeof(vwType));
    int *matchedu = (int*) calloc(G->nVrtx+1, sizeof(int));
    int *matchedv = (int*) calloc(G->nVrtx+1, sizeof(int));
    idxType *traversed = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "traversed");
    int *hasbeentraversed = (int*) calloc(G->nVrtx+1, sizeof(int));

    if (matchedu == NULL || matchedv == NULL || hasbeentraversed == NULL || weight == NULL)
        u_errexit ("failed allocation in clusteringCoCyc\n");
    ecType* edgesweight;
    switch (opt.co_eorder) {
        case CO_EORD_AGGWEIGHT :
        case CO_EORD_AGGRATIO :
            edgesweight = (ecType*) calloc(G->nVrtx+1, sizeof(ecType));
            break;
    }

    for (i=1; i<=G->nVrtx; i++) {
        leader[i] = i;
        weight[i] = G->vw[i];
    }
    t1 += u_wseconds();
    double t2  = -u_wseconds();
    for (n=1; n<=G->nVrtx; n++) {
        i = norder[n];
        if ((matchedu[i] == 1) || (matchedv[i] == 1))
            continue;

        int nbneighbors = lookNeighborsAggCheck(C, i, neighidx, isout, toplevels, leader, matchedu, matchedv, weight,
                                                opt);

        if (nbneighbors == 0)
            continue;

        idxType node;
        int nodeisout;

        selectBestNeighborsAgg(G, &node, &nodeisout, neighidx, nbneighbors, isout, leader, weight, edgesweight, opt);

        int cycle = 0;
        leader[i] = leader[node];

        if (nodeisout == 0){
            if (createCycleUp(G, i, leader, leader[node], hasbeentraversed, traversed, &nbtraversed, toplevels, 0))
                cycle = 1;
        }else {
            if (createCycleDown(G, i, leader, leader[node], hasbeentraversed, traversed, &nbtraversed, toplevels, 0))
                cycle = 1;
        }

        for (j = 0; j < nbtraversed; j++)
            hasbeentraversed[traversed[j]] = 0;
        nbtraversed = 0;

        if(cycle){leader[i] = i; continue;};

        weight[leader[node]] += G->vw[i];
        if (nodeisout == 1){
            matchedu[i] = 1;
            matchedv[node] = 1;
        }else{
            matchedv[i] = 1;
            matchedu[node] = 1;
        }
    }
    t2 += u_wseconds();

    free(matchedu);
    free(matchedv);
    free(weight);
    free(isout);
    free(neighidx);
    free(traversed);
    free(hasbeentraversed);

    switch (opt.co_eorder) {
        case CO_EORD_AGGWEIGHT :
        case CO_EORD_AGGRATIO :
            free(edgesweight);
            break;
    }

}

int isBig(dgraph* G, idxType node){
    idxType i;
    idxType limit = (idxType) ceil(0.5 * sqrt(1.0*G->nVrtx));
    idxType size_limit = 16 > limit ? 16 : limit;
    size_limit = size_limit < G->nVrtx ? size_limit : G->nVrtx;

    if (G->outEnd[node] - G->outStart[node] + 1 >= size_limit)
        return 1;

    if (G->inEnd[node] - G->inStart[node] + 1 >= size_limit)
        return 1;

    for (i = G->inStart[node]; i <= G->inEnd[node]; i++) {
        if (G->outEnd[G->in[i]] - G->outStart[G->in[i]] + 1 >= size_limit)
            return 1;
    }
    for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
        if (G->inEnd[G->out[i]] - G->inStart[G->out[i]] + 1 >= size_limit)
            return 1;
    }

    return 0;
}

void clusteringCoHyb(coarsen* C, idxType *toplevels, idxType *norder, idxType* leader, const MLGP_option opt)
{
    dgraph* G = C->graph;
    idxType i, n, k, j, nbtraversed = 0;

    double t1 = -u_wseconds();

    idxType *neighidx = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "neighidx");
    int *isout = (idxType*) umalloc(sizeof(int)*(G->nVrtx+1), "isout");
    vwType *weight = (vwType*) calloc(G->nVrtx+1, sizeof(vwType));
    int *matchedu = (int*) calloc(G->nVrtx+1, sizeof(int));
    int *matchedv = (int*) calloc(G->nVrtx+1, sizeof(int));
    idxType *nbbadoutput = (idxType *) calloc(G->nVrtx + 1, sizeof(idxType));
    idxType *badoutputnode = (idxType *) calloc(G->nVrtx + 1, sizeof(idxType));
    idxType *nbbadinput = (idxType *) calloc(G->nVrtx + 1, sizeof(idxType));
    idxType *badinputnode = (idxType *) calloc(G->nVrtx + 1, sizeof(idxType));
    idxType *traversed = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "traversed");
    int *hasbeentraversed = (int*) calloc(G->nVrtx+1, sizeof(int));

    if (weight == NULL || matchedu == NULL || matchedv == NULL || nbbadinput == NULL || nbbadoutput == NULL
        || badinputnode == NULL || badoutputnode == NULL || hasbeentraversed == NULL)
        u_errexit ("Failed allocation in clusteringCoHyb\n");
    ecType* edgesweight;
    switch (opt.co_eorder) {
        case CO_EORD_AGGWEIGHT :
        case CO_EORD_AGGRATIO :
            edgesweight = (ecType*) calloc(G->nVrtx+1, sizeof(ecType));
            break;
    }

    for (i=1; i<=G->nVrtx; i++) {
        leader[i] = i;
        weight[i] = G->vw[i];
    }
    t1 += u_wseconds();

    double t2  = -u_wseconds();
    for (n=1; n<=G->nVrtx; n++) {
        i = norder[n];
        if ((matchedu[i] == 1)||(matchedv[i] == 1))
            continue;

        int isbig = isBig(G, i);

        int nbneighbors;

        if (isbig == 1) {
            nbneighbors = lookNeighborsAggForbidden(C, i, neighidx, isout, toplevels, leader, matchedu, matchedv,
                                                    nbbadinput, badinputnode, nbbadoutput, badoutputnode, weight, opt);
        }else
            nbneighbors = lookNeighborsAggCheck(C, i, neighidx, isout, toplevels, leader, matchedu, matchedv, weight, opt);

        if (nbneighbors == 0)
            continue;

        idxType node;
        int nodeisout;
        selectBestNeighborsAgg(G, &node, &nodeisout, neighidx, nbneighbors, isout, leader, weight, edgesweight, opt);

        int cycle = 0;

        leader[i] = leader[node];

        if (isbig == 0){
            if (nodeisout == 0){
                if (createCycleUp(G, i, leader, leader[node], hasbeentraversed, traversed, &nbtraversed, toplevels, 0))
                    cycle = 1;
            }else {
                if (createCycleDown(G, i, leader, leader[node], hasbeentraversed, traversed, &nbtraversed, toplevels, 0))
                    cycle = 1;
            }

            for (j = 0; j < nbtraversed; j++)
                hasbeentraversed[traversed[j]] = 0;
            nbtraversed = 0;

            if(cycle){leader[i] = i; continue;};
        }

        weight[leader[node]] += G->vw[i];
        if (nodeisout == 1) {
            matchedu[i] = 1;
            if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT))
                for (k = G->outStart[i]; k <= G->outEnd[i]; k++) {
                    idxType child = G->out[k];
                    if (abs(toplevels[child] - toplevels[i]) <= 1) {
                        if ((nbbadinput[child] == 0) || (leader[badinputnode[child]] != leader[i]))
                            nbbadinput[child]++;
                        badinputnode[child] = i;
                    }
                }
            if (matchedv[node] == 0) {
                matchedv[node] = 1;
                if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT))
                    for (k = G->inStart[node]; k <= G->inEnd[node]; k++) {
                        idxType dad = G->in[k];
                        if (abs(toplevels[dad] - toplevels[node]) <= 1) {
                            if ((nbbadoutput[dad] == 0) || (leader[badoutputnode[dad]] != leader[node]))
                                nbbadoutput[dad]++;
                            badoutputnode[dad] = node;
                        }
                    }
            }
        }
        if (nodeisout == 0) {
            matchedv[i] = 1;
            if ((opt.co_dir == CO_DIR_OUT) || (opt.co_dir == CO_DIR_INOUT))
                for (k = G->inStart[i]; k <= G->inEnd[i]; k++) {
                    idxType dad = G->in[k];
                    if (abs(toplevels[dad] - toplevels[i]) <= 1) {
                        if ((nbbadoutput[dad] == 0) || (leader[badoutputnode[dad]] != leader[i]))
                            nbbadoutput[dad]++;
                        badoutputnode[dad] = i;
                    }
                }
            if (matchedu[node] == 0) {
                matchedu[node] = 1;
                if ((opt.co_dir == CO_DIR_IN) || (opt.co_dir == CO_DIR_INOUT))
                    for (k = G->outStart[node]; k <= G->outEnd[node]; k++) {
                        idxType child = G->out[k];
                        if (abs(toplevels[child] - toplevels[node]) <= 1) {
                            if ((nbbadinput[child] == 0) || (leader[badinputnode[child]] != leader[node]))
                                nbbadinput[child]++;
                            badinputnode[child] = node;
                        }
                    }
            }
        }
    }
    t2 += u_wseconds();

    free(matchedu);
    free(matchedv);
    free(weight);
    free(isout);
    free(neighidx);
    free(nbbadoutput);
    free(badoutputnode);
    free(nbbadinput);
    free(badinputnode);
    free(traversed);
    free(hasbeentraversed);

    switch (opt.co_eorder) {
        case CO_EORD_AGGWEIGHT :
        case CO_EORD_AGGRATIO :
            free(edgesweight);
            break;
    }
}

void matchIsolated(coarsen* C, idxType* norder, idxType* leader, const MLGP_option opt){
    //This function matched the isolated nodes two by two
    //We asume that thay have not been matched with anyone yet
    dgraph* G = C->graph;
    idxType i, vi, flag;
    double weight_lim = 0.1 * G->totvw / opt.nbPart;
    idxType* previous_isolated = (idxType*) umalloc((opt.nbPart+1) * sizeof(idxType), "previous_isolated");

    for (i=0; i<= opt.nbPart+1; i++)
        previous_isolated[i] = -1;

    for (vi = 1; vi <= G->nVrtx; vi++){
        i = norder[vi];
        flag = C->flag[i];
        if (G->inStart[i] != G->inEnd[i] + 1)
            continue;
        if (G->outStart[i] != G->outEnd[i] + 1)
            continue;
        if (leader[i] != i)
            u_errexit("matchIsolated: somebody touched the isolated node %d before!", i);
        if (previous_isolated[flag] == -1)
            previous_isolated[flag] = i;
        else{
            if (G->vw[i] + G->vw[previous_isolated[flag]] > weight_lim){
                previous_isolated[flag] = G->vw[i] <= G->vw[previous_isolated[flag]] ? i : previous_isolated[flag];
            } else {
                leader[i] = previous_isolated[flag];
                previous_isolated[flag] = -1;
            }
        }
    }
}

void computeClustering(const MLGP_option opt, int co_level, coarsen *C, idxType* matching)
{
    dgraph* G = C->graph;
    idxType* norder = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx + 1), "computeClustering: norder");
    computeMatchNOrder(G, opt, norder);

    idxType* levels = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx + 1), "levels");
    switch(opt.co_match_level){
        case CO_LEVEL_TL :
            computeToplevels(G, levels);
            break;

        case CO_LEVEL_BL :
            computeBottomlevels(G, levels);
            break;

        case CO_LEVEL_TLBL :
            if (co_level % 2 == 1) {
                computeToplevels(G, levels);
            }
            else {
                computeBottomlevels(G, levels);
            }
            break;

        default:
            u_errexit("computeClustering: unknown match_level = %d", opt.co_match_level);
    }

    switch(opt.co_match){
        case CO_ALG_MATCH_DIFF :
            matchingCoTop(C, levels, norder, matching, opt);
            break;

        case CO_ALG_AGG_DIFF :
            clusteringCoTop(C, levels, norder, matching, opt);
            break;

        case CO_ALG_AGG_CHECK :
            clusteringCoCyc(C, levels, norder, matching, opt);
            break;

        case CO_ALG_AGG_MIX :
            clusteringCoHyb(C, levels, norder, matching, opt);
            break;

        default:
            u_errexit("computeClustering: unknown matching = %d", opt.co_match);
    }

    if (opt.co_match_isolatedmatch)
        matchIsolated(C, norder, matching, opt);

    free(norder);
    free(levels);
}
