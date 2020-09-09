#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>

#include "rvcycle.h"
#include "option.h"
#include "utils.h"
#include "undirPartitioning.h"
#include "initialBisection.h"

#ifdef dagP_SCOTCH
#include "scotch.h"
#endif


#if defined(dagP_SCOTCH) || defined(dagP_METIS)
void dgraph_to_metis(dgraph* G, idx_t* in, idx_t* adj, idx_t* vw, idx_t* adjw)
{
    idxType i,j,k=0;

    for (i=1; i<=G->nVrtx; i++) {
        in[i-1] = k;
        vw[i-1] = G->vw[i];
        for (j = G->inStart[i]; j <= G->inEnd[i]; j++) {
            adj[k] = G->in[j]-1;
            adjw[k++] = G->ecIn[j];
        }
        for (j = G->outStart[i]; j <= G->outEnd[i]; j++) {
            adj[k] = G->out[j]-1;
            adjw[k++] = G->ecOut[j];
        }
    }
    in[G->nVrtx] = k;
}
#endif

void undirPartitioningLibrary(dgraph* G, idxType* part, const MLGP_option opt, int seed)
{
    //This function partition a graph using Metis or Scotch with fixing the acyclicity
    //It can be used with more than 2 partitions
    //It is meant to be use in Constraint-Partitioning but re-used as initial partitioning.
#if defined(dagP_SCOTCH) || defined(dagP_METIS)

    idx_t nVertices = (idx_t) G->nVrtx;
    idx_t nEdges = (idx_t) G->nEdge;
    idx_t nWeights = (idx_t) 1;
    idx_t nParts = (idx_t) opt.nbPart;
    idxType* list = (idxType*) malloc((G->nVrtx+1)*sizeof(idxType));

    idx_t objval;
    idx_t* xadj = (idx_t*) calloc(nVertices+1, sizeof(idx_t));
    idx_t* vw = (idx_t*) calloc(nVertices+1, sizeof(idx_t));
    idx_t* adjncy = (idx_t*) calloc(2*nEdges, sizeof(idx_t));
    idx_t* adjw = (idx_t*) calloc(2*nEdges, sizeof(idx_t));
    idx_t* part_metis = (idx_t*) calloc(nVertices+1, sizeof(idx_t)) ;
    idxType i;
    double t0=0.0, t1=0.0;

    dgraph_to_metis(G, xadj, adjncy, vw, adjw);

    int ret = -1;
    switch(opt.undir_alg){
        case UNDIR_SCOTCH: ;
#ifndef dagP_SCOTCH
            u_errexit("Please compile with Scotch!");
#else
            SCOTCH_Graph graph;   /* Graph to part */
            SCOTCH_Strat strat;  /* Strategy to achieve partitionning */
                /* Initialize data structures. */
            SCOTCH_graphInit (&graph);
            SCOTCH_stratInit (&strat);
            if(SCOTCH_graphBuild(&graph, 0, nVertices, xadj, xadj + 1, vw, NULL,
                                 xadj[nVertices] , adjncy, adjw)){
                printf("Scotch: error in building graph\n");
                return;
            }

            /*printf("check scotch graph %d\n",SCOTCH_graphCheck(&graph));*/

            SCOTCH_stratGraphMapBuild(&strat, SCOTCH_STRATQUALITY, opt.nbPart, opt.ratio-1);

            if (opt.debug > PD_NONE)
                t0 = -u_wseconds();
            ret = SCOTCH_graphPart (&graph, opt.nbPart, &strat, part_metis);

            if (opt.debug > PD_NONE){
                t0 += u_wseconds();
                printf("\t\t\tUndirPartitioning: Scotch (undirected) time %.2f\n", objval, t0);

            }
            if (ret!=0) {
                fprintf(stderr, "Scotch failed with error code:%d\n", ret);
            }

            SCOTCH_graphExit (&graph);
            SCOTCH_stratExit (&strat);

            if (opt.print >=PD_LOW)
                  printf("\t\tJust scotch:\t\t%d\n",objval);
#endif
        break;
        case UNDIR_METIS: ;
#ifndef dagP_METIS
            u_errexit("Please compile with Metis");
#else
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions(options);
            options[METIS_OPTION_SEED] = (idx_t) seed;
            options[METIS_OPTION_NCUTS] = 1;
            real_t* targetPartWeights = (real_t*) umalloc (nParts * sizeof(real_t), "targetPartWeights");
            real_t sumub=0.0;
            for (i=0; i < nParts; ++i){
                sumub += opt.ub[i];
            }
            for(i=0; i < nParts; ++i){
                targetPartWeights[i] = opt.ub[i] / sumub;
            }

            t0 = -u_wseconds();
            ret = METIS_PartGraphKway(&nVertices, &nWeights, xadj, adjncy,
                                      vw, NULL, adjw, &nParts, targetPartWeights,
                                      NULL, options, &objval, part_metis);
            if (opt.print >=PD_LOW)
                printf("\t\tJust metis:\t\t%d\n",objval);

            if (opt.debug > PD_NONE){
                t0 += u_wseconds();
                printf("\t\t\tUndirPartitioning: Metis (undirected) cut=%d, _time %.2f\n", objval, t0);

            }

            if (ret!=METIS_OK) {
                fprintf(stderr, "Metis failed with error code:%d\n", ret);
            }
#endif
        break;
    }

    if (opt.debug > PD_NONE )
        t1 = -u_wseconds();

    for (i=1; i<=G->nVrtx; i++)
        part[i] = part_metis[i - 1];

    if (opt.debug > PD_NONE){
        t1 += u_wseconds();
        printf("\t\t\tUndirPartitioning (conpar): time %.2f\n", t1);
    }
    free(list);
    free(xadj);
    free(adjncy);
    free(vw);
    free(adjw);
    free(part_metis);
#else
    UNUSED(G); UNUSED(part); UNUSED(opt); UNUSED(seed);
    u_errexit("Function 'undirPartitioningLibrary' is called, but code is not linked with MeTiS or Scotch!\n");
#endif
}

void undirBisectionRandom(dgraph* graph, idxType* part)
{
    idxType i;
    for (i=1; i <= graph->nVrtx; i++){
        part[i] = uRandom(2);
    }
}

void undirBisectionGGG(dgraph* G, idxType* part, MLGP_option opt)
{
    idxType i, vi, j, neigh, nVrtx = G->nVrtx;
    double* ub_pw = opt.ub;
    idxType *sources = (idxType *) malloc((G->nVrtx + 1) * sizeof(idxType));
    idxType *targets = (idxType *) malloc((G->nVrtx + 1) * sizeof(idxType));
    idxType nbsources = sourcesList(G, sources);
    idxType nbtargets = outputList(G, targets);
    idxType *extrem = (idxType *) malloc(2*(G->nVrtx + 1) * sizeof(idxType));
    for (i=0; i < nbsources; i++)
        extrem[i] = sources[i];
    for (i=0; i < nbtargets; i++)
        extrem[i + nbsources] = targets[i];

    shuffleArray(nbsources+nbtargets, extrem, nbsources+nbtargets);

    idxType firstnodeidx = 0;
    idxType firstnode = extrem[firstnodeidx++];

    idxType fromPart, otherPart;
    fromPart = 0;
    otherPart = 1;

    idxType  *inEnd = G->inEnd,  *inStart = G->inStart, *in = G->in;
    idxType *outEnd = G->outEnd, *outStart = G->outStart, *out = G->out;
    ecType  *ecIn = G->ecIn, *ecOut = G->ecOut;

    vwType maxvw;
    idxType heapsz = 0; /*will be from 1 to heapsz*/
    idxType maxind;

    idxType *heapIds, *placedinZero, moved, bestAt, *visitOrder, *ismoved, *Degrees;
    int *inheap;
    ecType *Gains; /*this will be keyval in the heap*/
    ecType cut, bestCut=ecType_MAX;
    int isValidBsctn;

    vwType *vw = G->vw;
    double pw[2], bestCutBal;

    Gains = (ecType *) malloc(sizeof(ecType) * (nVrtx+1));
    Degrees = (idxType*) malloc(sizeof(idxType) * (nVrtx+1));
    heapIds = (idxType *) malloc(sizeof(idxType) * (nVrtx+1));
    inheap = (int *) malloc(sizeof(int) * (nVrtx+1));
    placedinZero = (idxType *) malloc(sizeof(idxType) * (nVrtx+1));
    ismoved = (idxType *) malloc(sizeof(idxType) * (nVrtx+1));
    visitOrder = (idxType *) malloc(sizeof(idxType) * (nVrtx+1));

    shuffleTab(1, nVrtx, visitOrder+1);/*this is the randomness */

    maxvw = -1;
    heapIds[0] = -1;
    pw[otherPart] = 0;
    pw[fromPart] = 0;

    for (vi = 1; vi <= nVrtx; vi++)
    {
        ecType in_out = 0;
        i = visitOrder[vi];

        part[i] = fromPart;

        if (vw[i] > maxvw) maxvw = vw[i];

        pw[fromPart] += vw[i];

        Degrees[i] = inEnd[i] - inStart[i] + 1 + outEnd[i] - outStart[i] + 1;

        Gains[i] = 0;
        inheap[i] = 0;
        ismoved[i] = 0;

        for (j=inStart[i]; j <= inEnd[i]; j++)
            Gains[i] -= ecIn[j];

        for (j=outStart[i]; j <= outEnd[i]; j++)
            Gains[i] -= ecOut[j];

    }
    bestCutBal = pw[fromPart];
    heapIds[++heapsz] = firstnode;
    inheap[firstnode] = 1;

    maxBuildHeapKeyvals(heapIds, heapsz, Gains, vw);


    bestAt = -1;
    moved = 0;
    cut = 0;
    isValidBsctn = 0;

    while (pw[otherPart] <= ub_pw[otherPart] + maxvw && pw[fromPart]>=opt.lb[fromPart] - maxvw)
    {
        if (heapsz == 0){
            while(ismoved[extrem[firstnodeidx]] == 1)
                firstnodeidx++;
            firstnode = extrem[firstnodeidx];
            heapIds[++heapsz] = firstnode;
            inheap[firstnode] = 1;
        }

        maxind = maxHeapExtractKeyvals(heapIds, &heapsz, Gains, vw);
        if (part[maxind] == otherPart)
            u_errexit("undirBisectionGGG: already moved there?\n");

        part[maxind] = otherPart;
        pw[otherPart] += vw[maxind];
        pw[fromPart] -= vw[maxind];

        placedinZero[moved++] = maxind;
        ismoved[maxind] = 1;
        cut -= Gains[maxind];
        for (j = inStart[maxind]; j<= inEnd[maxind]; j++){
            neigh = in[j];
            Gains[neigh] += 2*ecIn[j];
            if (inheap[neigh] == 0) {
                maxHeapInsertKeyvals(heapIds, &heapsz, neigh, Gains, vw);
                inheap[neigh] = 1;
            }
        }
        for (j = outStart[maxind]; j<= outEnd[maxind]; j++){
            neigh = out[j];
            Gains[neigh] += 2*ecOut[j];
            if (inheap[neigh] == 0) {
                maxHeapInsertKeyvals(heapIds, &heapsz, neigh, Gains, vw);
                inheap[neigh] = 1;
            }
        }
        if (pw[otherPart] <= ub_pw[otherPart]+maxvw && pw[fromPart] <= ub_pw[fromPart]+ maxvw) /*if both are less than their upper bound, this is a valid partition*/
        {
            double maxpw = mymax(pw[0], pw[1]);

            if ((isValidBsctn == 0) || ( isValidBsctn == 1 && (cut < bestCut  ||
                                                               (cut == bestCut && maxpw < bestCutBal ))))
            {
                bestCut = cut;
                bestAt = moved - 1;
                bestCutBal = maxpw;
            }
            isValidBsctn = 1;
        }
    }

    if (bestAt == -1)
        u_errexit("undirBisectionGGG: cannot find any balanced bipartition\n");
    for (i = bestAt+1; i < moved; i++)
        part[placedinZero[i]] = fromPart;

    free(sources);
    free(targets);
    free(extrem);
    free(placedinZero);
    free(Gains);
    free(Degrees);
    free(heapIds);
    free(inheap);
    free(ismoved);
    free(visitOrder);
}



void conpar(dgraph* G, idxType* constraint, const MLGP_option opt, MLGP_info* info){
    idxType i;
    info->timing_conpar -= u_wseconds();
    MLGP_option* optconpar;
    rMLGP_info* infoconpar;
    rcoarsen* rcoars;

    switch (opt.conpar_fix){
        case 0 :
            undirPartitioningLibrary(G, constraint, opt, uGetSeed());
            break;

        case 1 :
            optconpar = copyOpt(opt, opt.nbPart);
            infoconpar = (rMLGP_info*)  umalloc (sizeof (rMLGP_info), "infoconpar");
            initRInfoPart(infoconpar);

            optconpar->co_stop_level = 0;
            optconpar->conpar = 0;
            optconpar->inipart = IP_UNDIR;
            optconpar->inipart_nrun = 1;
            rcoars = rVCycle(G, *optconpar, infoconpar);
            for (i = 1; i <= G->nVrtx; i++)
                constraint[i] = rcoars->coars->part[i];
            info->conpar_ec = infoconpar->info->final_ec;

            free_opt(optconpar);
            free(optconpar);
            freeRCoarsenHeader(rcoars);
            rcoars = (rcoarsen*) NULL;
            freeRInfoPart(infoconpar);
            break;
    }
    info->timing_conpar += u_wseconds();

}
