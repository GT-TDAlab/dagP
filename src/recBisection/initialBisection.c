#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "dgraph.h"
#include "utils.h"
#include "initialBisection.h"
#include "undirPartitioning.h"
#include "dgraphBisection.h"
#include "refinementBis.h"
#include "rvcycle.h"


void iniBisGGG(MLGP_option opt, coarsen*  coars, ecType *edgecut, idxType *partcopy, MLGP_info* info)
{
    ecType edgecut10, edgecut01;
    idxType i, *partcopy01, nVrtx = coars->graph->nVrtx;

    coars->nbpart = opt.nbPart;

    greedyGraphGrowingEqualTopOrder(coars->graph, opt.ub, partcopy, opt, 1, 0, &edgecut10);
    recBisRefinementStep(opt, info->coars_depth, coars->graph, partcopy, coars->nbpart, info);
    edgecut10 = edgeCut(coars->graph, partcopy);
    *edgecut = edgecut10;

    if (opt.debug)   {
        ecType ecutComputed =  edgeCut(coars->graph, partcopy);
        if (ecutComputed != edgecut10)
            u_errexit("iniPartGGG: cuts are not equal returned %d computed %d\n", (int) edgecut10, (int) ecutComputed);

    }
    partcopy01 = (idxType*) malloc(sizeof(idxType) * (nVrtx+1));

    reverseGraph(coars->graph);/*reverse the graph*/
    greedyGraphGrowingEqualTopOrder(coars->graph, opt.ub, partcopy01, opt, 0, 1, &edgecut01);
    recBisRefinementStep(opt, info->coars_depth, coars->graph, partcopy01, coars->nbpart, info);
    edgecut01 = edgeCut(coars->graph, partcopy01);

    if (opt.debug)   {
        ecType ecutComputed =  edgeCut(coars->graph, partcopy01);
        if (ecutComputed != edgecut01)
            u_errexit("iniPartGGG: cuts are not equal returned %d computed-R %d\n", (int) edgecut01, (int) ecutComputed);

    }

    reverseGraph(coars->graph);/*restore the graph*/
    if (edgecut01 < edgecut10) {
        *edgecut = edgecut01;
        for (i=1; i <= nVrtx; i++)
            partcopy[i] = partcopy01[i];
    }
    free(partcopy01);

}

void iniBisGGGTwoPhases(MLGP_option opt, coarsen*  coars, ecType *edgecut, idxType *partcopy, MLGP_info* info)
{
    ecType edgecut10, edgecut01;
    idxType i, *partcopy01, nVrtx = coars->graph->nVrtx;

    coars->nbpart = opt.nbPart;

    greedyGraphGrowingTwoPhases(coars->graph, opt.ub, partcopy, opt, 1, 0, &edgecut10);
    recBisRefinementStep(opt, info->coars_depth, coars->graph, partcopy, coars->nbpart, info);
    edgecut10 = edgeCut(coars->graph, partcopy);
    *edgecut = edgecut10;
    if (opt.debug)   {
        ecType ecutComputed =  edgeCut(coars->graph, partcopy);
        if (ecutComputed != edgecut10)
            u_errexit("iniPartGGG: cuts are not equal returned %.0f computed %.0f\n", (double) edgecut10, (double) ecutComputed);

    }

    partcopy01 = (idxType*) malloc(sizeof(idxType) * (nVrtx+1));

    reverseGraph(coars->graph);/*reverse the graph*/

    greedyGraphGrowingTwoPhases(coars->graph, opt.ub, partcopy01, opt, 0, 1, &edgecut01);
    recBisRefinementStep(opt, info->coars_depth, coars->graph, partcopy01, coars->nbpart, info);
    edgecut01 = edgeCut(coars->graph, partcopy01);

    if (opt.debug)   {
        ecType ecutComputed =  edgeCut(coars->graph, partcopy01);
        if (ecutComputed != edgecut01)
            u_errexit("iniPartGGG: cuts are not equal returned %.0f computed-R %.0f\n", (double) edgecut01, (double) ecutComputed);

    }
    reverseGraph(coars->graph);/*restore the graph*/

    if (edgecut01 < edgecut10) {
        *edgecut = edgecut01;
        for (i=1; i <= nVrtx; i++)
            partcopy[i] = partcopy01[i];
    }
    free(partcopy01);
}

void iniBisGGGIn(MLGP_option opt, coarsen*  coars, ecType *edgecut, idxType *partcopy, MLGP_info* info)
{
    ecType edgecut10, edgecut01;
    idxType i, *partcopy01, nVrtx = coars->graph->nVrtx;

    coars->nbpart = opt.nbPart;

    greedyGraphGrowingIn(coars->graph, opt.ub, partcopy, opt, 1, 0, &edgecut10);
    recBisRefinementStep(opt, info->coars_depth, coars->graph, partcopy, coars->nbpart, info);
    edgecut10 = edgeCut(coars->graph, partcopy);
    *edgecut = edgecut10;
    if (opt.debug)   {
        ecType ecutComputed =  edgeCut(coars->graph, partcopy);
        if (ecutComputed != edgecut10)
            u_errexit("iniPartGGG: cuts are not equal returned %.0f computed %.0f\n", (double) edgecut10, (double) ecutComputed);

    }
    partcopy01 = (idxType*) malloc(sizeof(idxType) * (nVrtx+1));

    reverseGraph(coars->graph);/*reverse the graph*/
    greedyGraphGrowingIn(coars->graph, opt.ub, partcopy01, opt, 0, 1, &edgecut01);
    recBisRefinementStep(opt, info->coars_depth, coars->graph, partcopy01, coars->nbpart, info);
    edgecut01 = edgeCut(coars->graph, partcopy01);

    if (opt.debug)   {
        ecType ecutComputed =  edgeCut(coars->graph, partcopy01);
        if (ecutComputed != edgecut01)
            u_errexit("iniPartGGG: cuts are not equal returned %.0f computed-R %.0f\n", (double) edgecut01, (double) ecutComputed);

    }

    reverseGraph(coars->graph);/*restore the graph*/
    if (edgecut01 < edgecut10) {
        *edgecut = edgecut01;
        for (i=1; i <= nVrtx; i++)
            partcopy[i] = partcopy01[i];
    }
    free(partcopy01);
}

void iniBisConstraintPartitioning(MLGP_option opt, coarsen*  coars, ecType *edgecut, idxType *partcopy, MLGP_info* info)
{
    coars->nbpart = opt.nbPart;
    idxType i;
    for (i = 1; i <= coars->graph->nVrtx; i++) {
        partcopy[i] = coars->flag[i];
    }
    recBisRefinementStep(opt, info->coars_depth, coars->graph, partcopy, coars->nbpart, info);
    *edgecut = edgeCut(coars->graph, partcopy);
}

void undirAcyclicBisectionLibrary(ecType* edgecut, coarsen* coars, idxType* part, MLGP_option opt, int seed, MLGP_info* info) {
    assert(opt.nbPart == 2);
    coars->nbpart = opt.nbPart;
    undirPartitioningLibrary(coars->graph, part, opt, seed);
    *edgecut = edgeCut(coars->graph, part);
    fixUndirBisection(edgecut, coars, part, opt, info);
}

void undirAcyclicBisectionRandom(ecType* edgecut, coarsen* coars, idxType* part, MLGP_option opt, MLGP_info* info) {
    assert(opt.nbPart ==2);
    coars->nbpart = opt.nbPart;
    undirBisectionRandom(coars->graph, part);
    fixUndirBisection(edgecut, coars, part, opt, info);
 }

void undirAcyclicBisectionGGG(ecType* edgecut, coarsen* coars, idxType* part, MLGP_option opt, MLGP_info* info) {
    assert(opt.nbPart ==2);
    coars->nbpart = opt.nbPart;
    undirBisectionGGG(coars->graph, part, opt);
    fixUndirBisection(edgecut, coars, part, opt, info);
}

void initialBisection(MLGP_option opt, coarsen*  coars, MLGP_info* info)
{
    int tryAll[] = {IP_UNDIR, IP_GGG_TWO};
    int ntryAll = sizeof(tryAll)/sizeof(tryAll[0]);


    idxType* partcopy = (idxType*) malloc(sizeof(idxType)*(coars->graph->nVrtx + 1));

    ecType edgecutmin = ecType_MAX, edgecut = ecType_MAX;

    int nruns = (opt.inipart == IP_CONPAR) ? 1 : opt.inipart_nrun;

    for (int run = 0; run < nruns; ++run) {
        int iniAlg = opt.inipart;

        if (iniAlg == IP_ALL )
            iniAlg = tryAll[run % ntryAll];

        info->timing_inipart_tab[iniAlg] -= u_wseconds();
        switch (iniAlg) {
            case IP_GGG:
                assert(opt.nbPart == 2);
                iniBisGGG(opt, coars, &edgecut, partcopy, info);
                break;

            case IP_UNDIR :
                assert(opt.nbPart == 2);
                undirAcyclicBisectionLibrary(&edgecut, coars, partcopy, opt, uGetSeed()+run, info);
                break;

            case IP_GGG_TWO:
                assert(opt.nbPart == 2);
                iniBisGGGTwoPhases(opt, coars, &edgecut, partcopy, info);
                break;

            case IP_GGG_IN:
                assert(opt.nbPart == 2);
                iniBisGGGIn(opt, coars, &edgecut, partcopy, info);
                break;

            case IP_CONPAR:
                iniBisConstraintPartitioning(opt, coars, &edgecut, partcopy, info);
                break;

            case IP_UNDIRRAND:
                assert(opt.nbPart == 2);
                undirAcyclicBisectionRandom(&edgecut, coars, partcopy, opt, info);
                break;

            case IP_UNDIRGGG:
                assert(opt.nbPart == 2);
                undirAcyclicBisectionGGG(&edgecut, coars, partcopy, opt, info);
                break;

            case IP_RECBISS:
            {
                coars->nbpart = opt.nbPart;
                MLGP_option* optini = copyOpt(opt, opt.nbPart);
                optini->inipart = IP_ALL;

                rMLGP_info infoini;
                infoini.depth=0;
                initRInfoPart(&infoini);

                rcoarsen * rcoars = rVCycle(coars->graph, *optini, &infoini);

                for (int j = 1; j <=coars->graph->nVrtx; j++) {
                    partcopy[j] = rcoars->coars->part[j];
                }

                free_opt(optini);
                edgecut = edgeCut(coars->graph, partcopy);
                /* UVC: TODO: there was no free for rcoars */
                break;
            }

            default:
                u_errexit("initialPartitioning: unknown initial partitioning = %d", opt.inipart);
        }
        info->timing_inipart_tab[iniAlg] += u_wseconds();

        if (opt.print>=PD_HIGH)
            printf("InitPart:%d:%s Cut=%d\n", run, iniAlgName[iniAlg], (int)edgecut);
        info->ec_inipart_tab[iniAlg] = edgecut;
        if (edgecut < edgecutmin) {
            edgecutmin = edgecut;
            info->best_inipart = iniAlg;
            info->ini_ec = edgecut;
            for (int j = 0; j <= coars->graph->nVrtx; ++j)
                coars->part[j] = partcopy[j];
        }
    }
    if (opt.print>=PD_MEDIUM)
        printf("InitPart:Final:%s Cut=%d\n", iniAlgName[opt.inipart], (int)edgecutmin);
    if (opt.debug){
        int isAcyclic = checkAcyclicity(coars->graph, coars->part, opt.nbPart);
        if(isAcyclic == 0)
            u_errexit("initialPartitioning: the initial partition is not acyclic\n");
        else
            printf("initialPartitioning: the initial partition is acyclic\n");
    }

    free(partcopy);
    info->current_edgecut = edgecutmin;
}
