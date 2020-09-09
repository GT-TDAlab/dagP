#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "initialBisection.h"
#include "vcycle2way.h"
#include "option.h"
#include "utils.h"
#include "dgraphReader.h"
#include "debug.h"
#include "clustering.h"
#include "refinementBis.h"
#include "info.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wunused-parameter"


coarsen *initializeCoarsen(dgraph* graph)
{
    coarsen * coars = (coarsen*) umalloc(sizeof(coarsen),"coars");
    coars->graph = graph;
    coars->previous_coarsen = NULL;
    coars->next_coarsen = NULL;
    coars->leader = (idxType * ) umalloc(sizeof(idxType) * (graph->nVrtx +1), "coars->leader");
    coars->new_index = (idxType * ) umalloc(sizeof(idxType) * (graph->nVrtx +1), "coars->new_index");
    coars->part = (idxType * ) umalloc(sizeof(idxType) * (graph->nVrtx +1), "coars->part");
    coars->flag = (idxType * ) umalloc((graph->nVrtx +1)*sizeof(idxType), "flag");
    coars->nbpart = 0;
    return coars;
}

/* don't call this with HEAD; or finest graph!!! */
void freeCoarsen(coarsen *coars)
{
    free(coars->leader);
    free(coars->new_index);
    free(coars->part);
    free(coars->flag);
    freeDGraphData(coars->graph);
    free(coars->graph);
    free(coars);
}

void freeCoarsenHeader(coarsen *coars)
{
    free(coars->leader);
    free(coars->new_index);
    free(coars->part);
    free(coars->flag);

    free(coars);
}


dgraph *initializeNextGraph(coarsen *coars, MLGP_option opt)
{
    /*Allocate and compute the next graph based on the matching in coars*/
    /*Modify coars->new_index accordingly*/
    idxType i,j,k;
    dgraph *next_graph = (dgraph*) umalloc(sizeof(dgraph), "initializeNextGraph: next_graph");

    dgraph* graph = coars->graph;
    idxType* leader = coars->leader;
    idxType nbnewnode = 0, nbnode = graph->nVrtx;
    idxType* newnodes = coars->new_index;
    double totalec;

    next_graph->vw = NULL;
    next_graph->ecIn = next_graph->ecOut = NULL;

    /*Compute new_index table*/
    for (i=1; i <= nbnode; i++) {
#if _dGP_DEBUG >= PD_HIGH
        if (opt.debug >= PD_HIGH )
            if (leader[leader[i]] != leader[i])
                u_errexit("initializeNextGraph: node %d is the leader of node %d and not its own leader\n", (int) leader[i], (int) i);
#endif
        if (leader[i] == i)
            newnodes[i] = ++nbnewnode;
    }

    idxType *cluster = NULL, *xcluster = NULL, idx;

    cluster = (idxType *) umalloc(sizeof(idxType)*(nbnode+1), "cluster");
    xcluster = (idxType *) calloc(nbnewnode+2, sizeof(idxType));
    next_graph->vw = (vwType *) calloc(nbnewnode+1, sizeof(vwType));
    if (cluster == NULL || xcluster == NULL || next_graph->vw == NULL) {
        printf("%d %d %d\n",nbnode+1,nbnewnode+2,nbnewnode+1 );
        u_errexit("initializeNextGraph: cluster, xcluster, or next_graph->vw could not be initialized.\n");
    }

    for (i=1; i<= nbnode; ++i) {
        newnodes[i] = newnodes[leader[i]];
        ++xcluster[newnodes[i]];
    }

    for (i=1, idx=0; i<= nbnewnode; ++i) {
        idxType sz=xcluster[i];
        xcluster[i] = idx;
        idx += sz;
    }
    xcluster[i] = idx;

    next_graph->maxVW = 0;
    for (i=1; i<= nbnode; ++i) {
        idxType ni=newnodes[i];
        cluster[xcluster[ni]++] = i;
        next_graph->vw[ni] += graph->vw[i];
        next_graph->maxVW = next_graph->maxVW < next_graph->vw[ni] ? next_graph->vw[ni] : next_graph->maxVW;
    }

    for (i=nbnewnode; i>0; --i) {
        xcluster[i] = xcluster[i-1];
    }

    next_graph->inStart = next_graph->inEnd = next_graph->in = NULL;
    next_graph->ecIn = NULL;
    next_graph->inStart = (idxType*) umalloc (sizeof(idxType) * (nbnewnode + 2),"next_graph->inStart" );
    next_graph->inEnd = (idxType*) umalloc (sizeof(idxType) * (nbnewnode + 2),"next_graph->inEnd" );
    next_graph->in = (idxType*) umalloc (sizeof(idxType) * (graph->nEdge+1),"next_graph->in" );
    next_graph->ecIn = (ecType*) umalloc (sizeof(ecType) * (graph->nEdge+1),"next_graph->ecIn" );

    next_graph->sources = (idxType *) umalloc(sizeof(idxType) * (nbnewnode+2), "next_graph->sources");
    next_graph->targets = (idxType *) umalloc(sizeof(idxType) * (nbnewnode+2), "next_graph->targets");
    if (next_graph->inStart == NULL || next_graph->inEnd == NULL || next_graph->in == NULL || next_graph->ecIn == NULL)
        u_errexit("vcycle.c: next_graph->inStart, next_graph->inEnd, next_graph->in, or next_graph->ecIn could not be initialized.\n");
    next_graph->maxindegree = 0;

    idxType *marker = NULL;
    ecType* work = NULL;
    idxType *in=next_graph->in;
    ecType *ecIn=next_graph->ecIn;

    work = (ecType *) malloc(sizeof(ecType) * (nbnewnode+1));
    marker = (int *) calloc(nbnewnode+1, sizeof(idxType));

    if (work == NULL || marker == NULL )
        u_errexit("initializeNextGraph: inedgeidx, work, marker could not be initialized.\n");

    next_graph->maxEC = 0;
    totalec = 0;
    for (i=1, idx=0; i <= nbnewnode; ++i) {
        next_graph->inStart[i] = idx;
        for (j = xcluster[i]; j < xcluster[i+1]; ++j) {
            idxType node = cluster[j];
            for (k=graph->inStart[node]; k <= graph->inEnd[node]; ++k) {
                idxType pred_node = newnodes[graph->in[k]];
                if (pred_node != i) { /*this is an in-coming neighbor*/
                    if (marker[pred_node] != i) {/*this is a new in-coming neighbor*/
                        marker[pred_node] = i;
                        work[pred_node] = 0;
                        in[idx++] = pred_node;
                    }
                    work[pred_node] += graph->ecIn[k];
                }
            }
        }
        for (k = next_graph->inStart[i]; k < idx; ++k) {
            ecIn[k] = work[in[k]];
            totalec += work[in[k]];
            next_graph->maxEC = next_graph->maxEC < ecIn[k] ? ecIn[k] : next_graph->maxEC;
        }

        next_graph->inEnd[i] = idx-1;
        idxType degree = next_graph->inEnd[i] - next_graph->inStart[i] + 1;
        next_graph->maxindegree = next_graph->maxindegree < degree ? degree : next_graph->maxindegree;

    }
    int nbnewedge = idx;

    next_graph->inStart[nbnewnode+1] = nbnewedge;
    next_graph->inEnd[nbnewnode+1] = nbnewedge-1;

    next_graph->in = (idxType*) realloc (next_graph->in, sizeof(idxType) * (nbnewedge+1));
    next_graph->ecIn = (ecType *) realloc(next_graph->ecIn, (nbnewedge+1) * sizeof(ecType));
    next_graph->ecOut = (ecType *) umalloc(sizeof(ecType) * (nbnewedge+1),"next_graph->ecOut" );
    if (next_graph->in == NULL || next_graph->ecIn == NULL || next_graph->ecOut == NULL)
        u_errexit("vcycle.c: next_graph->in, next_graph->ecIn, or next_graph->ecOut could not be initialized.\n");

    /*Fill next graph*/
    next_graph->frmt = graph->frmt;
    next_graph->nVrtx = nbnewnode;
    next_graph->nEdge = nbnewedge;
    next_graph->totvw = graph->totvw;
    next_graph->totec = totalec;
    next_graph->hollow = next_graph->outStart = next_graph->outEnd = next_graph->out = NULL;
    next_graph->hollow  = (int * ) calloc((nbnewnode + 2),sizeof(int));
    next_graph->outStart  = (idxType * ) umalloc(sizeof(idxType) * (nbnewnode + 2),"next_graph->outStart" );
    next_graph->outEnd  = (idxType * ) umalloc(sizeof(idxType) * (nbnewnode + 2),"outEnd" );
    next_graph->out = (idxType * ) umalloc(sizeof(idxType) * (nbnewedge+1),"out" );
    next_graph->nbsources = 0;
    if (next_graph->hollow == NULL || next_graph->outStart == NULL || next_graph->outEnd == NULL || next_graph->out == NULL)
        u_errexit("vcycle.c: next_graph->hollow, next_graph->outStart, next_graph->outEnd, or next_graph->out could not be initialized.\n");
    fillOutFromIn(next_graph);

    next_graph->nbtargets = 0;
    for (i=1; i<=nbnewnode; i++){
        if (next_graph->inEnd[i] < next_graph->inStart[i])
            next_graph->sources[next_graph->nbsources++] = i;
        if (next_graph->outEnd[i] < next_graph->outStart[i])
            next_graph->targets[next_graph->nbtargets++] = i;
    }
    free(work);
    free(marker);

    free(cluster);
    free(xcluster);

    return next_graph;
}
#pragma GCC diagnostic pop

coarsen *projectBack(coarsen* coars)
{
    /*Compute the partition for coars given the one in coars->next_coarsen*/
    coarsen* prev_coars = coars->previous_coarsen;
    idxType i;

    prev_coars->nbpart = coars->nbpart;
    for (i=1; i<=prev_coars->graph->nVrtx; i++) {
        prev_coars->part[i] = coars->part[prev_coars->new_index[i]];
    }
    return prev_coars;
}


coarsen *VCycle2way(dgraph *igraph, idxType* flag, const MLGP_option opt, MLGP_info* info)
{
    /*Proceed a vcycle starting from graph*/
    info->timing_global -= u_wseconds();
    idxType i;

    coarsen *current, *head = initializeCoarsen(igraph);

    int level = 0;
    idxType old_nbnode = igraph->nVrtx;

    current = head;

    for (i = 1; i <= igraph->nVrtx; i++){
        current->flag[i] = flag[i];
    }
    if (opt.nbPart == 1) {
        for (i = 1; i<= igraph->nVrtx; ++i){
            current->part[i] = 0;
        }
        return head;
    }
    info->timing_coars -= u_wseconds();
    while ( (current->graph->nVrtx > opt.co_stop_size) &&  (level < opt.co_stop_level)) {
        if (opt.print >= PD_HIGH) {
            dgraph *g=current->graph;
            printf("Coarsening:%2d: G(%6d, %6d)  maxDeg= <%3d, %3d>   maxW= [%4d, %d]\n",
                   level, g->nVrtx, g->nEdge, g->maxindegree, g->maxoutdegree, g->maxVW, (int) g->maxEC);
        }

        level++;
        info->timing_coars_tab[level] -= u_wseconds();
        info->timing_matching_tab[level] -= u_wseconds();
        computeClustering(opt, level, current, current->leader);
        idxType nodetoprint;
        dgraph* G = current->graph;

        info->timing_matching_tab[level] += u_wseconds();

        info->timing_buildgraph_tab[level] -= u_wseconds();
        dgraph* nextgraph = initializeNextGraph(current, opt);

        if (opt.toggle){
            printf("Ploting coars_%d.pdf ..... ", level);
            char fna[PATH_MAX];
            sprintf(fna, "plot/coars_%d.dot", level);
            dgraph_to_dot_with_matching(current->graph, current->leader, fna);
            sprintf(fna, "dot -Tpdf plot/coars_%d.dot -o plot/coars_%d.pdf", level, level);
            system(fna);
            printf("Done\n");
        }

        if (opt.debug >= PD_MEDIUM) {
            char fna[PATH_MAX];
            sprintf(fna, "plot/coars_%d.dot", level+1);
            dgraph_to_dot_with_matching(nextgraph, NULL, fna);

            printf("Check acyclicity depth %d level %d size %d edge %d\n", info->depth, level, nextgraph->nVrtx, nextgraph->nEdge);
            checkGraphAcyclicity(nextgraph);
            printf("Graph is acyclic\n");
        }

        info->nbnodes_coars_tab[level] = nextgraph->nVrtx;
        info->nbedges_coars_tab[level] = nextgraph->nEdge;
        coarsen* next = initializeCoarsen(nextgraph);
        for (i = 1; i <= G->nVrtx; i++)
            next->flag[current->new_index[i]] = current->flag[i];
        info->timing_buildgraph_tab[level] += u_wseconds();


        current->next_coarsen = next;
        next->previous_coarsen = current;

        if (opt.debug >= PD_HIGH) {
            char fn[PATH_MAX];
            sprintf(fn, "%s.CL.%02d.dot", opt.file_name, level);
            writeDGraphDot(next->graph, fn, NULL);
        }

        if (opt.debug > 100) {
            char fnn[100];

            sprintf(fnn,"matching_%d.txt", level);
            FILE *qq=fopen(fnn,"w");
            for (int i=1;i<=current->graph->nVrtx;++i) {
                if (i!=current->leader[i])
                    fprintf(qq, "%d %d\n", i, current->leader[i]);
            }
            fclose(qq);
        }

        current = next;
        info->timing_coars_tab[level] += u_wseconds();
        if ((double) current->graph->nVrtx / (double) old_nbnode > opt.co_stop_ratio)
            break;
        old_nbnode = current->graph->nVrtx;
    }

    if (opt.print >= PD_MEDIUM) {
        dgraph *g=current->graph;
        printf("Outside Coarsening Loop\nCoarsening:%2d: G(%6d, %6d)  maxDeg= <%3d, %3d>   maxW= [%4d, %d]\n",
               level, g->nVrtx, g->nEdge, g->maxindegree, g->maxoutdegree, g->maxVW, (int) g->maxEC);
    }
    info->coars_depth = level;
    info->timing_coars += u_wseconds();

    /*Initial partitioning*/
    info->timing_inipart -= u_wseconds();
    ecType coarsew = 0;
    for (int i=0; i<current->graph->nEdge; i++) {
        coarsew += current->graph->ecIn[i];
    }
    info->coarse_ew = coarsew;

    initialBisection(opt, current, info);

    info->timing_inipart += u_wseconds();

    if (opt.debug > 0)
        printPartWeights(current->graph, current->part);

    if (opt.toggle){
        printf("Ploting part_%d.pdf ..... ", level+1);
        char fna[PATH_MAX];
        sprintf(fna, "plot/part_%d.dot", level+1);
        dgraph_to_dot(current->graph, current->part, fna);
        sprintf(fna, "dot -Tpdf plot/part_%d.dot -o plot/part_%d.pdf", level+1, level+1);
        system(fna);
        printf("Done\n");
    }

    info->timing_uncoars -= u_wseconds();

    /* Uncoarsening */
    for ( ; level>=0 ;  level--) {
        info->timing_uncoars_tab[level] -= u_wseconds();

        info->timing_refinement_tab[level] -= u_wseconds();
        recBisRefinementStep(opt, level, current->graph, current->part, current->nbpart, info);
        info->timing_refinement_tab[level] += u_wseconds();

        /* UVC TODO: refinment should update edgecut/nbcombb, partsizes and we should print them without recompute*/
        if (opt.print >= PD_HIGH) {
            ecType edgecut = edgeCut(current->graph, current->part);
            printf("Refinment:%2d: G(%6d, %6d) Cut=%d\n", level, current->graph->nVrtx, current->graph->nEdge, (int) edgecut);
        }

        if (current->previous_coarsen) {
            info->timing_project_tab[level] -= u_wseconds();
            coarsen *prev = current;
            current = projectBack(current);
            freeCoarsen(prev);
            info->timing_project_tab[level] += u_wseconds();
        }

        if (opt.toggle){
            printf("Ploting part_%d.pdf ..... ", level);
            char fna[PATH_MAX];
            sprintf(fna, "plot/part_%d.dot", level);
            dgraph_to_dot(current->graph, current->part, fna);
            sprintf(fna, "dot -Tpdf plot/part_%d.dot -o plot/part_%d.pdf", level, level);
            system(fna);
            printf("Done\n");
        }

        if (opt.debug > 0)
            printPartWeights(current->graph, current->part);


        info->timing_uncoars_tab[level] += u_wseconds();
    }
    if (opt.print > PD_NONE) {
        info->final_ec = edgeCut(head->graph, head->part);
    }
    if (opt.print >= PD_MEDIUM) {
        printf("Refinment:Final: Cut=%d\n", (int) info->final_ec);
    }

    info->timing_uncoars += u_wseconds();
    info->timing_global += u_wseconds();

    return head;
}
