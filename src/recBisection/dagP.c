#include "utils.h"
// #include "dgraph.h"
#include "dgraphReader.h"
#include "rvcycle.h"
#include "dagP.h"

int dagP_init_parameters(MLGP_option* opt, const int nbPart) {
    initMLGPoptions(opt, nbPart);
    return 0;
}

int dagP_init_filename(MLGP_option* opt, char* file_name) {
    strncpy(opt->file_name, file_name, PATH_MAX);
    strncpy(opt->out_file, opt->file_name, PATH_MAX);
    strcat(opt->out_file, "_dagP.out");
    return 0;
}
int dagP_opt_reallocUBLB(MLGP_option* opt, const int nbPart) {

    free(opt->ub);
    free(opt->lb);
    opt->ub = (double*) umalloc((opt->nbPart+1)*sizeof(double), "opt->ub");
    opt->lb = (double*) umalloc((opt->nbPart+1)*sizeof(double), "opt->lb");
    int i=0;
    for (i=0; i<opt->nbPart; ++i) {
        opt->ub[i] = -1;
        opt->lb[i] = -1;
    }
    if (opt->ub == NULL || opt->lb == NULL)
        u_errexit("Could not reallocate ub, lb\n");
    return 0;
}

int dagP_read_graph(char* file_name, dgraph *G, const MLGP_option *opt) {
    readDGraph(G, file_name, opt->use_binary_input);
    set_dgraph_info(G);
    int maxindegree, minindegree, maxoutdegree, minoutdegree;
    double aveindegree, aveoutdegree;
    dgraph_info(G, &maxindegree, &minindegree, &aveindegree, &maxoutdegree, &minoutdegree, &aveoutdegree);
    G->maxindegree = maxindegree;
    G->maxoutdegree = maxoutdegree;
    return 0;
}
ecType dagP_partition_from_dgraph(dgraph* G, const MLGP_option *opt, idxType* parts) {
    set_dgraph_info(G);

    int maxindegree, minindegree, maxoutdegree, minoutdegree;
    double aveindegree, aveoutdegree;
    dgraph_info(G, &maxindegree, &minindegree, &aveindegree, &maxoutdegree, &minoutdegree, &aveoutdegree);
    G->maxindegree = maxindegree;
    G->maxoutdegree = maxoutdegree;

    int i;
    for (i=0; i<opt->nbPart; i++) {
        if (opt->lb[i] < 0)
            opt->lb[i] = 1;
        if (opt->ub[i] < 0)
            opt->ub[i] = opt->ratio * (double)G->totvw/(double)opt->nbPart; /*TODO:BU2JH: if no opt->ratio *, then first run uses ub as exact balance. Is that so?*/

        if ((floor(opt->lb[i]) < opt->lb[i])&&(floor(opt->lb[i]) == floor(opt->ub[i]))) {
            printf("WARNING: The balance can not be matched!!!\n");
            opt->lb[i] = floor(opt->lb[i]);
            opt->ub[i] = floor(opt->ub[i])+1;
        }
    }

    ecType * edgecut = (ecType *) malloc(sizeof(ecType) * opt->runs);
    ecType bestEdgecut = ecType_MAX;
    int* nbcomm = (int*) malloc(sizeof(int) * opt->runs);
    double* latencies = (double*) malloc(sizeof(double) * opt->runs);
    int r;
    rcoarsen * rcoars;
    idxType nbcomp = 0;

    if(opt->seed == 0) {
        usRandom((int) time(NULL));
    }
    else
        usRandom(opt->seed);

    for (r = 0; r<opt->runs; r++) {
        int isAcyclic;
        rMLGP_info* info = (rMLGP_info*)  malloc (sizeof (rMLGP_info));
        initRInfoPart(info);
        rcoars = rVCycle(G, *opt, info);
        edgecut[r] = edgeCut(G, rcoars->coars->part);
        int maxsize = printPartWeights(G, rcoars->coars->part);
        // printf("Partition:\n\tEdgecut: %d\n\tBalance: %f\n\tVCycle depth: %d\n\tVertex Contraction: %.3f\n\tEdge Contraction: %.3f\n\tEdge Weight Contraction: %.3f\nTimes in seconds:\n\tCoarsening: %.3lf\n\tInitial Partition: %.3lf\n\tUncoarsening: %.3lf\n\tTotal: %.3lf\n", (int) edgecut[r], (double) maxsize / (G->totvw/opt->nbPart), info->info->coars_depth, (double) info->info->nbnodes_coars_tab[info->info->coars_depth] / (double) G->nVrtx, (double) info->info->nbedges_coars_tab[info->info->coars_depth] / (double) G->nEdge,  (double) ((double) (int) info->info->coarse_ew/ (double) G->nEdge),  info->timing_coars, info->timing_inipart, info->timing_uncoars,info->timing_global);

        if (edgecut[r] < bestEdgecut) {
            bestEdgecut = edgecut[r];
            for (i=1; i <= G->nVrtx; ++i) {
                parts[i] = rcoars->coars->part[i];
            }
        }

        if (opt->print > 0) {
            printf("######################## OUTPUT %d ########################\n\n", r);
            printRInfoPart(info, G->nVrtx, G->nEdge, opt->print);
            printf("\n###########################################################\n\n");
        }
        freeRCoarsenHeader(rcoars);
        rcoars = (rcoarsen*) NULL;
        freeRInfoPart(info);
    }

    free(edgecut);
    free(nbcomm);
    free(latencies);
    return bestEdgecut;
}
int dagP_free_graph(dgraph* G) {
    freeDGraphData(G);
    return 0;
}

int dagP_free_option(MLGP_option* opt) {
    free_opt(opt);
    return 0;
}
