#include "utils.h"
#include "dgraph.h"
#include "dgraphReader.h"
#include "rvcycle.h"


int processArgs_rMLGP(int argc, char **argv, MLGP_option* opt)
{
    if (argc < 3) {
        printMLGPusage(argv[0]);
        u_errexit("%s: There is a problem with the arguments.\n", argv[0]);
    }

    initMLGPoptions(opt, atoi(argv[2]));
    processMLGPargs(argc, argv, opt);

    // if ((opt->nbPart<1) || !isPowOfTwo(opt->nbPart)) {
    if ((opt->nbPart<1)) {
        printMLGPusage(argv[0]);
        u_errexit("%s: This number of partitions is not supported for now. Try a power of 2.\n", argv[0]);
    }

    return 0;
}


void run_rMLGP(char* file_name, const MLGP_option opt)
{
    dgraph G;
    double graphRead =- u_wseconds();
    readDGraph(&G, file_name, opt.use_binary_input);
    graphRead += u_wseconds();
    set_dgraph_info(&G);
    if (opt.print > 0) {
        printf ("Reading input file: %.2lf seconds\n", graphRead);
    }

    int maxindegree, minindegree, maxoutdegree, minoutdegree;
    double aveindegree, aveoutdegree;
    dgraph_info(&G, &maxindegree, &minindegree, &aveindegree, &maxoutdegree, &minoutdegree, &aveoutdegree);
    G.maxindegree = maxindegree;
    G.maxoutdegree = maxoutdegree;

    int i;
    for (i=0; i<opt.nbPart; i++) {
        if (opt.lb[i] < 0)
            opt.lb[i] = 1;
        if (opt.ub[i] < 0)
            opt.ub[i] = opt.ratio * (double)G.totvw/(double)opt.nbPart; /*TODO:BU2JH: if no opt.ratio *, then first run uses ub as exact balance. Is that so?*/

        if ((floor(opt.lb[i]) < opt.lb[i])&&(floor(opt.lb[i]) == floor(opt.ub[i]))) {
            printf("WARNING: The balance can not be matched!!!\n");
            opt.lb[i] = floor(opt.lb[i]);
            opt.ub[i] = floor(opt.ub[i])+1;
        }
    }
    // enable if you wish to see the parameters in the run
    // printOptions(&opt, ';');

    printf("Graph Information:\n\tNb node: %d\n\tNb Edges: %d\n\tMax in-degree: %d\n\tMax out-degree: %d\n\tAv. in-degree: %.2f\n\tAv. out-degree: %.2f\nProblem Information:\n\tNb part: %d\n\tLower Bound[0]: %.1f\n\tUpper Bound[0]: %.1f\n\n", G.nVrtx, G.nEdge, maxindegree, maxoutdegree,aveindegree,aveoutdegree,opt.nbPart, opt.lb[0], opt.ub[0]);
    ecType * edgecut = (ecType *) malloc(sizeof(ecType) * opt.runs);
    int* nbcomm = (int*) malloc(sizeof(int) * opt.runs);
    double* latencies = (double*) malloc(sizeof(double) * opt.runs);
    int r;
    rcoarsen * rcoars;
    idxType nbcomp = 0;

    if(opt.seed == 0) {
        usRandom((int) time(NULL));
    }
    else
        usRandom(opt.seed);

    for (r = 0; r<opt.runs; r++) {
        int isAcyclic;
        rMLGP_info* info = (rMLGP_info*)  malloc (sizeof (rMLGP_info));
        initRInfoPart(info);
        printf("########################## RUN %d (seed=%d) ########################\n", r, uGetSeed());

        rcoars = rVCycle(&G, opt, info);

        edgecut[r] = edgeCut(&G, rcoars->coars->part);

        int maxsize = printPartWeights(&G, rcoars->coars->part);
        printf("Partition:\n\tEdgecut: %d\n\tBalance: %f\n\tVCycle depth: %d\n\tVertex Contraction: %.3f\n\tEdge Contraction: %.3f\n\tEdge Weight Contraction: %.3f\nTimes in seconds:\n\tCoarsening: %.3lf\n\tInitial Partition: %.3lf\n\tUncoarsening: %.3lf\n\tTotal: %.3lf\n", (int) edgecut[r], (double) maxsize / (G.totvw/opt.nbPart), info->info->coars_depth, (double) info->info->nbnodes_coars_tab[info->info->coars_depth] / (double) G.nVrtx, (double) info->info->nbedges_coars_tab[info->info->coars_depth] / (double) G.nEdge,  (double) ((double) (int) info->info->coarse_ew/ (double) G.nEdge),  info->timing_coars, info->timing_inipart, info->timing_uncoars,info->timing_global);


        if (opt.debug) {
            isAcyclic = checkAcyclicity(&G, rcoars->coars->part, opt.nbPart);
            if (isAcyclic == 0)
                u_errexit("rMLGP: the partition obtained is not acyclic\n");
        }

        // print the parts to file.
        if(opt.write_parts>0)
        {
            char name_tmp[200];
            sprintf(name_tmp,".partsfile.part_%d.seed_%d.txt", opt.nbPart, opt.seed);
            char res_name[200];
            strcpy(res_name, file_name);
            strcat(res_name, name_tmp);
            writePartsToTxt(&G,res_name,rcoars->coars->part);
            sprintf(name_tmp,".partitioned.part_%d.seed_%d.dot", opt.nbPart, opt.seed);
            strcpy(res_name, file_name);
            strcat(res_name, name_tmp);
            writeDGraphDot(&G, res_name, rcoars->coars->part);

        }

        printf("######################## END RUN %d ########################\n\n", r);
        if (opt.print > 0) {
            printf("######################## OUTPUT %d ########################\n\n", r);
            printRInfoPart(info, G.nVrtx, G.nEdge, opt.print);
            printf("\n###########################################################\n\n");
        }
        freeRCoarsenHeader(rcoars);
        rcoars = (rcoarsen*) NULL;
        freeRInfoPart(info);
    }

    ecType edgecutave = 0.0, nbcommave = 0.0, edgecutsd = 0.0, nbcommsd = 0.0;
    for (r = 0; r<opt.runs; r++) {
        edgecutave += edgecut[r];
        nbcommave += nbcomm[r];
    }
    edgecutave = edgecutave / opt.runs;
    nbcommave = nbcommave / opt.runs;
    for (r = 0; r<opt.runs; r++) {
        edgecutsd = (edgecut[r] - edgecutave) < 0 ? edgecutave - edgecut[r] : edgecut[r] - edgecutave;
        nbcommsd = (nbcomm[r] - nbcommave) < 0 ? nbcommave - nbcomm[r] : nbcomm[r] - nbcommave;
    }
    edgecutsd = edgecutsd / opt.runs;
    nbcommsd = nbcommsd / opt.runs;
    printf("Average Edgecut:%d\tStandard Deviation: %d\n", (int) edgecutave, (int) edgecutsd);

    free(edgecut);
    free(nbcomm);
    free(latencies);
    freeDGraphData(&G);
}

int main(int argc, char *argv[])
{
    MLGP_option opt;

    processArgs_rMLGP(argc, argv, &opt);
    run_rMLGP(opt.file_name, opt);

    free_opt(&opt);
    return 0;
}
