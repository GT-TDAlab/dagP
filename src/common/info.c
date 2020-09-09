#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "info.h"


int significantTime(double t, int printdetail)
{
    return ((t > 0.0051) || (printdetail > PD_HIGHEST));
}

void initInfoPart(MLGP_info* info)
{
    idxType i;

    info->timing_global = 0.0;
    info->timing_coars = 0.0;
    info->timing_inipart = 0.0;
    info->timing_uncoars = 0.0;
    info->timing_conpar = 0.0;

    info->timing_coars_tab = (double*) calloc(INFO_LENGTH, sizeof(double));
    info->timing_matching_tab = (double*) calloc(INFO_LENGTH, sizeof(double));
    info->timing_buildgraph_tab = (double*) calloc(INFO_LENGTH, sizeof(double));

    info->timing_inipart_tab = (double*) calloc(INFO_LENGTH, sizeof(double));
    info->ec_inipart_tab = (ecType*) calloc(INFO_LENGTH, sizeof(ecType));

    info->timing_uncoars_tab = (double*) calloc(INFO_LENGTH, sizeof(double));
    info->timing_refinement_tab = (double*) calloc(INFO_LENGTH, sizeof(double));
    info->timing_project_tab = (double*) calloc(INFO_LENGTH, sizeof(double));
    info->timing_forced_tab = (double*) calloc(INFO_LENGTH, sizeof(double));

    info->ec_afterforced_tab = (ecType*) calloc(INFO_LENGTH, sizeof(ecType));
    info->nbref_tab = (int*) calloc(INFO_LENGTH, sizeof(int));
    info->ec_afterref_tab = (ecType**) calloc(INFO_LENGTH, sizeof(ecType*));
    info->timing_refinement_tab_tab = (double**) calloc(INFO_LENGTH, sizeof(double*));
    for (i=0; i < INFO_LENGTH; i++) {
        info->ec_afterref_tab[i] = (ecType *) calloc(INFO_LENGTH, sizeof(ecType));
        info->timing_refinement_tab_tab[i] = (double *) calloc(INFO_LENGTH, sizeof(double));
    }

    info->nbnodes_coars_tab = (idxType*) calloc(INFO_LENGTH, sizeof(idxType));
    info->nbedges_coars_tab = (idxType*) calloc(INFO_LENGTH, sizeof(idxType));

    info->coars_depth = -1;
    info->best_inipart = -1;
    info->conpar_ec = -1;
    info->ini_ec = -1;
    info->final_ec = -1;
    info->current_edgecut = 0.0;

    info->depth = 0;
}

void freeInfoPart(MLGP_info* info)
{
    idxType i;

    free(info->timing_coars_tab);
    free(info->timing_matching_tab);
    free(info->timing_buildgraph_tab);

    free(info->timing_inipart_tab);
    free(info->ec_inipart_tab);

    free(info->timing_uncoars_tab);
    free(info->timing_refinement_tab);
    free(info->timing_project_tab);
    free(info->timing_forced_tab);

    free(info->ec_afterforced_tab);
    free(info->nbref_tab);
    for (i=0; i < INFO_LENGTH; i++) {
        free(info->ec_afterref_tab[i]);
        free(info->timing_refinement_tab_tab[i]);
    }
    free(info->timing_refinement_tab_tab);
    free(info->ec_afterref_tab);

    free(info->nbnodes_coars_tab);
    free(info->nbedges_coars_tab);
    free(info);
}

void initRInfoPart(rMLGP_info* info)
{
    info->timing_global = 0.0;
    info->timing_coars = 0.0;
    info->timing_inipart = 0.0;
    info->timing_uncoars = 0.0;
    info->depth = 0;

    info->info = (MLGP_info*) umalloc (sizeof (MLGP_info), "info->info");
    initInfoPart(info->info);
    info->rec1 = NULL;
    info->rec2 = NULL;
}

void freeRInfoPart(rMLGP_info* info)
{
    if (info->rec1 != NULL) {
        freeRInfoPart(info->rec1);
        info->rec1 = (rMLGP_info*) NULL;
    }
    if (info->rec2 != NULL) {
        freeRInfoPart(info->rec2);
        info->rec2 = (rMLGP_info *) NULL;
    }
    freeInfoPart(info->info);
    free(info);
}

void printInfoPart(MLGP_info* info, int nVrtx, int nEdge, int printdetail)
{
    int i,j;

    printf("Initial Graph : %d nodes %d edges\n", nVrtx, nEdge);
    printf("Final edgecut : %d \n", (int) info->final_ec);
    printf("Total time    : %.2f seconds\n", info->timing_global);

    if (significantTime(info->timing_conpar, 0)){
      printf("\nConPar\t :\t %.2f seconds (Edgecut %d)\n", info->timing_conpar, (int) info->conpar_ec);
    }

    printf("\nCoarsening [depth : %d]\t :\t %.2f seconds\t (%5.2f %%)\n", info->coars_depth, info->timing_coars, 100*info->timing_coars/info->timing_global);
    double timing_matching = 0.0, timing_buildgraph = 0.0;
    for (i=0; i < INFO_LENGTH; i++) {
        timing_matching += info->timing_matching_tab[i];
        timing_buildgraph += info->timing_buildgraph_tab[i];
    }
    printf("\tMatching\t :\t %.2f seconds\t (%5.2f %%)\n", timing_matching, 100*timing_matching/info->timing_coars);
    printf("\tBuilding Graph\t :\t %.2f seconds\t (%5.2f %%)\n", timing_buildgraph, 100*timing_buildgraph/info->timing_coars);
    if (printdetail > PD_HIGH) {
        for (i=1; i <= info->coars_depth; i++) {
            if (significantTime(info->timing_coars_tab[i], printdetail))
                printf("\t\tCoarsening Step %d\t :\t %.2f seconds (%5.2f %%)\t %d nodes %d edges\n", i, info->timing_coars_tab[i], 100*info->timing_coars_tab[i]/info->timing_coars, (int) info->nbnodes_coars_tab[i], (int) info->nbedges_coars_tab[i]);
        }
    }

    printf("Initial Partitioning\t :\t %.2f seconds\t (%5.2f %%)\n", info->timing_inipart, 100*info->timing_inipart/info->timing_global);

    if (printdetail > PD_HIGH) {
        for (i=0; i<IP_SIZE; ++i) {
            if (significantTime(info->timing_inipart_tab[i], printdetail)) {
                printf("\t%s           \t :\t %.2f seconds\t (%5.2f %%)\t(EdgeCut %d)\n", iniAlgName[i], info->timing_inipart_tab[i],
                        100 * info->timing_inipart_tab[i] / info->timing_inipart, (int) info->ec_inipart_tab[i]);
            }
        }
        printf("Best Initial Partition : %s with edgecut %d\n", iniAlgName[info->best_inipart], (int) info->ini_ec);
    }

    printf("Uncoarsening [depth : %d] :\t %.2f seconds\t (%5.2f %%)\n", info->coars_depth, info->timing_uncoars, 100*info->timing_uncoars/info->timing_global);
    double timing_refinement = 0.0, timing_project = 0.0, timing_forced = 0.0;
    for (i=0; i < INFO_LENGTH; i++) {
        timing_forced += info->timing_forced_tab[i];
        timing_project += info->timing_project_tab[i];
        for (j= 1; j <= info->nbref_tab[i]; j++)
            timing_refinement += info->timing_refinement_tab_tab[i][j];
    }
    printf("\tForce Balance\t :\t %.2f seconds\t (%5.2f %%)\n", timing_forced, 100*timing_forced/info->timing_uncoars);
    printf("\tRefinement\t :\t %.2f seconds\t (%5.2f %%)\n", timing_refinement, 100*timing_refinement/info->timing_uncoars);
    printf("\tProject Back\t :\t %.2f seconds\t (%5.2f %%)\n", timing_project, 100*timing_project/info->timing_uncoars);
    if (printdetail > PD_HIGH) {
        for (j = info->coars_depth; j >= 0; j--)
            if (significantTime(info->timing_uncoars_tab[j], printdetail)) {
                printf("\t\tUncoarsening Step %d\t :\t %.2f seconds (%5.2f %%)\n", j,
                       info->timing_uncoars_tab[j],
                       100 * info->timing_uncoars_tab[j] / info->timing_uncoars);
                if (significantTime(info->timing_forced_tab[j], printdetail))
                    printf("\t\t\tForce Balance pass\t :\t %.2f seconds (%5.2f %%)\t (edgecut %d)\n",
                       info->timing_forced_tab[j],
                       100 * info->timing_forced_tab[j] / info->timing_uncoars_tab[j],
                       (int) info->ec_afterforced_tab[j]);
                for (i = 1; i <= info->nbref_tab[j]; i++)
                    if (significantTime(info->timing_refinement_tab_tab[j][i], printdetail))
                        printf("\t\t\tRefinement pass %d\t :\t %.2f seconds (%5.2f %%)\t (edgecut %d)\n", i,
                           info->timing_refinement_tab_tab[j][i],
                           100 * info->timing_refinement_tab_tab[j][i] / info->timing_uncoars_tab[j],
                           (int) info->ec_afterref_tab[j][i]);

        }
    }
}

void printRInfoPart(rMLGP_info* info, int nVrtx, int nEdge, int printdetail)
{
    printInfoPart(info->info,  nVrtx,  nEdge, printdetail);
}



