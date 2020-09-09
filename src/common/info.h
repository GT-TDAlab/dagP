#ifndef INFO_H_
#define INFO_H_

#include "option.h"
#include "utils.h"

#define IMBAL_IS_MAXIMBAL
#define INFO_LENGTH 100


typedef struct MLGP_info{
    double timing_global;
    double timing_coars;
    double timing_inipart;
    double timing_uncoars;
    double timing_conpar;

    double* timing_coars_tab;
    double* timing_matching_tab;
    double* timing_buildgraph_tab;

    double* timing_inipart_tab;
    ecType* ec_inipart_tab;

    double* timing_uncoars_tab;
    double* timing_refinement_tab;
    double* timing_project_tab;
    double* timing_forced_tab;

    ecType* ec_afterforced_tab;
    int* nbref_tab;
    ecType** ec_afterref_tab;
    double** timing_refinement_tab_tab;

    idxType* nbnodes_coars_tab;
    idxType* nbedges_coars_tab;

    int coars_depth;
    int best_inipart;
    ecType conpar_ec;
    ecType ini_ec;
    ecType final_ec;
    ecType coarse_ew;
    //some info used in the recursive bisection for debugging
    int depth;
    ecType current_edgecut;
} MLGP_info;

typedef struct rMLGP_info{
    double timing_global;
    double timing_coars;
    double timing_inipart;
    double timing_uncoars;

    MLGP_info* info;
    struct rMLGP_info* rec1;
    struct rMLGP_info* rec2;

    int depth;
} rMLGP_info;



void initInfoPart(MLGP_info* info);
void freeInfoPart(MLGP_info* info);
void printInfoPart(MLGP_info* info, int nVrtx, int nEdge, int printdetail);
void printRInfoPart(rMLGP_info* info, int nVrtx, int nEdge, int printdetail);

void initRInfoPart(rMLGP_info* info);
void freeRInfoPart(rMLGP_info* info);






#endif
