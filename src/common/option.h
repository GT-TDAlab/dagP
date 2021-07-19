#ifndef OPTION_H_
#define OPTION_H_

#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <limits.h>

#include "dgraph.h"


#ifndef PATH_MAX
#define PATH_MAX 300
#endif


/* Print/Debug Details >= value means those output will be printed
   Use Print option to print program's progress with a clean output,
   Use Debug option to do debug checks and print additional debugging related output */
#define PD_NONE      0  /* lowest level of output, only final result  */
#define PD_LOW       1  /* Print: output of each run, Debug: Constant time debug output  */
#define PD_MEDIUM    2  /* Pring: only at coarsening, initial partitioning, refinement, recursive bisection level
                           Debug: same complexity with the current output  */
#define PD_HIGH      3  /* Print: at each level of coarsening, uncoarsening, or initial partitioning
                           Debug: higher complexity debugging  */
#define PD_HIGHEST   4  /* Print: print even non-significant times
                           Debug: highest complexity debugging  */

/* this is to enable/disable debug in compile time  */
#define _dGP_DEBUG   PD_NONE

#define CO_STOP_SIZE_MULT   50  /* default stop size k * CO_STOP_SIZE_MULT */
#define CO_STOP_LEVEL_DEF   20
#define CO_STOP_RATIO_DEF   0.9


#define CO_ALG_MATCH_DIFF       0  /* the matching is the one described in Coarsening 2 */
#define CO_ALG_AGG_DIFF         1  /* the matching is the aggregative with diff of TL 1 */
#define CO_ALG_AGG_CHECK        2  /* the matching is the aggregative with diff of TL 1 but we check cycle each time*/
#define CO_ALG_AGG_MIX          3  /* the matching is a mix between 2 and 1 (2 by default and 1 for nodes with large degree)*/

#define CO_ALG_DEF     CO_ALG_AGG_MIX

#define CO_DIR_IN       0  /* In the matching we only look at in-edges */
#define CO_DIR_OUT      1  /* In the matching we only look at out-edges */
#define CO_DIR_INOUT    2  /* In the matching we look at in-edges and out-edges */

#define CO_DIR_DEF      CO_DIR_INOUT

#define CO_LEVEL_TL             0  /* Matchings are based on top level values */
#define CO_LEVEL_BL             1  /* Matching are based on bottom level values */
#define CO_LEVEL_TLBL           2  /* Alternates the two */

#define CO_LEVEL_DEF    CO_LEVEL_TLBL

#define CO_NORD_RAND            0 /* Random traversal for the matching */
#define CO_NORD_RANDDFSTOP      1 /* random DFS top order traversal for the matching */
#define CO_NORD_RANDBFSTOP      2 /* random BFS top order traversal for the matching */
#define CO_NORD_RANDTOP         3 /* Random traversal top order for the matching */
#define CO_NORD_DFSTOP          4 /* DFS top order traversal for the matching */
#define CO_NORD_BFSTOP          5 /* BFS top order for the matching */
#define CO_NORD_RANDDFS         6 /* random DFS traversal for the matching */
#define CO_NORD_RANDBFS         7 /* random BFS traversal for the matching */
#define CO_NORD_DFS             8 /* DFS traversal for the matching */
#define CO_NORD_BFS             9 /* BFS traversal for the matching */

#define CO_NORD_DEF     CO_NORD_RANDBFSTOP

#define CO_EORD_NONE            0 /* Traverse the edges in the construction order */
#define CO_EORD_WEIGHT          1 /* Traverse the edges in the non-increasing order of their weight */
#define CO_EORD_RATIO           2 /* Traverse the edges in the non-increasing order of the ratio weight / node weight */
#define CO_EORD_AGGWEIGHT          3 /* Traverse the clusters in the non-increasing order of their sum edge weight */
#define CO_EORD_AGGRATIO           4 /* Traverse the clusters in the non-increasing order of the ratio sum edge weight / sum node weight */

#define CO_EORD_DEF     CO_EORD_RATIO

#define IP_ALL          0 /* try all initial partitioning and pick the best one */
#define IP_KERDFS       1 /* Kernighan algo with DFS traversal */
#define IP_KERBFS       2 /* Kernighan algo with BFS traversal */
#define IP_KERDFSBFS    3 /* Kernighan algo with mixed BFS-DFS traversal */
#define IP_GGG          4 /* greedy graph growing */
#define IP_FORCED       5 /* Forced Balance algo */
#define IP_UNDIR        6 /* Metis, then fixing acyclicity and forced balancing */
#define IP_BFS          7 /* BFS Growing */
#define IP_CS           8 /* CocktailShake */
#define IP_DIA_GGG      9 /* IP_GGG but starting with a node on the diameter of the graph */
#define IP_RECBISS      10 /* recursive bissectioh */
#define IP_GGG_TWO      11 /* greedy graph growing with two phases */
#define IP_GGG_IN       12 /* greedy graph growing with only in edges */
#define IP_CONPAR       13 /* initial partitioning computed in constraint partitioning */
#define IP_UNDIRRAND    14 /* undir random + fix */
#define IP_UNDIRGGG     15 /* undir GGG + fix */
#define IP_DEF          IP_CONPAR

#define IP_SIZE          16
extern char *iniAlgName[IP_SIZE];

#define CO_KERNIGHAN_DEFAULT 0 /*We use the regular Kernighan algo with increasing bound until we have the correct number of partition */
#define CO_KERNIGHAN_EXACT   1 /*We use the variant of Kernighan algo where exact number of component is ensured */

#define REF_NONE            0 /* No refinement */
#define REF_bFM             1 /* bisection FM */
#define REF_bFM_MAXW        2 /* bisection FM with vertex selected from maximum weighted part */
#define REF_KL              3 /* Kernighan-Lin (bisection with swap) */
#define REF_KL_bFM_MAXW     4 /* one pass Kernighan-Lin + one pass bFM_MAXW */
#define REF_kFM_PO          5 /* the refinement toward the greatest partition in post order */
#define REF_kFM_TB          6 /* the refinement toward the best partition if it creats no cycle */
#define REF_kFM_POTB        7 /* performes combination of refinement_po and then refinement_tb */
#define REF_kFM_RANDV       8 /* the refinement where the node to be moved to another partition is selected randomly among boundary vertices. */

#define REF_DEF         REF_bFM

#define CONPAR_DEF          1

#define CONPAR_FIX_DEF      1

#define ANCHORED_DEF        0

#define CO_OBJ_EC 0 /* edge-cut is the objective to minimize */
#define CO_OBJ_CV 1 /* communication volume is the objective to minimize */

#define CO_OBJ_COMMVOL_EXP      0 /* communication volume is computed using expected value */
#define CO_OBJ_COMMVOL_BEST     1 /* communication volume is computed using best case scenario */
#define CO_OBJ_COMMVOL_WORST    2 /* communication volume is computed using worst case scenario */
#define CO_OBJ_COMMVOL_AVE      3 /* communication volume is computed using average case scenario */

#define LIVE_CUT_MAX            0 /* Graph is cut is recursively cut in half until the live set value is smaller than opt.maxliveset */
#define LIVE_CUT_RECUT          1 /* We cut in half and recut if 2/3 1/3 if one part as a liveset smaller than opt.maxliveset / 2 */
#define UNDIR_METIS		0
#define UNDIR_SCOTCH		1

#define LIVE_TRAV_NAT            0 /* Traverse the nodes in the construction order to compute liveset */
#define LIVE_TRAV_DFS            1 /* Traverse the nodes in a DFS order to compute liveset */
#define LIVE_TRAV_BFS            2 /* Traverse the nodes in a BFS order to compute liveset */
#define LIVE_TRAV_MIX            3 /* Traverse the nodes in a mix DFS-BFS order to compute liveset */
#define LIVE_TRAV_RAND            4 /* Traverse the nodes in a random topological order to compute liveset */

#define LIVE_TRAV_DEF     LIVE_TRAV_NAT

typedef struct{
    char file_name[PATH_MAX];
    char file_name_dot[PATH_MAX];
    int nbPart;

    int co_stop_size; /* size of graph for stoping criteria */
    int co_stop_level; /* number of level for stoping criteria */
    double co_stop_ratio; /* if coarsening ratio higher than this, coarsening will stop */

    int co_match; /* describe the matching we use for the coarsening */
    int co_dir; /* describe wich edges we look at */
    int co_match_level; /* describe if the matching are based on toplevel or bottomlevel values */
    int co_norder; /* describe the node traversal order for the matching */
    int co_norder_reverse; /* do we reverse norder? */
    int co_eorder; /* describe the edge traversal order for the matching */
    int co_match_onedegreefirst; /* do we traverse one degree nodes first */
    int co_match_isolatedmatch; /* do we match the isolated nodes together in coarsening */
    int inipart; /* describe the initial partitioning */
    int inipart_nrun; /* number of initial partitioning */
    int conpar; /* describe if we do a constraint partitioning */
    int conpar_fix; /* describe if we fix constraint partitioning */
    int anchored; /* describe if we anchor sources and targets after constraint partitioning */
    int refinement; /* describe the refinement */
    int obj; /* describe the objective to minimize */
    int obj_commvol; /* describe how to compute communication volume */
    int kernighan; /* describe the variant of the Kernighan algorithm */
    int live_cut; /* describe the cuting strategy for rMLGPlive */
    int live_cut_merge; /* describe if we merge parts after recursive cut */
    // int live_comp; /*do we keep the output data in the cache when computing live set*/
    int cache_size; /* For rMLGPcache application: part size and cache size could be different*/
    int live_traversal; /*what traversal do we use to compute the liveset*/
    idxType maxliveset;
    int use_binary_input; /*  if not present, generate a binary file using graph data.
                        otherwise, use binary version of the input graph to read */
    //double part_ub; /* upper bound for the partitions */
    //double part_lb; /* lower bound for the partitions */

    double* ub; /* upper bound for the partitions */
    double* lb; /* lower bound for the partitions */

    int nbProc; /*Number of available processor in runtime*/
    double ccr; /*When randomizing weights what is the comm/comp ratio?*/
    int ref_step; /* Number of reffinement steps */
    int undir_alg; /* undirected initial partitioning algorithm */
    int debug; /* Debug value */
    int toggle; /* Toggle for debugging */
    int print; /* Print value */
    int seed; /* Random seed */
    double ratio; /* Imbalance ratio */
    int runs; /* Number of runs */
    char out_file[PATH_MAX];
    int write_parts;
    int ignore_livesize;
} MLGP_option;




MLGP_option* copyOpt(const MLGP_option opt, int new_nbPart);
void reallocUBLB(MLGP_option* opt);
void printOptions(const MLGP_option* opt, char lineDelim);
void print_opt(MLGP_option* opt);
void free_opt(MLGP_option* opt);


void printMLGPusage(char *exeName);
int initMLGPoptions(MLGP_option *opt, int nbPart);
int processMLGPargs(int argc, char **argv, MLGP_option* opt);

#endif
