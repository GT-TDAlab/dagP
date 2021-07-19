#include <stdio.h>
#include <stdlib.h>

#include <limits.h>
#include <getopt.h>
#include <sys/stat.h>

#include "utils.h"
#include "dgraph.h"
#include "option.h"


/* for unified naming of functions */
char *iniAlgName[IP_SIZE] = {
        "ALL",
        "KERDFS",
        "KERBFS",
        "KERDFSBFS",
        "GGG",
        "FORCED",
        "Undir+Fix",
        "BFS",
        "COCKSHAKE",
        "DIA_GGG",
        "RECBISS",
        "GGG_TWO",
        "GGG_IN",
        "CONPAR",
        "UndirRand + fix",
        "UndirGGG + fix"
	};


MLGP_option* copyOpt(const MLGP_option opt,int new_nbPart)
{
    /* Creates a new MLGP_option object and copies all the variables to the copy */
    MLGP_option *newOpt = (MLGP_option*) umalloc(sizeof(MLGP_option), "newOpt");
    newOpt->nbPart = new_nbPart;
    newOpt->co_stop_size = opt.co_stop_size;
    newOpt->co_stop_level = opt.co_stop_level;
    newOpt->co_stop_ratio = opt.co_stop_ratio;
    newOpt->co_match = opt.co_match;
    newOpt->co_dir = opt.co_dir;
    newOpt->co_match_level = opt.co_match_level;
    newOpt->co_norder = opt.co_norder;
    newOpt->co_norder_reverse = opt.co_norder_reverse;
    newOpt->co_eorder = opt.co_eorder;
    newOpt->co_match_onedegreefirst = opt.co_match_onedegreefirst;
    newOpt->co_match_isolatedmatch = opt.co_match_isolatedmatch;
    newOpt->inipart = opt.inipart;
    newOpt->inipart_nrun = opt.inipart_nrun;
    newOpt->refinement = opt.refinement;
    newOpt->conpar = opt.conpar;
    newOpt->anchored = opt.anchored;
    newOpt->conpar_fix = opt.conpar_fix;
    newOpt->obj = opt.obj;
    newOpt->nbProc = opt.nbProc;
    newOpt->kernighan = opt.kernighan;
    newOpt->undir_alg = opt.undir_alg;

    newOpt->ub = (double*) umalloc((newOpt->nbPart+1) * sizeof(double), "newOpt->ub");
    newOpt->lb = (double*) umalloc((newOpt->nbPart+1) * sizeof(double), "newOpt->lb");
    int i;
    for (i=0; i<newOpt->nbPart; ++i) {
        newOpt->ub[i] = opt.ub[i];
        newOpt->lb[i] = opt.lb[i];
    }
    newOpt->ref_step = opt.ref_step;
    newOpt->ccr = opt.ccr;
    // newOpt->live_comp = opt.live_comp;
    newOpt->live_traversal = opt.live_traversal;
    newOpt->cache_size = opt.cache_size;
    newOpt->maxliveset = opt.maxliveset;
    newOpt->debug = opt.debug;
    newOpt->toggle = opt.toggle;
    newOpt->print = opt.print;
    newOpt->seed = opt.seed;
    newOpt->ratio = opt.ratio;
    newOpt->runs = opt.runs;
    newOpt->ignore_livesize = opt.ignore_livesize;
    strcpy(newOpt->out_file,opt.out_file);
    strcpy(newOpt->file_name,opt.file_name);
    if (opt.print > 0 || opt.debug > 0)
        strcpy(newOpt->file_name_dot,opt.file_name_dot);

    return newOpt;
}

void reallocUBLB(MLGP_option* opt)
{
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
    return;
}

void printOptions(const MLGP_option* opt, char lineDelim)
{
    printf("Options:\n", lineDelim);
    printf("file_name %s%c", opt->file_name, lineDelim);
    printf("nbPart %d%c", opt->nbPart, lineDelim);

    printf("co_stop_size %d%c", opt->co_stop_size, lineDelim);
    printf("co_stop_level %d%c", opt->co_stop_level, lineDelim);
    printf("co_stop_ratio %lf%c", opt->co_stop_ratio, lineDelim);

    printf("co_match %d%c", opt->co_match, lineDelim);
    printf("co_dir %d%c", opt->co_dir, lineDelim);
    printf("co_match_level %d%c", opt->co_match_level, lineDelim);
    printf("co_norder %d%c", opt->co_norder, lineDelim);
    printf("co_norder_reverse %d%c", opt->co_norder_reverse, lineDelim);
    printf("co_eorder %d%c", opt->co_eorder, lineDelim);
    printf("co_match_onedegreefirst %d%c", opt->co_match_onedegreefirst, lineDelim);
    printf("co_match_isolatedmatch %d%c", opt->co_match_isolatedmatch, lineDelim);
    printf("inipart %d%c", opt->inipart, lineDelim);
    printf("inipart_nrun %d%c", opt->inipart_nrun, lineDelim);
    printf("conpar %d%c", opt->conpar, lineDelim);
    printf("conpar_fix %d%c", opt->conpar_fix, lineDelim);
    printf("anchored %d%c", opt->anchored, lineDelim);
    printf("refinement %d%c", opt->refinement, lineDelim);
    printf("obj %d%c", opt->obj, lineDelim);
    // printf("kernighan %d%c", opt->kernighan, lineDelim);
    printf("live_cut %d%c", opt->live_cut, lineDelim);
    printf("live_cut_merge %d%c", opt->live_cut_merge, lineDelim);
    printf("cache_size %d%c", opt->cache_size, lineDelim);
    printf("live_traversal %d%c", opt->live_traversal, lineDelim);
    printf("maxliveset %d%c", opt->maxliveset, lineDelim);
    printf("use_binary_input %d%c", opt->use_binary_input, lineDelim);

    printf("nbProc %d%c", opt->nbProc, lineDelim);
    printf("ccr %lf%c", opt->ccr, lineDelim);
    printf("ref_step %d%c", opt->ref_step, lineDelim);
    printf("undir_alg %d%c", opt->undir_alg, lineDelim);
    printf("debug %d%c", opt->debug, lineDelim);
    printf("print %d%c", opt->print, lineDelim);
    printf("seed %d%c", opt->seed, lineDelim);
    printf("ratio %lf%c", opt->ratio, lineDelim);
    printf("runs %d%c", opt->runs, lineDelim);
    printf("out_file %s%c", opt->out_file, lineDelim);
    printf("write_parts %d%c", opt->write_parts, lineDelim);
    printf("ignore_livesize %d%c", opt->ignore_livesize, lineDelim);
    printf("%c\n\n", lineDelim);
    return;
}

void print_opt(MLGP_option* opt)
{
    idxType i;
    printf("### opt.nbPart = %d\n", opt->nbPart);
    printf("### opt.ub = ");
    for (i=0; i < opt->nbPart; i++)
            printf("%f ", opt->ub[i]);
    printf("\n");
    printf("### opt.lb = ");
    for (i=0; i < opt->nbPart; i++)
            printf("%f ", opt->lb[i]);
    printf("\n");
    printf("### opt.ratio = %f\n", opt->ratio);
    return;
}

void free_opt(MLGP_option* opt)
{
    free(opt->ub);
    free(opt->lb);
}

void printMLGPusage(char *exeName)
{
    printf("Usage: %s <input-file-name> <#parts> [options...]\n", exeName);
    printf("Options:\n");

    printf("\nPartition-related Parameters\n");
    printf("============================\n\n");
    //seed
    printf("\t--seed int \t(default = 0 [= initialize with time])\n");
    printf("\t\tSeed for the random generator\n");

    //ratio
    printf("\t--ratio float \t(default = 1.03)\n");
    printf("\t\tMaximum imbalance ratio\n");

    //debug
    printf("\t--debug int \t(default = 0)\n");
    printf("\t\tDebug value\n");

    //toggle
    printf("\t--toggle int \t(default = 0)\n");
    printf("\t\tToggle for generating plots for debugging\n");

    //print
    printf("\t--print int \t(default = 0)\n");
    printf("\t\tValue to set output details\n");
    printf("\t\t= 0 : None\n");
    printf("\t\t= 1 : Low\n");
    printf("\t\t= 2 : Medium\n");
    printf("\t\t= 3 : High\n");

    //runs
    printf("\t--runs int \t(default = 1)\n");
    printf("\t\tNumber of runs for the graph\n");

    //part_ub
    printf("\t-part_ub double \t(default = (tot_weight/nbpart) * 1.03)\n");
    printf("\t\tUpper bound for the partition\n");

    //part_lb
    printf("\t-part_lb double \t(default = (tot_weight/nbpart) / 1.03)\n");
    printf("\t\tLower bound for the partition\n");

    //use_binary_input
    printf("\t--use_binary_input (-x) int \t(default = 1)\n");
    printf("\t\tCreate/Read binary version of the input graph.\n");

    //write_parts
    printf("\t--write_parts (-5) int \t(default = 0)\n");
    printf("\t\tCreate a file containing partition number for each node. Zero for False, 1 for True.\n");
    printf("\t\tThe created file is in the same folder as the input file, named: <input>.partsfile.part_<part>.seed_<seed>.txt\n");


    printf("\nPartition-related Parameters: Coarsening/Matching\n");
    printf("=================================================\n\n");
    //co_stop_size
    printf("\t--co_stop_size int \t(default = %d * #parts)\n", CO_STOP_SIZE_MULT);
    printf("\t--co_stop_level int \t(default = %d)\n", CO_STOP_LEVEL_DEF);
    printf("\t--co_stop_ratio double \t(default = %.3lf)\n", CO_STOP_RATIO_DEF);

    //co_match
    printf("\t--co_match int \t(default = %d)\n", CO_ALG_DEF);
    printf("\t\tSet the clustering technic for the coarsening step\n");
    printf("\t\t= 0 : Pairwise matching with level difference 1\n");
    printf("\t\t= 1 : Aggregative clustering with level difference 1 (CoTop)\n");
    printf("\t\t= 2 : Aggregative clustering with level difference 1 with cycle detection (CoCyc)\n");
    printf("\t\t= 3 : Hybrid aggregative clustering with level difference 1 with cycle detection by default except for large degree nodes (CoHyb)\n");

    //co_match_level
    printf("\t--co_match_level int \t(default = %d)\n", CO_LEVEL_DEF);
    printf("\t\tDetermine if the matchings are based on toplevel or bottom level values\n");
    printf("\t\t= 0 : Top-level values\n");
    printf("\t\t= 1 : Bottom-level values\n");
    printf("\t\t= 2 : Alternates the two\n");

    //co_dir
    printf("\t--co_dir int \t(default = %d)\n", CO_DIR_DEF);
    printf("\t\tWich edges do we look at when coarsening?\n");
    printf("\t\t= 0 : We only look at in-edges\n");
    printf("\t\t= 1 : We only look at out-edges\n");
    printf("\t\t= 2 : We look at in-edges and out-edges\n");

    //co_norder
    printf("\t--co_norder int \t(default = %d)\n", CO_NORD_DEF);
    printf("\t\tSet the node traversal order for the matching phase\n");
    printf("\t\t= 0 : Random traversal\n");
    printf("\t\t= 1 : Random Depth First Topological order\n");
    printf("\t\t= 2 : Random Breadth First Topological order\n");
    printf("\t\t= 3 : Random Topological order\n");
    printf("\t\t= 4 : Depth First Topological order\n");
    printf("\t\t= 5 : Breadth First Topological order\n");
    printf("\t\t= 6 : Random Depth First traversal\n");
    printf("\t\t= 7 : Random Breadth First traversal\n");
    printf("\t\t= 8 : Depth First traversal\n");
    printf("\t\t= 9 : Breadth First traversal\n");

    //co_norder_reverse
    printf("\t--co_norder_reverse int \t(default = 0)\n");
    printf("\t\tDo we reverse the node order previously defined in the coarsening phase\n");
    printf("\t\t= 0 : No\n");
    printf("\t\t= 1 : Yes\n");

    //co_eorder
    printf("\t--co_eorder int \t(default = %d)\n", CO_EORD_DEF);
    printf("\t\tSet the edge traversal order for the matching phase\n");
    printf("\t\t= 0 : Natural construction order\n");
    printf("\t\t= 1 : Edge weight priority\n");
    printf("\t\t= 2 : Edge weight / node weight priority\n");
    printf("\t\t= 3 : Sum of edge weight coming from a cluster\n");
    printf("\t\t= 4 : Sum of edge weight coming from a cluster / weight of the cluster\n");

    //co_match_onedegreefirst
    printf("\t--co_match_onedegreefirst int \t(default = 1)\n");
    printf("\t\tDo we traverse the nodes with degree 1 first in the matching?\n");
    printf("\t\t= 0 : No\n");
    printf("\t\t= 1 : Yes\n");

    //co_match_isolatedmatch
    printf("\t--co_match_isolatedmatch int \t(default = 1)\n");
    printf("\t\tDo we match isolated nodes two by two in coarsening?\n");
    printf("\t\t= 0 : No\n");
    printf("\t\t= 1 : Yes\n");

    //conpar
    printf("\t--conpar int \t(default = %d)\n", CONPAR_DEF);
    printf("\t\tDo we run a Constraint-Partitioning to guide the coarsening phase?\n");
    printf("\t\t= 0 : No\n");
    printf("\t\t= 1 : Yes\n");

    //conpar_fix
    printf("\t--conpar_fix int \t(default = %d)\n", CONPAR_FIX_DEF);
    printf("\t\tDetermine if we fix the partitioning in the conpar\n");
    printf("\t\t= 0 : We use UndirPartitioning for conpar\n");
    printf("\t\t= 1 : We use UndirBisection + fix for conpar\n");

    //anchored
    printf("\t--anchored int \t(default = %d)\n", ANCHORED_DEF);
    printf("\t\tDetermine if we anchor the sources and the targets after running conpar\n");
    printf("\t\t= 0 : No\n");
    printf("\t\t= 1 : Yes\n");


    printf("\nPartition-related Parameters: Initial Partitioning\n");
    printf("==================================================\n\n");
    //co_init
    printf("\t--inipart int \t(default = %d (works when --conpar 1))\n", IP_CONPAR);
    printf("\t\tDetermine the initial partitioning\n");
    printf("\t\t= 0 : try all partitioning and pick the best one\n");
    printf("\t\t= 1 : Kernighan with random DFS\n");
    printf("\t\t= 2 : Kernighan with random BFS\n");
    printf("\t\t= 3 : Kernighan with random mix BFS-DFS\n");
    printf("\t\t= 4 : Greedy graph growing approach\n");
    printf("\t\t= 5 : Forced balance initial partitioning\n");
    printf("\t\t= 6 : Undirected Partitioning with fixing acyclicity and forced balance (follow with --undir_alg parameter)\n");
    printf("\t\t= 7 : BFS Growing - Fan Out\n");
    printf("\t\t= 8 : Cocktail Shake Growing - Fan In Fan Out\n");
    printf("\t\t= 9 : Greedy graph growing approach starting from diagonale\n");
    printf("\t\t= 10 : Recursive bisectiont\n");
    printf("\t\t= 11 : Greedy graph growing approach with two phases\n");
    printf("\t\t= 12 : Greedy graph growing approach looking only at in edges\n");
    printf("\t\t= 13 : Initial partitioning computed in Constraint-Partitioning (works only with --conpar 1)\n");
    printf("\t\t= 14 : Undirected Random Partitioning with fixing acyclicity and forced balance\n");
    printf("\t\t= 15 : Undirected GGG Partitioning with fixing acyclicity and forced balance\n");

    //undir_alg
    printf("\t--undir_alg int \t(default = %d)\n", UNDIR_METIS);
    printf("\t\t= 0 : Metis Undirected partitioning\n");
    printf("\t\t= 1 : Scotch Undirected partitioning\n");

    //co_nbinit
    printf("\t--inipart_nrun int \t(default = 5)\n");
    printf("\t\tNumber of runs for the initial partitioning\n");


    printf("\nPartition-related Parameters: Refinement\n");
    printf("========================================\n\n");

    //refinement
    printf("\t--refinement int \t(default = %d)\n", REF_DEF);
    printf("\t\tDetermine the refinement in uncoarsening step\n");
    printf("\t\t= 0 : no refinement\n");
    printf("\t\t= 1 : bFM: bisection FM\n");
    printf("\t\t= 2 : bFM_MAXW: bisection FM with vertex selected from maximum weighted part\n");
    printf("\t\t= 3 : KL: Kernighan-Lin (bisection with swap)\n");
    printf("\t\t= 4 : KL_bFM_MAXW: one pass Kernighan-Lin + one pass bFM_MAXW\n");
    printf("\t\t= 5 : kFM_PO: the refinement toward the greatest partition in post order\n");
    printf("\t\t= 6 : kFM_TB: the refinement toward the best partition if it creats no cycle\n");
    printf("\t\t= 7 : kFM_POTB: performes combination of kFM_PO and then kFM_TB\n");
    printf("\t\t= 8 : kFM_RANDV: refinement randomized selection of boundary vertices to move\n");



    printf("\nApplication-Related (Cache/Live-Run/Scheduling) Parameters\n");
    printf("==========================================================\n\n");

    //ccr
    printf("\t--ccr double \t(default = 0)\n");
    printf("\t\tWhen randomizing the weights, what is the comm/comp ratio?\n");

    //live_traversal
    printf("\t--live_traversal int \t(default = 0)\n");
    printf("\t\tWhat traversal do we use to compute the live set?\n");
    printf("\t\t= 0 : natural ordering (Warning natural ordering has to be a topologic order)\n");
    printf("\t\t= 1 : DFS topological order\n");
    printf("\t\t= 2 : BFS topological order\n");
    printf("\t\t= 3 : mix DFS-BFS topological order\n");
    printf("\t\t= 4 : random topological order\n");

    //obj
    printf("\t--obj int \t(default = 0)\n");
    printf("\t\tDetermine the objective we are trying to minimize\n");
    printf("\t\t= 0 : minimize edge-cut\n");
    printf("\t\t= 1 : minimize communication volume\n");

    //
    printf("\t--ignore_livesize (-0) int \t(default = 0)\n");

    //nbproc
    printf("\t--nbproc int \t(default = 2)\n");
    printf("\t\tNumber of processors available in runtime computation.\n");
}

int initMLGPoptions(MLGP_option *opt, int nbPart)
{
    opt->nbPart = nbPart;
    // opt->maxliveset = nbPart;
    opt->seed = 0;
    opt->co_stop_size = CO_STOP_SIZE_MULT*opt->nbPart;
    opt->co_stop_level = CO_STOP_LEVEL_DEF;
    opt->co_stop_ratio = CO_STOP_RATIO_DEF;

    opt->co_match = CO_ALG_DEF;
    opt->co_dir = CO_DIR_DEF;
    opt->co_match_level = CO_LEVEL_DEF;
    opt->co_norder = CO_NORD_DEF;
    opt->co_norder_reverse = 0;
    opt->co_eorder = CO_EORD_DEF;
    opt->co_match_onedegreefirst = 0;
    opt->co_match_isolatedmatch = 1;
    opt->inipart = IP_DEF;
    opt->undir_alg = UNDIR_METIS;
    opt->inipart_nrun = 5;
    opt->refinement = REF_DEF;
    opt->ccr = 0;
    opt->kernighan = -1;
    opt->cache_size = -1;
    opt->live_traversal = 0;
    opt->conpar = CONPAR_DEF;
    opt->anchored = ANCHORED_DEF;
    opt->conpar_fix = CONPAR_FIX_DEF;
    opt->write_parts = 0;
    opt->ignore_livesize = 0;
    opt->ub = (double*) umalloc((opt->nbPart+1)*sizeof(double), "opt->ub");
    opt->lb = (double*) umalloc((opt->nbPart+1)*sizeof(double), "opt->lb");
    int i=0;
    for (i=0; i<opt->nbPart; ++i) {
        opt->ub[i] = -1;
        opt->lb[i] = -1;
    }
    opt->ref_step = 10;
    opt->nbProc = 2;

    opt->obj = CO_OBJ_EC;
    opt->ratio = 1.03;
    opt->runs = 1;
    opt->debug = 0;
    opt->toggle = 0;
    opt->print = 0;
    opt->use_binary_input = 1;
    opt->live_cut = 0;
    opt->live_cut_merge = 1;
    return 0;
}

int processMLGPargs(int argc, char **argv, MLGP_option* opt)
{
    int inipart_is_default = 1;

    strncpy(opt->file_name, argv[1], PATH_MAX);
    strncpy(opt->out_file, opt->file_name, PATH_MAX);
    strcat(opt->out_file, "_MLGP.out");

    // opt->nbPart = atoi(argv[2]);

    static struct option long_options[] =
    {
        {"seed",  required_argument, 0, 's'},
        {"co_stop_size",  required_argument, 0, 'z'},
        {"co_stop_level",  required_argument, 0, 'p'},
        {"co_stop_ratio",  required_argument, 0, 't'},
        {"co_match",  required_argument, 0, 'm'},
        {"co_dir",  required_argument, 0, 'h'},
        {"co_match_level",  required_argument, 0, '2'},
        {"co_norder",  required_argument, 0, 'o'},
        {"co_norder_reverse",  required_argument, 0, 'v'},
        {"co_eorder",  required_argument, 0, 'g'},
        {"co_match_onedegreefirst",  required_argument, 0, '1'},
        {"co_match_isolatedmatch",  required_argument, 0, '6'},
        {"inipart",  required_argument, 0, 'i'},
        {"inipart_nrun",  required_argument, 0, 'b'},
        {"refinement",  required_argument, 0, 'r'},
        {"conpar",  required_argument, 0, 'a'},
        {"conpar_fix",  required_argument, 0, 'u'},
        {"anchored",  required_argument, 0, 'e'},
        {"obj",  required_argument, 0, 'j'},
        {"nbproc",  required_argument, 0, '9'},
        {"part_ub",  required_argument, 0, '3'},
        {"part_lb",    required_argument, 0, 'l'},
        {"ratio",    required_argument, 0, '4'},
        {"runs",    required_argument, 0, 'n'},
        {"debug",    required_argument, 0, 'd'},
        {"toggle",    required_argument, 0, 'k'},
        {"print",    required_argument, 0, 'q'},
        {"ref_step",    required_argument, 0, 'f'},
        {"ccr",    required_argument, 0, '*'},
        {"use_binary_input", required_argument, 0, 'x'},
        {"cache_size",    required_argument, 0, '7'},
        {"live_traversal",    required_argument, 0, '8'},
        {"undir_alg",	required_argument, 0, 'y'},
        {"write_parts", required_argument, 0, '5'},
        {"ignore_livesize", required_argument, 0, '0'},
        {NULL, 0, NULL, 0}
    };

    const char *delims = "s:t:z:p:m:e:o:v:g:1:i:b:r:j:u:l:a:n:d:k:q:f:x:c:w:h:y:4:3:2:5:6:7:8:9:0:*";

    char o;
    int option_index = 0;
    while ((o = (char) getopt_long(argc, argv, delims, long_options, &option_index)) != -1) {
        switch (o) {
            case 's':
                opt->seed = atoi(optarg);
                break;
            case 'z':
                opt->co_stop_size = atoi(optarg);
                break;
            case 'p':
                opt->co_stop_level = atoi(optarg);
                break;
            case 't':
                opt->co_stop_ratio = atof(optarg);
                break;
            case 'm':
                opt->co_match = atoi(optarg);
                break;
            case 'h':
                opt->co_dir = atoi(optarg);
                break;
            case '2':
                opt->co_match_level = atoi(optarg);
                break;
            case 'o':
                opt->co_norder = atoi(optarg);
                break;
            case 'v':
                opt->co_norder_reverse = atoi(optarg);
                break;
            case 'g':
                opt->co_eorder = atoi(optarg);
                break;
            case '1':
                opt->co_match_onedegreefirst = atoi(optarg);
                break;
            case '6':
                opt->co_match_isolatedmatch = atoi(optarg);
                break;
            case 'i':
                inipart_is_default = 0;
                opt->inipart = atoi(optarg);
                break;
            case 'b':
                opt->inipart_nrun = atoi(optarg);
                break;
            case 'r':
                opt->refinement = atoi(optarg);
                break;
            case 'a':
                opt->conpar = atoi(optarg);
                break;
            case 'u':
                opt->conpar_fix = atoi(optarg);
                break;
            case 'e':
                opt->anchored = atoi(optarg);
                break;
            case 'j':
                opt->obj = atoi(optarg);
                break;
            case '*':
                opt->ccr = (double) atoi(optarg);
                break;
            case '9':
                opt->nbProc = atoi(optarg);
                break;
            case '3':
                opt->ub[0] = atof(optarg);
                break;
            case 'l':
                opt->lb[0] = atof(optarg);
                break;
            case '4':
                opt->ratio = atof(optarg);
                break;
            case 'd':
                opt->debug = atoi(optarg);
                break;
            case 'k':
                opt->toggle = atoi(optarg);
                break;
            case 'q':
                opt->print = atoi(optarg);
                break;
            case 'n':
                opt->runs = atoi(optarg);
                break;
            case 'f':
                opt->ref_step = atoi(optarg);
                break;
            case '8':
                opt->live_traversal = atoi(optarg);
                break;
            case '7':
                opt->cache_size = atoi(optarg);
                break;
            case 'x':
                opt->use_binary_input = atoi(optarg);
                break;
            case '5':
                opt->write_parts = atoi(optarg);
                break;
            case '0':
                opt->ignore_livesize = atoi(optarg);
                break;
            case 'y':
        		opt->undir_alg = atoi(optarg);
        		break;
            default:
                printMLGPusage(argv[0]);
                u_errexit("Option %c not recognized.", optopt);
                return 1;
        }
    }

    for(int i=1; i<opt->nbPart; ++i) {
        opt->ub[i]=opt->ub[0];
        opt->lb[i]=opt->lb[0];
    }

    if (opt->print > 0 || opt->debug > 0) {
        char suf[PATH_MAX];
        char fname[PATH_MAX];
        char* split;
        char fname_sav[PATH_MAX];

        sprintf(fname_sav, "%s", opt->file_name);

        split = strtok(fname_sav, "/");
        while(split != NULL) {
            sprintf(fname, "%s", split);
            split = strtok(NULL, "/");
        }
        sprintf(opt->file_name_dot, "out/");
        mkdir(opt->file_name_dot, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH);
        strcat(opt->file_name_dot, fname);
        sprintf(suf, "__%d_%d", opt->nbPart, opt->seed);
        strcat(opt->file_name_dot, suf);
        mkdir(opt->file_name_dot, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH);
        sprintf(suf, "/graph");
        strcat(opt->file_name_dot, suf);
    }

    if ((opt->conpar == 1)&&(opt->conpar_fix == 1)&&(inipart_is_default == 1))
        opt->inipart = IP_CONPAR;

#if !defined(dagP_SCOTCH) && !defined(dagP_METIS)
        printf("Warning: ");
        printf("Code is not linked with MeTiS or Scotch, setting initial partitioning to inipart=IP_GGG_TWO and conpar=0\n");
        opt->inipart = IP_GGG_TWO;
        opt->conpar = 0;
#endif

    if (((opt->conpar == 0)||((opt->conpar == 1)&&(opt->conpar_fix == 0)))&&(opt->inipart == IP_CONPAR))
        u_errexit("Initial partitioning %d only works with --conpar 1\n", IP_CONPAR);
    return 0;
}

