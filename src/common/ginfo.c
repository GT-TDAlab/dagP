#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <getopt.h>


#include "dgraph.h"
#include "utils.h"
#include "debug.h"
#include "dgraphReader.h"

/*********************************************************************************/
void printUsage_ginfo(char *exeName)
{
    printf("Usage: %s <input-file-name> [options...]\n", exeName);
    printf("Option:\n");


    printf("\t--tabular int \t(default = 0)\n");
    printf("\t\tTabular format outputp\n");

}


int main(int argc, char *argv[])
{
    
    idxType numSourcesAndTargets;
    if (argc < 2){
    printUsage_ginfo(argv[0]);
    u_errexit("%s: There is a problem with the arguments.\n", argv[0]);
    }

    char fname[PATH_MAX];
    strncpy(fname, argv[1], PATH_MAX);


    int tabular = 0;

    static struct option long_options[] =
    {
    {"tabular",  required_argument, 0, 't'},
    };

    const char *delims = "t";
    char o;
    int option_index = 0;
    while ((o = (char) getopt_long(argc, argv, delims, long_options, &option_index)) != -1) {
        switch (o) {
            case 't':
                tabular = 1;
                break;
            default:
                printUsage_ginfo(argv[0]);
        u_errexit("Option %s not recognized.", o);
                return 1;
        }
    }

    dgraph G;

    readDGraph(&G, fname, 1);
    int maxindegree, minindegree, maxoutdegree, minoutdegree;
    double aveindegree, aveoutdegree;
    dgraph_info(&G, &maxindegree, &minindegree, &aveindegree, &maxoutdegree, &minoutdegree, &aveoutdegree);
    idxType* sourceOutput = (idxType*) umalloc( sizeof(idxType)*(G.nVrtx+1),"sourceOutput");
    int sCount = sourcesList(&G, sourceOutput);
    int tCount = outputList(&G, sourceOutput);
    computeToplevels(&G, sourceOutput);
    idxType i, critpath = 0;
    for (i = 1; i < G.nVrtx; i++)
        critpath = critpath < sourceOutput[i] ? sourceOutput[i] : critpath;

    numSourcesAndTargets = 0;
    for (i = 1; i <= G.nVrtx; i++)
    {
        if(G.inEnd[i] < G.inStart[i]  && G.outEnd[i] < G.outStart[i] )
            numSourcesAndTargets += 1;
    }
    free(sourceOutput);
    if (tabular)
        printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%d\t%d\t%d\t%d\n", fname, G.nVrtx, G.nEdge, minindegree, maxindegree, minoutdegree, maxoutdegree, aveoutdegree,sCount,tCount,numSourcesAndTargets,critpath);
    else {
        printf("Nb node: %d\n", G.nVrtx);
        printf("Nb Edges: %d\n", G.nEdge);
        printf("Min in degree: %d\n", minindegree);
        printf("Max in degree: %d\n", maxindegree);
        printf("Min out degree: %d\n", minoutdegree);
        printf("Max out degree: %d\n", maxoutdegree);
        printf("Average in/out degree: %f\n", aveoutdegree);
        printf("Source Node Count: %d\n", sCount);
        printf("Target Node Count: %d\n", tCount);
        printf("Isoloated vertices %d\n", numSourcesAndTargets);
        printf("Critical path: %d\n\n", critpath);

    }
    return 0;
}

