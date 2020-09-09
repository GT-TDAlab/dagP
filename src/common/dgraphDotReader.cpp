
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <fcntl.h>

#include "dgraphReader.h"
#include <map>
#include <string>
#include <vector>



#define DOT_LINE_MAX 1000

/* note that this ignores whole like if it has any of these chars */
static inline int dotCommentLine(char *line)
{
    return (line[0] == '#') || strstr(line, "//") || strstr(line, "/*") || strstr(line, "*/");
}

int readDGraphDot(dgraph *G, char* file_name)
{
    char delims[] = " \t->[,=];\n";
    int graphStatus=0; /* 0 not started, 1 saw '{' */

    char line[DOT_LINE_MAX];
    idxType nVrtx = 0, nEdge = 0;
    idxType i, j;
    idxType* nbneighbors=NULL;
    idxType** inneighbors=NULL;
    ecType** ecInNeighbors = NULL;
    std::map <std::string, int> nodeNameMap;
    std::vector <int> vweights;
    vweights.reserve(1000);
    int readForNodes = 2;
    int ignoreLines = 0;
    while(readForNodes){
        --readForNodes;
        FILE* file = fopen(file_name, "r");

        if (file == NULL)
            u_errexit("readDGraphDot: Cannot open file %s\n", file_name);
        while (!feof(file)) {
            do {
                fgets(line, DOT_LINE_MAX, file);
            } while (!feof(file) && dotCommentLine(line));
            if (feof(file))
                break;

            if (strlen(line)>=LINE_MAX-2)
                u_errexit("readDGraphDot: line is too long! \n'%'\n", line);

            char *pos = strstr(line, "{");
            if (pos != NULL) {
                graphStatus = 1;
                continue;
            }
	    if (graphStatus ==0) {
		continue;
	    }
            pos = strstr(line, "}");
            if (pos != NULL) { // end of graph no need to continue
                break;
            }
            if(ignoreLines){
                pos = strstr(line,"];");
                if (pos != NULL){
                    ignoreLines = 0;
                    continue;
                }
            }
            if (ignoreLines){
                // printf("line ignored %s\n", line);
                continue;
            }

            if (graphStatus<1) {
                u_errexit("readDGRaphDot: graph open curly not seen line='%s' in file %s\n", line,file_name);
            }
            /* it is an edge */
            pos = strstr(line, "--");
            if (pos != NULL) {
                u_errexit("readDGraphDot: graph is undirected.");
                return -1;
            }
            pos = strstr(line, "->");
            if (pos != NULL) {
                if (readForNodes)
                    continue;
                // printf("%s is edge line\n",line);
                nEdge++;
                if (nEdge == 1) {
                    nbneighbors = (idxType*) calloc(nVrtx+1, sizeof(idxType));
                    inneighbors = (idxType**) malloc((nVrtx+1) * sizeof(idxType*));
                    ecInNeighbors = (ecType**) malloc((nVrtx+1) * sizeof(ecType*));
                    for (i = 1; i <= nVrtx; i++){
                        inneighbors[i] = (idxType*) malloc(100 * sizeof(idxType));
                        ecInNeighbors[i] = (ecType*) malloc(100 * sizeof(ecType));
                    }
                }
                char* token = strtok(line, delims);
                std::string two(token);
                if (nodeNameMap.count(two)==0){
                    u_errexit("readDGraphDot: node name not seen in the input file before an edge %s",token);
                }
                int v2 = nodeNameMap[two];
                token = strtok(NULL, delims);
                std::string one(token);
                if (nodeNameMap.count(one)==0){
                    u_errexit("readDGraphDot: node name not seen in the input file before an edge %s",token);
                }
                int v1 = nodeNameMap[one];
                // int v1 = atoi(token);
                int w=1;
                token = strtok(NULL, delims);
                if (token!=NULL) {
                    if (!strcmp("style", token) || !strcmp("fillcolor", token)) /* ignore syle and fill color */
                        token = strtok(NULL, delims);
                    else if (!strcmp("weight", token) || !strcmp("Weight",token)) {
                        token = strtok(NULL, delims);
                        w = atoi(token);
                    } else
                        u_errexit("readDGraphDot: while reading edge unknown token '%s' line='%s'", token, line);
                }
                if ((nbneighbors[v1+1] + 1) % 100 == 0){
                    inneighbors[v1+1] = (idxType*) realloc (inneighbors[v1+1], (nbneighbors[v1+1] + 101) * sizeof(idxType));
                    if (inneighbors == NULL)
                        printf("realloc failed for inneighbors\n");
                    ecInNeighbors[v1+1] = (ecType*) realloc (ecInNeighbors[v1+1], (nbneighbors[v1+1] + 101) * sizeof(ecType));
                    if (ecInNeighbors == NULL)
                        printf("realloc failed for ecInNeighbors\n");
                }
                inneighbors[v1+1][nbneighbors[v1+1]] = v2+1;
                ecInNeighbors[v1+1][nbneighbors[v1+1]] = w;
                nbneighbors[v1+1]++;
            } else { /* it is a vertex */
                if(readForNodes == 0)
                    continue;
                // printf("%s is vertex line\n",line);
                if (strstr(line,"graph")!=NULL){
                    ignoreLines = 1;
                    if (strstr(line, "];")!=NULL){
                        ignoreLines = 0;
                    }
                    continue;
                }
                char* token = strtok(line, delims); /* vertex name; we assume they are 0-based */
                std::string one(token);
                if (nodeNameMap.count(one)==0){
                    nodeNameMap[one] = nVrtx;
                }
                else{
                    u_errexit("we added this node before! %s", token);
                }
                nVrtx++;
                // printf("%d\n", nodeNameMap.size());
                token = strtok(NULL, delims);
                int w=1;
                if (token!=NULL) {
                    if (!strcmp("style", token) || !strcmp("fillcolor", token) || !strcmp("color", token)) /* ignore syle and fill color */
                        token = strtok(NULL, delims);
                    else if (!strcmp("weight", token) || !strcmp("Weight",token)) {
                        token = strtok(NULL, delims);
                        w = atoi(token);
                    } else
                        u_errexit("readDGraphDot: while reading vertex unknown token '%s' line='%s'", token, line);
                }
                vweights.push_back(w);
            }
        }
        fclose(file);

    }
    char nodeMappingsFile[PATH_MAX];
    sprintf(nodeMappingsFile, "%s.nodemappings", file_name);
    FILE *f = fopen(nodeMappingsFile,"w");
    for (auto& s: nodeNameMap)
        fprintf(f,"%s %d\n",s.first.c_str(), s.second+1);
    fclose(f);
    G->frmt = 3;
    G->nVrtx = nVrtx;
    G->nEdge = nEdge;
    G->totvw = nVrtx;
    G->maxindegree = 0;
    G->hollow  = (int * ) calloc(nVrtx + 2, sizeof(int));
    G->inStart  = (idxType * ) calloc(nVrtx + 2, sizeof(idxType));
    G->inEnd= (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->in = (idxType * ) malloc(sizeof(idxType) * nEdge);
    idxType idx = 0, degree;
    G->ecIn = (ecType *) malloc(sizeof(ecType) * nEdge);
    G->ecOut = (ecType *) malloc(sizeof(ecType) * nEdge);

    for (i=1; i<=nVrtx; i++) {
        G->inStart[i] = idx;
        G->inEnd[i-1] = idx-1;
        if (i>1) {
            degree = G->inEnd[i-1] - G->inStart[i-1] + 1;
            G->maxindegree = G->maxindegree < degree ? degree : G->maxindegree;
        }
        for (j=0; j< nbneighbors[i]; j++) {
            G->ecIn[idx] = ecInNeighbors[i][j];
            G->in[idx++] = inneighbors[i][j];
        }
    }
    G->inStart[0] = 0;
    G->inEnd[0] = -1;
    G->inEnd[nVrtx] = idx-1;
    G->inEnd[nVrtx+1] = nEdge;
    G->inStart[nVrtx+1] = nEdge;
    degree = G->inEnd[nVrtx] - G->inStart[nVrtx] + 1;
    G->maxindegree = G->maxindegree < degree ? degree : G->maxindegree;
    if (nEdge <= 0)
        u_errexit("allocateDGraphData: empty edge set\n");

    if (nVrtx <=-1)
        u_errexit("allocateDGraphData: empty vertex set\n");

    G->outStart  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->outEnd  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 2));
    G->out = (idxType * ) malloc(sizeof(idxType) * nEdge);

    G->vw = (vwType *) malloc(sizeof(vwType) * (nVrtx+1));
    for (i=1; i<=nVrtx; i++)
        G->vw[i] = vweights[i-1];

    fillOutFromIn(G);

    G->sources  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 1));
    G->targets  = (idxType * ) malloc(sizeof(idxType) * (nVrtx + 1));
    G->nbsources = 0;
    G->nbtargets = 0;
    for (i=1; i<=nVrtx; i++){
        if (G->inEnd[i] < G->inStart[i])
            G->sources[G->nbsources++] = i;
        if (G->outEnd[i] < G->outStart[i])
            G->targets[G->nbtargets++] = i;
    }

    for (i = 1; i <= nVrtx; i++)
        free(inneighbors[i]);
    free(nbneighbors);
    free(inneighbors);
    return 0;
}
