#include "debug.h"

void createDGraph_example(dgraph* G)
{
    allocateDGraphData(G, 9, 14, 3);
    G->inStart[1] = 0;
    G-> inStart[2] = 0;
    G->inStart[3] = 0;
    G-> inStart[4] = 0;
    G->inStart[5] = 2;
    G-> inStart[6] = 5;
    G->inStart[7] = 7;
    G-> inStart[8] = 9;
    G->inStart[9] = 12;
    G->inStart[10] = 14;
    int i;
    for (i=1; i<=9; i++)
        G->inEnd[i] = G->inStart[i+1]-1;
    G->inEnd[10] = 13;
    G->in[0] = 1;
    G->in[1] = 2;
    G->in[2] = 1;
    G->in[3] = 2;
    G->in[4] = 3;
    G->in[5] = 2;
    G->in[6] = 3;
    G->in[7] = 4;
    G->in[8] = 5;
    G->in[9] = 4;
    G->in[10] = 5;
    G->in[11] = 6;
    G->in[12] = 5;
    G->in[13] = 6;
    for (i=1; i<=9; i++)
	G->vw[i] = 1;
    for (i=0; i<14; i++){
	G->ecIn[i] = 1;
	G->ecOut[i] = 1;
    }

    fillOutFromIn(G);
}


/********************************************************************************/
void print_dgraph(dgraph G)
{
    printf("\n### Begin print dgraph ###\n");
    printf("Nunber of nodes: %d\n", (int) G.nVrtx);
    printf("Number of edges: %d\n", (int) G.nEdge);
    int i,j;
    int frmt = G.frmt;

    if (frmt != 2 && frmt != 3){
	printf("Print edges from out tables\n");
	for (i=1; i<=G.nVrtx; i++)
	    for (j=G.outStart[i]; j<=G.outEnd[i]; j++)
		printf("%d -> %d\n", (int) i, (int) G.out[j]);
    }

    if (frmt == 1 || frmt == 3){
	printf("Print node weights\n");
	for (i=1; i<=G.nVrtx; i++)
	    printf("Weight %i = %d\n", i, G.vw[i]);
    }

    if (frmt == 2 || frmt == 3){
	printf("Print edge weights from out tables\n");
	for (i=1; i<=G.nVrtx; i++)
	    for (j=G.outStart[i]; j<=G.outEnd[i]; j++)
		printf("%d -> %d (ceil(w) = %d)\n", (int) i, (int) G.out[j], (int)ceil(G.ecOut[j]));
    }
    printf("### End print dgraph ###\n\n");
}


void dgraph_to_dot(dgraph* G, idxType* part, char* file_name)
{
    FILE *file;
    char* col[30];
    col[0] = "red";
    col[1] = "blue";
    col[2] = "green";
    col[3] = "purple";
    col[4] = "orange";
    col[5] = "brown";
    col[6] = "chocolate";
    col[7] = "crimson";
    col[8] = "salmon";
    col[9] = "grey";
    col[10] = "gold";
    col[11] = "indigo";

    int i,j;
    file = (FILE *) fopen(file_name, "w");
    if (file == NULL)
        u_errexit("dgraph_to_dot: cannot open file %s\n", file_name);

    fprintf(file, "digraph G {\n");
    for (i=1; i<=G->nVrtx; i++){
	   fprintf(file, "%d [style=filled, fillcolor=", i-1);
       if (part)
	       fprintf(file, "%s", col[part[i] % 12]);
       else
           fprintf(file, "red");
        fprintf(file, ", weight=%d];\n", G->vw[i]);
    }
    for (i=1; i<=G->nVrtx; i++){
    	for (j=G->outStart[i]; j<=G->outEnd[i]; j++){
    	    idxType outnode = G->out[j];
    	    fprintf(file, "%d->%d [weight=%d];\n", i-1, outnode-1, (int) ceil(G->ecOut[j]));
    	}
    }
    fprintf(file, "}\n");
    fclose(file);
}

void dgraph_to_dot_with_schedule(dgraph* G, idxType* part, idxType* sched, char* file_name)
{
    //graph node index are 1-based
    FILE *file;
    char* col[30];
    col[0] = "red";
    col[1] = "blue";
    col[2] = "green";
    col[3] = "purple";
    col[4] = "orange";
    col[5] = "brown";
    col[6] = "chocolate";
    col[7] = "crimson";
    col[8] = "salmon";
    col[9] = "grey";
    col[10] = "gold";
    col[11] = "indigo";

    int i,j;
    idxType* order = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    for (i=0; i<G->nVrtx; i++)
	order[sched[i]] = i;

    file = fopen(file_name, "w");
    if (file == NULL)
        u_errexit("dgraph_to_dot_with_schedule: cannot open file %s\n", file_name);

    
    fprintf(file, "digraph G {\n");
    for (i=1; i<=G->nVrtx; i++){
	fprintf(file, "%d [style=filled, fillcolor=", order[i]);
	fprintf(file, "%s", col[part[i]%12]);
	fprintf(file, ", weight=%d];\n", G->vw[i]);
    }
    for (i=1; i<=G->nVrtx; i++){
	for (j=G->outStart[i]; j<=G->outEnd[i]; j++){
	    idxType outnode = G->out[j];
	    fprintf(file, "%d->%d [weight=%d];\n", i-1, outnode-1, (int) ceil(G->ecOut[j]));
	}
    }
    fprintf(file, "}\n");
    free(order);
    fclose(file);
}

void dgraph_to_dot_with_matching(dgraph* G, idxType* leader, char* file_name)
{
    //graph node index are 1-based
    FILE *file;
    char* col[30];
    col[0] = "red";
    col[1] = "blue";
    col[2] = "green";
    col[3] = "purple";
    col[4] = "orange";
    col[5] = "brown";
    col[6] = "chocolate";
    col[7] = "crimson";
    col[8] = "salmon";
    col[9] = "grey";
    col[10] = "gold";
    col[11] = "indigo";

    int i,j;
    int leaderisnull = (leader == NULL);
    if (leaderisnull)
        leader = (idxType*) calloc(G->nVrtx+1, sizeof(idxType));
    idxType* leaderpart = (idxType*) calloc(G->nVrtx+1, sizeof(idxType));
    idxType nbleaderpart = 0;
    file = fopen(file_name, "w");
    if (file == NULL)
        u_errexit("dgraph_to_dot_with_matching: cannot open file %s\n", file_name);
    
    
    fprintf(file, "digraph G {\n");
    
    for (i=1; i<=G->nVrtx; i++){
        if (leaderpart[leader[i]] == 0)
            leaderpart[leader[i]] = nbleaderpart++;
        //fprintf(file, "%d [style=filled, fillcolor=", i-1);
        fprintf(file, "%d [style=filled, fillcolor=", i);
        fprintf(file, "%s", col[leaderpart[leader[i]] % 12]);
        fprintf(file, ", weight=%d];\n", G->vw[i]);
    }
    for (i=1; i<=G->nVrtx; i++){
        for (j=G->outStart[i]; j<=G->outEnd[i]; j++){
            idxType outnode = G->out[j];
            //fprintf(file, "%d->%d [weight=%d];\n", i-1, outnode-1, G->ecOut[j]);
            fprintf(file, "%d->%d [weight=%d];\n", (int) i, (int) outnode, (int) G->ecOut[j]);
        }
    }
    if (leaderisnull){
        free(leader);
        leader = NULL;
    }
    fprintf(file, "}\n");
    fclose(file);
    free(leaderpart);
}


void printPartMatrix(ecType** partmatrix, idxType nbpart){
    idxType i,j;
    idxType cut = 0, tot = 0;
    printf("Partmatrix:\n");
    for (i=0; i<nbpart; i++){
        for (j=0; j<nbpart; j++){
            printf("%.0f ", (double)partmatrix[i][j]);
            tot += partmatrix[i][j];
            if (i!=j)
                cut += partmatrix[i][j];
        }
        printf("\n");
    }
    printf("Tot edge: %d\n", tot);
    printf("Edge cut: %d\n", cut);
}

void print_graph(dgraph* G, idxType* part, char* path, int arg1, int arg2)
{
    char file_name[200];
    sprintf(file_name, "%s", path);
    char tmp[100];

    sprintf(tmp, "_%d_%d.dot", arg1, arg2);
    strcat(file_name, tmp);
    dgraph_to_dot(G, part, file_name);
}

