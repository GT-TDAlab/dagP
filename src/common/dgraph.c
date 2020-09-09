#include "dgraph.h"
#include "utils.h"

int computePartmatrix(dgraph *G, idxType *part, ecType** partmatrix, idxType node, idxType movetopart)
{
    idxType j;
    int flag = 0;
    for (j=G->outStart[node]; j<=G->outEnd[node]; j++) {
        idxType out = G->out[j];

        if (partmatrix[part[node]][part[out]] - G->ecOut[j]== 0 )
            flag=1;
        partmatrix[part[node]][part[out]] -= G->ecOut[j];

        if (partmatrix[movetopart][part[out]] == 0)
            flag=1;
        partmatrix[movetopart][part[out]] += G->ecOut[j];
    }
    for (j=G->inStart[node]; j<=G->inEnd[node]; j++) {
        idxType in = G->in[j];
        if (partmatrix[part[in]][part[node]] - G->ecIn[j] == 0)
            flag = 1;
        partmatrix[part[in]][part[node]] -= G->ecIn[j];
        if (partmatrix[part[in]][movetopart] == 0)
            flag =1;
        partmatrix[part[in]][movetopart] += G->ecIn[j];
    }
    return flag;
}

void buildPartmatrix(dgraph *G, idxType * part, idxType nbpart, ecType** partmatrix)
{
    idxType i,j;
    for (i=0; i< nbpart; i++)
        for (j=0; j<nbpart; j++)
            partmatrix[i][j] = 0;
    for (i=1; i<= G->nVrtx; i++){
        for (j=G->outStart[i]; j<=G->outEnd[i]; j++){
            idxType node = G->out[j];
            partmatrix[part[i]][part[node]] += G->ecOut[j];
        }
    }
}

int thereIsCycle(idxType nbpart, ecType** partmatrix)
{
    idxType* ready = (idxType*) malloc(nbpart * sizeof(idxType));
    idxType* traversed = (idxType*) calloc(nbpart, sizeof(idxType));
    idxType nbready = 0, i, j, nbtraversed = 0;

    for (i = 0; i < nbpart; i++){
        int nbin = 0;
        for (j = 0; j < nbpart; j++)
            if ((partmatrix[j][i] > 0)&&(j!=i))
                nbin++;
        if (nbin == 0)
            ready[nbready++] = i;
    }
    while (nbready > 0){
        idxType node = ready[--nbready];
        traversed[node] = 1;
        nbtraversed++;
        for (i = 0; i <nbpart; i++){
            if (partmatrix[node][i] == 0)
                continue;
            if (traversed[i] == 1)
                continue;
            int is_ready = 1;
            for (j = 0; j < nbpart; j++){
                if ((partmatrix[j][i] == 0)||(j==i))
                    continue;
                if (traversed[j] == 0)
                    is_ready = 0;
            }
            if (is_ready == 1)
                ready[nbready++] = i;
        }
    }
    free(ready);
    free(traversed);
    return (nbtraversed != nbpart);
}

int partSizeChecker(idxType* partsize, double ub_pw[], int nbpart, vwType maxvw){
    //as long as at least one of them is bigger than bound...
    int i;
    for ( i=0; i<nbpart; ++i )
        if (partsize[i] > ub_pw[i]+maxvw)
            return 1;
    return 0;
}


void allocateDGraphData(dgraph *G, idxType nVrtx, idxType nEdge, int frmt)
{
    /*assumes *G is allocated*/
    G->frmt = frmt;
    G->nVrtx = nVrtx;
    G->nEdge = nEdge;

    G->vw = NULL;
    G->ecIn = G->ecOut = NULL;

    if (nEdge <= 0)
        u_errexit("allocateDGraphData: empty edge set\n");

    if (nVrtx <=-1)
        u_errexit("allocateDGraphData: empty vertex set\n");

    G->inStart  = (idxType * ) umalloc(sizeof(idxType) * (nVrtx + 2), "G->inStart");
    G->inEnd= (idxType * ) umalloc(sizeof(idxType) * (nVrtx + 2), "G->inEnd");
    G->in = (idxType * ) umalloc(sizeof(idxType) * nEdge, "G->in");
    G->hollow = (int *) calloc(nVrtx+2, sizeof(int));
    G->outStart  = (idxType * ) umalloc(sizeof(idxType) * (nVrtx + 2), "G->outStart");
    G->outEnd  = (idxType * ) umalloc(sizeof(idxType) * (nVrtx + 2), "G->outEnd");
    G->out = (idxType * ) umalloc(sizeof(idxType) * nEdge, "G->out");
    G->sources = (idxType *) umalloc(sizeof(idxType) * (nVrtx+2), "G->sources");
    G->targets = (idxType *) umalloc(sizeof(idxType) * (nVrtx+2), "G->targets");

    if(G->hollow == NULL)
        u_errexit("Failed allocation of G->hollow");

    if (frmt == 1 || frmt == 3)
        G->vw = (vwType *) umalloc(sizeof(vwType) * (nVrtx+1), "G->vw");

    if (frmt == 2 || frmt == 3)
    {
        G->ecIn = (ecType *) umalloc(sizeof(ecType) * nEdge, "G->ecIn");
        G->ecOut = (ecType *) umalloc(sizeof(ecType) * nEdge, "G->ecOut");
    }
    G->maxVW = G->maxEC = -1;
}

void freeDGraphData(dgraph *G)
{
    free(G->inStart);
    free(G->inEnd);
    free(G->in);
    free(G->outStart);
    free(G->outEnd);
    free(G->out);
    free(G->hollow);
    free(G->sources);
    free(G->targets);

    if (G->frmt == 1 || G->frmt == 3)
    {
        if (!G->vw)
            u_errexit("free a graph with frmt %d, but vw is not allocated", G->frmt);
        else
            free(G->vw);
    }

    if (G->frmt == 2 || G->frmt == 3)
    {
        if (!G->ecIn) {
            u_errexit("free a graph with frmt %d, but ecIn is not allocated", G->frmt);
        }
        else
            free(G->ecIn);

        if (!G->ecOut)
            u_errexit("free a graph with frmt %d, but ecOut is not allocated", G->frmt);
        else
            free(G->ecOut);
    }

    G->frmt = -1;
    G->nVrtx = 0;
    G->nEdge = 0;
    G->maxVW = G->maxEC = -1;

    G->out = G->outEnd = G->outStart = G->in = G->inEnd = G->inStart = NULL;

    G->vw = NULL;
    G->ecIn = G->ecOut = NULL;
}

void dgraph_info(dgraph* G, int* maxindegree, int* minindegree, double* aveindegree, int* maxoutdegree, int* minoutdegree, double* aveoutdegree)
{
    idxType i;
    int indegree, outdegree;
    *aveindegree = 0.0;
    *aveoutdegree = 0.0;
    *maxindegree = 0;
    *minindegree = 0;
    *maxoutdegree = 0;
    *minoutdegree = 0;
    for (i=1; i<=G->nVrtx; i++) {
	indegree = G->inEnd[i]-G->inStart[i]+1;
	outdegree = G->outEnd[i]-G->outStart[i]+1;
	if (indegree > *maxindegree)
	    *maxindegree = indegree;
	if (indegree < *minindegree)
	    *minindegree = indegree;
	*aveindegree += indegree;
	if (outdegree > *maxoutdegree)
	    *maxoutdegree = outdegree;
	if (outdegree < *minoutdegree)
	    *minoutdegree = outdegree;
	*aveoutdegree += outdegree;
    }
    *aveindegree = *aveindegree / G->nVrtx;
    *aveoutdegree = *aveoutdegree / G->nVrtx;
}

void set_dgraph_info(dgraph* G)
{
    idxType i,j;

    double totvw = 0.0;
    double totec = 0.0;
    idxType maxindegree = 0;
    idxType maxoutdegree = 0;
    vwType maxVW = 0;
    ecType maxEC = 0;

    for (i=1; i<=G->nVrtx; i++) {
        if (maxindegree < G->inEnd[i] - G->inStart[i]+1)
            maxindegree = G->inEnd[i] - G->inStart[i]+1;

        totvw += G->vw[i];
        if (maxVW < G->vw[i])
            maxVW = G->vw[i];

        for (j=G->inStart[i]; j<=G->inEnd[i]; j++) {
            totec += G->ecIn[j];
            if (maxEC < G->ecIn[j])
                maxEC = G->ecIn[j];
        }
        if (maxoutdegree < G->outEnd[i] - G->outStart[i]+1)
            maxoutdegree = G->outEnd[i] - G->outStart[i]+1;

    }
    G->totvw = totvw;
    G->totec = totec;
    G->maxindegree = maxindegree;
    G->maxoutdegree = maxoutdegree;
    G->maxVW = maxVW;
    G->maxEC = maxEC;
}


typedef struct{
    idxType ngh;
    ecType cst;
} sortData;

int dgqsortfnct(const void *i1, const void *i2)
{
    idxType ii1=((sortData*)i1)->ngh, ii2 = ((sortData*)i2)->ngh;
    if ( ii1 < ii2) return -1;
    else if (ii1 > ii2) return 1;
    else return 0;
}

void sortNeighborLists(dgraph *G)
{
    /*sorts in and out in the increasing order of neighbor indices*/
    idxType i, j, at, nVrtx = G->nVrtx;
    idxType *out=G->out, *outStart = G->outStart, *outEnd = G->outEnd;
    idxType *in=G->in, *inStart = G->inStart, *inEnd = G->inEnd;
    ecType *ecIn=G->ecIn, *ecOut = G->ecOut;

    sortData *sd = (sortData*) umalloc(nVrtx * sizeof(sortData), "sd");

    /*sort in*/
    for (j=1; j <= nVrtx; j++)
    {
        at = 0;
        for (i = inStart[j]; i <= inEnd[j]; i++)
        {
            sd[at].ngh = in[i];
            if (G->frmt == 2 || G->frmt == 3)
                sd[at].cst = ecIn[i];
            at ++;
        }
        qsort(sd, inEnd[j]-inStart[j]+1, sizeof(sortData), dgqsortfnct);
        at = 0;
        for (i = inStart[j]; i <= inEnd[j]; i++)
        {
            in[i] =  sd[at].ngh ;
            if (G->frmt == 2 || G->frmt == 3)
                ecIn[i] = sd[at].cst ;
            at ++;
        }
    }

    /*sort out*/
    for (j=1; j <= nVrtx; j++)
    {
        at = 0;
        for (i = outStart[j]; i <= outEnd[j]; i++)
        {
            sd[at].ngh = out[i];
            if (G->frmt == 2 || G->frmt == 3)
                sd[at].cst = ecOut[i];
            at ++;
        }
        qsort(sd, outEnd[j]-outStart[j]+1, sizeof(sortData), dgqsortfnct);
        at = 0;
        for (i = outStart[j]; i <= outEnd[j]; i++)
        {
            out[i] =  sd[at].ngh ;
            if (G->frmt == 2 || G->frmt == 3)
                ecOut[i] = sd[at].cst ;
            at ++;
        }
    }

    free(sd);
}

void fillOutFromIn(dgraph *G)
{
    idxType i, at, j, nVrtx = G->nVrtx;
    idxType *out=G->out, *outStart = G->outStart, *outEnd=G->outEnd;
    idxType *in=G->in, *inStart = G->inStart, *inEnd = G->inEnd;
    ecType *ecIn=G->ecIn, *ecOut = G->ecOut;

    for (i = 0; i <= G->nVrtx+1; i++)
        outStart[i] = 0;

    /*count*/
    for (j = 1; j<= nVrtx; j++) {
        for (i = inStart[j]; i <= inEnd[j]; i++) {
		  outStart[G->in[i]]++;
	    }
    }
    /*prefix sum*/
    for (j = 1; j<= nVrtx+1; j++) {
        outStart[j] += outStart[j-1];
    }

    /*write*/
    for (j = 1; j<= nVrtx; j++) {
        for (i = inStart[j]; i <= inEnd[j]; i++) {
            at =  --outStart[in[i]];
            out[at] = j;
            if (G->frmt & DG_FRMT_EC)
                ecOut[at] = ecIn[i];
        }
    }
    G->maxoutdegree = 0;
    for (j = 1; j<= nVrtx; j++) {
        outEnd[j] = outStart[j+1]-1;
	    idxType degree = outEnd[j] - outStart[j] + 1;
	    G->maxoutdegree = G->maxoutdegree < degree ? degree : G->maxoutdegree;
    }

    if (inStart[1] != 0)
        u_errexit("fillOutFromIn: the first index not right");

    if (inStart[nVrtx+1] != G->nEdge)
        u_errexit("fillOutFromIn: the last index not right");
}

void fillInFromOut(dgraph *G)
{
    /*Fill G->in, G->inStart, G->inEnd, G->ecIn and G->maxindegree
     * based on G->out, G->outStart, G->outEnd, G->ecOut and G->maxoutdegree
     */
    idxType i, at, j, nVrtx = G->nVrtx;
    idxType *out=G->out, *outStart = G->outStart, *outEnd=G->outEnd;
    idxType *in=G->in, *inStart = G->inStart, *inEnd = G->inEnd;
    ecType *ecIn=G->ecIn, *ecOut = G->ecOut;
    for (i = 0; i <= G->nVrtx+1; i++)
        inStart[i] = 0;

    /*count*/
    for (j = 1; j<= nVrtx; j++) {
        for (i = outStart[j]; i <= outEnd[j]; i++) {
            inStart[G->out[i]]++;
        }
    }
    /*prefix sum*/
    for (j = 1; j<= nVrtx+1; j++) {
        inStart[j] += inStart[j-1];
    }

    /*write*/
    for (j = 1; j<= nVrtx; j++) {
        for (i = outStart[j]; i <= outEnd[j]; i++) {
            at =  --inStart[out[i]];
            in[at] = j;
            if (G->frmt & DG_FRMT_EC)
                ecIn[at] = ecOut[i];
        }
    }
    G->maxindegree = 0;

    for (j = 1; j<= nVrtx; j++) {
        inEnd[j] = inStart[j+1]-1;
        idxType degree = inEnd[j] - inStart[j] + 1;
        G->maxindegree = G->maxindegree < degree ? degree : G->maxindegree;
    }

    if (outStart[1] != 0)
        u_errexit("fillInFromOut: the first index not right");

    if (outStart[nVrtx+1] != G->nEdge)
        u_errexit("fillInFromOut: the last index not right");
}


void setVertexWeights(dgraph *G, vwType *vw)
{
    /*assuming n weights in vw*/
    idxType i, nVrtx = G->nVrtx;
    G->totvw = 0.0;
    G->maxVW = -1;

    if (G->frmt == 2 || G->frmt == 0)
    {
        if (G->vw == NULL)
            G->vw = (vwType * ) malloc(sizeof(vwType) * (nVrtx+1));
        G->frmt +=1;/*if ec, then ec+vw, if no w no cost, then vw*/
    }
    for (i = 1; i <= nVrtx; i++)
    {
        G->vw[i] = vw[i-1];
        G->totvw +=vw[i-1];
        G->maxVW = G->maxVW < vw[i-1] ? vw[i-1] : G->maxVW;
    }

}

int checkAcyclicity(dgraph *G, idxType *part, idxType nbpart)
{
    idxType i, j, k, ip;
    idxType** outpart = (idxType**) malloc(sizeof(idxType*)*nbpart);
    idxType* nbout = (idxType*) malloc(sizeof(idxType)*nbpart);
    idxType** inpart = (idxType**) malloc(sizeof(idxType*)*nbpart);
    idxType* nbin = (idxType*) malloc(sizeof(idxType)*nbpart);
    for (i=0; i<nbpart; i++) {
        outpart[i] = (idxType*) malloc(sizeof(idxType)*nbpart);
        nbout[i] = 0;
        inpart[i] = (idxType*) malloc(sizeof(idxType)*nbpart);
        nbin[i] = 0;
    }
    int isAcyclic = 1;

    for (i=1; i<=G->nVrtx; i++) {
        for (j=G->outStart[i]; j<=G->outEnd[i]; j++) {
            idxType outnode = G->out[j];
            if (part[i] == part[outnode])
                continue;

            int is_new = 1;
            for (k=0; k<nbout[part[i]]; k++)
                if (outpart[part[i]][k] == part[outnode]) {
                    is_new = 0;
                    break;
                }

            if (is_new == 1) {
                outpart[part[i]][nbout[part[i]]] = part[outnode];
                nbout[part[i]]++;
            }

            is_new = 1;
            for (k=0; k<nbin[part[outnode]]; k++)
                if (inpart[part[outnode]][k] == part[i]) {
                    is_new = 0;
                    break;
            }

            if (is_new == 1) {
                inpart[part[outnode]][nbin[part[outnode]]] = part[i];
                nbin[part[outnode]]++;
            }
        }
    }

    idxType* ready = (idxType*) malloc(sizeof(idxType)*(nbpart));
    idxType* nbinleft = (idxType*) malloc(sizeof(idxType)*(nbpart));
    int nbready = 0;
    for (i=0; i<nbpart; i++) {
        nbinleft[i] = nbin[i];
        if (nbin[i] == 0)
            ready[nbready++] = i;
    }

    int to = 0;
    while (nbready > 0) {
        idxType pno = ready[nbready-1];
        nbready--;
        to++;
        for (ip = 0; ip < nbout[pno]; ip++) {
            idxType succ = outpart[pno][ip];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0) {
                ready[nbready++] = succ;
            }
        }
    }
    for (i=0; i<nbpart; i++) {
        free(outpart[i]);
        free(inpart[i]);
    }
    free(nbinleft);
    free(ready);
    free(nbin);
    free(nbout);
    free(inpart);
    free(outpart);
    if (to != nbpart) {
        isAcyclic = 0;
    }
    return isAcyclic;
}

void removeMarkedEdges(dgraph* G, int* marked, idxType nbmarked){
    /* Modify G to remove marked edges
     */
    idxType i,j, outidx = 0;
    G->maxEC = G->totec = 0;

    for (i = 1; i <= G->nVrtx; i++){
        idxType formerstart = G->outStart[i];
        G->outStart[i] = outidx;
        G->outEnd[i-1] = outidx-1;
        for (j = formerstart; j <= G->outEnd[i]; j++){
            if (marked[j] == 1) {
                continue;
            }
            G->out[outidx] = G->out[j];
            G->ecOut[outidx] = G->ecOut[j];
            G->maxEC = G->maxEC < G->ecOut[j] ? G->ecOut[j] : G->maxEC;
            G->totec += G->ecOut[j];
            outidx++;
        }
    }
    G->outEnd[G->nVrtx] = outidx-1;
    G->outStart[G->nVrtx+1] = G->nEdge - nbmarked;
    G->outEnd[G->nVrtx+1] = G->nEdge - nbmarked - 1;
    G->maxoutdegree = 0;
    for (i = 1; i <= G->nVrtx; i++) {
        idxType degree = G->outEnd[i] - G->outStart[i] + 1;
        G->maxoutdegree = G->maxoutdegree < degree ? degree : G->maxoutdegree;
    }
    G->ecOut = (ecType*) realloc(G->ecOut, (G->nEdge - nbmarked + 1)*sizeof(ecType));
    G->out = (idxType*) realloc(G->out, (G->nEdge - nbmarked + 1)*sizeof(idxType));
    G->in = (idxType*) realloc(G->in, (G->nEdge - nbmarked + 1)*sizeof(idxType));
    G->ecIn = (ecType*) realloc(G->ecIn, (G->nEdge - nbmarked + 1)*sizeof(ecType));
    G->nEdge = G->nEdge - nbmarked;

    fillInFromOut(G);
}

idxType findPath(dgraph* G, idxType source, int* istarget, idxType* path, idxType* toplevel, idxType maxtoplevel){
    /* Find a path in G from the source to one of the target and return its length
     * Source is not considered as one of the targets in this function
     */
    if (toplevel[source] >= maxtoplevel)
        return 0;

    for (idxType i = G->outStart[source]; i <= G->outEnd[source]; i++){
        idxType outnode = G->out[i];
        if (istarget[outnode]){
            path[0] = i;
            return 1;
        }
        idxType length = findPath(G, outnode, istarget, path+1, toplevel, maxtoplevel);
        if (length > 0){
            path[0] = i;
            return length + 1;
        }
    }
    return 0;
}

void transitiveReduction(dgraph* G){
    /*Modify G to make it transitive irreductible
     * (It uses toplevel to optimize the traversal)
     */
    idxType* toplevel = (idxType*) malloc(sizeof(idxType) * (G->nVrtx + 1));
    computeToplevels(G, toplevel);

    int* toRemove = (int*) calloc(G->nEdge+1, sizeof(int));
    int* istarget = (int*) calloc(G->nVrtx+1, sizeof(int));
    idxType* path = (idxType*) malloc(sizeof(idxType)*(G->nEdge+1));
    idxType i,j, nbmarked = 0;

    for (i = 1; i <= G->nVrtx; i++){
        idxType maxtoplevel = -1, length = 0;
        for (j = G->outStart[i]; j<= G->outEnd[i]; j++){
            istarget[G->out[j]] = 1;
            maxtoplevel = maxtoplevel < toplevel[G->out[j]] ? toplevel[G->out[j]] : maxtoplevel;
        }
        for (j = G->outStart[i]; j<= G->outEnd[i]; j++) {
            length = findPath(G, G->out[j], istarget, path, toplevel, maxtoplevel);
            if (length > 0)
                break;
        }
        for (j = G->outStart[i]; j<= G->outEnd[i]; j++)
            istarget[G->out[j]] = 0;
        if (length > 0 ){
            idxType edgeweigth = -1;
            for (j = G->outStart[i]; j<= G->outEnd[i]; j++)
                if (G->out[j] == G->out[path[length-1]]) {
                    toRemove[j] = 1;
                    nbmarked++;
                    edgeweigth = G->ecOut[j];
                    break;
                }
            for (j = 0; j <= length-1; j++)
                G->ecOut[j] += edgeweigth;
        }
    }

    removeMarkedEdges(G, toRemove, nbmarked);

    free(toRemove);
    free(istarget);
    free(path);
    free(toplevel);
}


idxType sourcesList(dgraph* G, idxType* sources)
{
    /*Assume that sources is already allocated and fill the tab sources with the sources nodes of G
      Return the number of sources in G*/
    idxType i, ind = 0;
    for (i=1; i<=G->nVrtx; i++)
        if (G->inStart[i] > G->inEnd[i])
            sources[ind++] = i;
    return ind;
}

idxType sourcesListPart(dgraph* G, idxType* sources, idxType *part, int part_idx)
{
    /*Assume that sources is already allocated and fill the tab sources with the sources nodes of part_idx
      Return the number of sources in G*/
    idxType i, ind = 0, is_source, j;
    for (i=1; i<=G->nVrtx; i++) {
        if (part[i] != part_idx)
            continue;
        is_source = 1;
        for (j=G->inStart[i]; j <= G->inEnd[i]; j++)
            if (part[G->in[j]] == part_idx)
                is_source = 0;
        if (is_source == 1)
            sources[ind++] = i;
    }
    return ind;
}

idxType outputList(dgraph* G, idxType* output)
{
    /*Assume that output is already allocated and fill the tab output with the output nodes of G
      Return the number of output in G*/
    idxType i, ind = 0;
    for (i=1; i<=G->nVrtx; i++)
        if (G->outStart[i] > G->outEnd[i])
            output[ind++] = i;
    return ind;
}

void oneDegreeFirst(dgraph* G, idxType* order)
{
    idxType i,j;
    u_errexit("oneDegreeFirst: this function, as implemented, should not be called at all\n");
    for (i=1; i<=G->nVrtx; i++) {
	    idxType node = order[i];
	    int nboutnodes = G->outEnd[node] - G->outStart[node] + 1;
        int nbinnodes = G->inEnd[node] - G->inStart[node] + 1;
	    if ((nboutnodes == 1)&&(nbinnodes == 1)) {
	        for (j=i; j>=2; j--)
		        order[j] = order[j-1];
	        order[1] = node;
	    }
    }
}

idxType farthestNode(dgraph* G, idxType startnode)
{
    idxType to = 1, succ, node, i;
    idxType* ready = (idxType*) malloc(sizeof(idxType) * (G->nVrtx+1));
    int* touched = (int*) calloc(G->nVrtx+1, sizeof(int));
    ready[0] = startnode;
    touched[startnode] = 1;
    idxType endready = 0;
    idxType beginready = 0;

    while (endready >= beginready) {
        if (to > G->nVrtx)
            u_errexit("Proof to = %d\n", to);
        node = ready[beginready];
        beginready++;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            succ = G->out[i];
            if (touched[succ] == 0) {
                ready[++endready] = succ;
                touched[succ] = 1;
            }
        }
        for (i = G->inStart[node]; i <= G->inEnd[node]; i++) {
            succ = G->in[i];
            if (touched[succ] == 0) {
                ready[++endready] = succ;
                touched[succ] = 1;
            }
        }
    }
    return node;
}

void topsort(dgraph* G, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a topological order of nodes in G
    with DFS*/
    idxType to = 1;
    idxType* ready = (idxType*) malloc(sizeof(idxType) * (G->nVrtx+1));
    idxType i;
    idxType nbready = G->nbsources;
    for (i = 0; i < nbready; i++)
        ready[i] = G->sources[i];
    sourcesList(G, ready);
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    while (nbready > 0) {
        if (to > G->nVrtx)
            u_errexit("Proof topsort to = %d\n", to);
        idxType node = ready[nbready-1];
        nbready--;
        toporder[to++] = node;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
            else if (nbinleft[succ]<0)
                u_errexit("topsort: negative indegree\n");
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("topsort : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(ready);
}

void computeToplevels(dgraph* G, idxType* toplevels)
{
    /*We assume that toplevels is already allocated*/
    int i,j;
    idxType* toporder = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    topsort(G, toporder);
    for (i=1; i<=G->nVrtx; i++) {
        idxType node = toporder[i];
        idxType tl = 0;

        for (j=G->inStart[node]; j<= G->inEnd[node]; j++)
            tl = tl > toplevels[G->in[j]]+1 ? tl : toplevels[G->in[j]]+1;
        toplevels[node] = tl;
    }
    free(toporder);
}
void computeWeightedToplevels(dgraph* G, ecType* toplevels)
{
    /*We assume that toplevels is already allocated*/
    int i,j;
    idxType* toporder = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    topsort(G, toporder);
    for (i=1; i<=G->nVrtx; i++) {
        idxType node = toporder[i];
        idxType tl = 0;

        for (j=G->inStart[node]; j<= G->inEnd[node]; j++)
            // tl = tl > toplevels[G->in[j]] + G->vw[G->in[j]] + G->ecIn[j] ? tl : toplevels[G->in[j]] + G->vw[G->in[j]] + G->ecIn[j];
            if (tl < toplevels[G->in[j]] + G->vw[G->in[j]] + G->ecIn[j])
                tl = toplevels[G->in[j]] + G->vw[G->in[j]] + G->ecIn[j];
        toplevels[node] = tl;
    }
    free(toporder);
}

void computeBottomlevels(dgraph* G, idxType* bottomlevels)
{
    /*We assume that bottomlevels is already allocated*/
    int i,j;
    idxType* toporder = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    topsort(G, toporder);
    for (i=1; i<=G->nVrtx; i++) {
        idxType node = toporder[G->nVrtx-i+1];
        idxType bl = 0;

        for (j=G->outStart[node]; j<= G->outEnd[node]; j++)
            bl = bl > bottomlevels[G->out[j]]+1 ? bl : bottomlevels[G->out[j]]+1;
        bottomlevels[node] = bl;
    }
    free(toporder);
}

void computeWeightedBottomlevels(dgraph* G, idxType* bottomlevels)
{
    /*We assume that bottomlevels is already allocated*/
    int i,j;
    idxType* toporder = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    topsort(G, toporder);
    for (i=1; i<=G->nVrtx; i++) {
        idxType node = toporder[G->nVrtx-i+1];
        idxType bl = 0;

        for (j=G->outStart[node]; j<= G->outEnd[node]; j++)
            bl = bl > bottomlevels[G->out[j]]+ G->ecOut[j] ? bl : bottomlevels[G->out[j]]+ G->ecOut[j];
        bottomlevels[node] = bl + G->vw[node];
    }
    free(toporder);
}

void computeToplevelsWithTopOrder(dgraph* G, idxType* toplevels, idxType* toporder)
{
    /*We assume that toplevels is already computed and toprder allocated*/
    idxType i,j;
    for (i=1; i<=G->nVrtx; i++)
    {
        idxType node = toporder[i];
        int tl = 0;
        for (j=G->inStart[node]; j<= G->inEnd[node]; j++)
            tl = tl > toplevels[G->in[j]]+1 ? tl : toplevels[G->in[j]]+1;
        toplevels[node] = tl;
    }
}


double computeLatency(dgraph* G, idxType* part, double l1, double l2)
{
    int i,j;
    double max = 0.0;
    double* latencies = (double*) malloc(sizeof(double)*(G->nVrtx + 1));
    idxType* toporder = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    topsort(G, toporder);
    for (i=1; i<=G->nVrtx; i++)
    {
        idxType node = toporder[i];
	double lat = 1.0;
	for (j=G->inStart[node]; j<= G->inEnd[node]; j++) {
	  double currentLat = part[node] == part[G->in[j]] ? l1 : l2;
	  lat = lat > latencies[G->in[j]]+currentLat+1.0 ? lat : latencies[G->in[j]]+currentLat+1.0;
	}
	latencies[node] = lat;
	max = max > lat ? max : lat;
    }
    free(latencies);
    free(toporder);
    return max;
}

void computeDistances(dgraph* G, idxType sourceNode, idxType* dist)
{
    /*Fill the distance array with the distance toward sourceNode in G*/
    idxType i,j;
    for (i=1; i<=G->nVrtx; i++)
        dist[i] = G->nVrtx+1;
    idxType* ready = (idxType*) malloc(sizeof(idxType)*(G->nVrtx + 1));
    ready[0] = sourceNode;
    idxType nbready = 1;
    dist[sourceNode] = 0;
    while(nbready > 0){
        idxType node = ready[nbready-1];
        nbready--;
        for (j=G->inStart[node]; j<= G->inEnd[node]; j++) {
            idxType inNode = G->in[j];
            if (dist[inNode] == G->nVrtx+1){
                dist[inNode] = dist[node] + 1;
                ready[nbready++] = inNode;
            }
        }
        for (j=G->outStart[node]; j<= G->outEnd[node]; j++) {
            idxType outNode = G->out[j];
            if (dist[outNode] == G->nVrtx+1){
                dist[outNode] = dist[node] + 1;
                ready[nbready++] = outNode;
            }
        }
    }
    free(ready);
}

void connectedComponents(dgraph* G, idxType** components, idxType* sizes, idxType* nbcomp){
    /*Fill the components[i] array with the list of nodes in this connected component*/
    idxType i,j;
    int* marked = (int*) calloc(G->nVrtx + 1, sizeof(int));
    int nbmarked = 0;
    idxType* ready = (idxType*) malloc(sizeof(idxType)*(G->nVrtx + 1));
    idxType nbready = 1;
    int lastsource = 0;
    idxType source = 1;

    *nbcomp = 0;

    while (nbmarked < G->nVrtx){
        components[*nbcomp] = (idxType*) malloc((G->nVrtx - nbmarked+1)*sizeof(idxType));
        sizes[*nbcomp] = 0;
        while (marked[source] == 1)
            source++;
        marked[source] = 1;
        nbmarked++;
        ready[0] = source;
        nbready = 1;
        while(nbready > 0){
            idxType node = ready[--nbready];
            components[*nbcomp][sizes[*nbcomp]++] = node;
            for (j=G->inStart[node]; j<= G->inEnd[node]; j++) {
                idxType inNode = G->in[j];
                if (marked[inNode] == 0){
                    marked[inNode] = 1;
                    nbmarked++;
                    ready[nbready++] = inNode;
                }
            }
            for (j=G->outStart[node]; j<= G->outEnd[node]; j++) {
                idxType outNode = G->out[j];
                if (marked[outNode] == 0){
                    marked[outNode] = 1;
                    nbmarked++;
                    ready[nbready++] = outNode;
                }
            }
        }
        *nbcomp = *nbcomp + 1;
    }

    free(ready);
    free(marked);
}

void checkGraphAcyclicity(dgraph* G)
{
    /*Try to run a topolotical order. If cannot, then error.*/
    idxType to = 1;
    idxType* ready = (idxType*) malloc(sizeof(idxType) * (G->nVrtx+1));
    idxType i;
    idxType nbready = sourcesList(G, ready);
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    while (nbready > 0) {
        idxType node = ready[nbready-1];
        nbready--;
        to++;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
            else if (nbinleft[succ]<0)
                u_errexit("Graph is not acyclic\n");
        }
    }
    to = 1;
    nbready = outputList(G, ready);
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->outEnd[i] - G->outStart[i] + 1;
    while (nbready > 0) {
        idxType node = ready[nbready-1];
        nbready--;
        to++;
        for (i = G->inStart[node]; i <= G->inEnd[node]; i++) {
            idxType pred = G->in[i];
            nbinleft[pred]--;
            if (nbinleft[pred] == 0)
                ready[nbready++] = pred;
            else if (nbinleft[pred]<0)
                u_errexit("Graph is not acyclic\n");
        }
    }

    if (to != G->nVrtx+1)
        u_errexit("Graph is not acyclic to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(ready);
}

void addSingleSourceTarget(dgraph* G, idxType* flag) {
    /*Let's add single source and single target*/
    /*And we change flag*/
    idxType i,j;
    idxType flagin = 0;
    idxType flagout = 1;
    for (i = 0; i < G->nVrtx; i++)
        for (j = G->inStart[i]; j <= G->inEnd[i]; j++)
            if (flag[i] != flag[G->in[j]]){
                flagin = flag[G->in[j]];
                flagout = flag[i];
            }

    //We compute sources and issource
    idxType *sources = G->sources;
    idxType nbsources = G->nbsources;
    int *issource = (int *) calloc(G->nVrtx + 3, sizeof(int));
    for (i = 0; i < nbsources; i++)
        issource[sources[i]] = 1;
    //We compute targets and istargets
    idxType *targets = G->targets;
    idxType nbtargets = G->nbtargets;
    int *istarget = (int *) calloc(G->nVrtx + 1, sizeof(int));
    for (i = 0; i < nbtargets; i++)
        istarget[targets[i]] = 1;
    //Index for the source and target
    idxType source = G->nVrtx + 1;
    idxType target = G->nVrtx + 2;
    //Update nVrtx, nEdge, totvw and maxindegree
    G->nVrtx += 2;
    G->nEdge += nbsources + nbtargets;
    G->totvw += 2;
    G->maxindegree = G->maxindegree < nbsources ? nbsources : G->maxindegree;
    //Realloc pointers
    idxType* newinStart = (idxType *) malloc((G->nVrtx + 2)*sizeof(idxType));
    int* newhollow = (int *) calloc(G->nVrtx + 2, sizeof(int));
    idxType* newinEnd = (idxType *) malloc((G->nVrtx + 2)*sizeof(idxType));
    idxType* newin = (idxType *) malloc(G->nEdge * sizeof(idxType));
    ecType* newecIn = (ecType *) malloc(G->nEdge * sizeof(ecType));
    vwType* newvw = (vwType *) malloc((G->nVrtx+1)*sizeof(vwType));
    //Fill new vw
    for (i = 1; i <= G->nVrtx - 2; i++)
        newvw[i] = G->vw[i];
    newvw[source] = 1;
    newvw[target] = 1;
    //Fil new inStart
    newinStart[0] = 0;
    idxType sourceidx = 0;
    for (i = 1; i <= G->nVrtx - 2; i++){
        newinStart[i] = G->inStart[i] + sourceidx;
        if (issource[i] == 1)
            sourceidx++;
    }
    newinStart[source] = G->inStart[source] + sourceidx;
    newinStart[target] = G->inStart[source] + sourceidx;
    newinStart[target+1] = G->nEdge;
    //Fill new inEnd
    newinEnd[0] = -1;
    for (i = 1; i <= G->nVrtx; i++)
        newinEnd[i] = newinStart[i+1] - 1;
    //newinEnd[G->nVrtx] = G->nEdge;
    newinEnd[G->nVrtx + 1] = G->nEdge;
    //Fill new in and ecIn
    idxType oldidx = 0;
    idxType newidx = 0;
    for (i = 1; i <= G->nVrtx-2; i++){
        if (issource[i] == 1) {
            if (flag[i] == flagin)
                newecIn[newidx] = G->nEdge / 2;
            else
                newecIn[newidx] = 1;
            newin[newidx++] = source;
        }
        for (j = G->inStart[i]; j <= G->inEnd[i]; j++) {
            newecIn[newidx] = G->ecIn[j];
            newin[newidx++] = G->in[j];
        }
    }
    for (j = newinStart[target]; j <= newinEnd[target]; j++) {
        if (flag[targets[j - newinStart[target]]] == flagout)
            newecIn[j] = G->nEdge / 2;
        else
            newecIn[j] = 1;
        newin[j] = targets[j - newinStart[target]];
    }
    //Realloc outStart, outEnd and out to fill them

    idxType* newoutStart  = (idxType * ) malloc((G->nVrtx + 2)*sizeof(idxType));
    idxType* newoutEnd  = (idxType * ) malloc((G->nVrtx + 2)*sizeof(idxType));
    idxType* newout = (idxType * ) malloc(G->nEdge * sizeof(idxType));
    ecType* newecOut = (ecType * ) malloc(G->nEdge * sizeof(ecType));

    //free everything
    free(G->hollow);
    G->hollow = newhollow;
    free(G->in);
    G->in = newin;
    free(G->ecIn);
    G->ecIn = newecIn;
    free(G->inEnd);
    G->inEnd = newinEnd;
    free(G->inStart);
    G->inStart = newinStart;
    free(G->out);
    G->out = newout;
    free(G->ecOut);
    G->ecOut = newecOut;
    free(G->outEnd);
    G->outEnd = newoutEnd;
    free(G->outStart);
    G->outStart = newoutStart;
    free(G->vw);
    G->vw = newvw;

    checkDgraph(G);
    fillOutFromIn(G);

    G->sources[0] = source;
    G->targets[0] = target;
    G->nbsources = 1;
    G->nbtargets = 1;

    free(issource);
    free(istarget);

    flag = (idxType *) realloc(flag, (G->nVrtx + 1)*sizeof(idxType));
    flag[source] = flagin;
    flag[target] = flagout;
}


void copyDgraph(dgraph* Gcopy, dgraph* G) {
    idxType i;
    allocateDGraphData(G, Gcopy->nVrtx, Gcopy->nEdge, Gcopy->frmt);
    G->nVrtx = Gcopy->nVrtx;
    G->nEdge = Gcopy->nEdge;
    G->frmt = Gcopy->frmt;
    G->nbsources = Gcopy->nbsources;
    G->nbtargets = Gcopy->nbtargets;
    G->totvw = Gcopy->totvw;
    G->totec = Gcopy->totec;
    G->maxindegree = Gcopy->maxindegree;
    G->maxoutdegree = Gcopy->maxoutdegree;
    G->maxVW = Gcopy->maxVW;
    G->maxEC = Gcopy->maxEC;

    for(i = 0; i <= G->nVrtx + 1; i++)
        G->inStart[i] = Gcopy->inStart[i];

    for(i = 0; i <= G->nVrtx + 1; i++)
        G->inEnd[i] = Gcopy->inEnd[i];

    for(i = 0; i <= G->nEdge -1; i++)
        G->in[i] = Gcopy->in[i];

    for(i = 0; i <= G->nVrtx + 1; i++)
        G->hollow[i] = Gcopy->hollow[i];

    for(i = 0; i <= G->nVrtx + 1; i++)
        G->outStart[i] = Gcopy->outStart[i];

    for(i = 0; i <= G->nVrtx + 1; i++)
        G->outEnd[i] = Gcopy->outEnd[i];

    for(i = 0; i <= G->nEdge -1; i++)
        G->out[i] = Gcopy->out[i];


    if (G->frmt == 1 || G->frmt == 3)
        for(i = 0; i <= G->nVrtx; i++)
            G->vw[i] = Gcopy->vw[i];


    if (G->frmt == 2 || G->frmt == 3)
        for(i = 0; i <= G->nEdge - 1; i++) {
            G->ecIn[i] = Gcopy->ecIn[i];
            G->ecOut[i] = Gcopy->ecOut[i];
        }

    for(i = 0; i <= G->nVrtx; i++)
        G->sources[i] = Gcopy->sources[i];

    for(i = 0; i <= G->nVrtx; i++)
        G->targets[i] = Gcopy->targets[i];
}

void checkDgraph(dgraph* G){
    idxType i;
    for (i = 1; i <= G->nVrtx-1; i++) {
        if (G->inStart[i] > G->inStart[i + 1])
            u_errexit("inStart not correct");
        if (G->inEnd[i] > G->inEnd[i + 1])
            u_errexit("inEnd not correct");
    }
    for (i = 1; i <= G->nVrtx; i++) {
        if (G->inStart[i] > G->inEnd[i] + 1)
            u_errexit("inEnd and inStart not coherent (inStart[%d] = %d, inEnd[%d] = %d", i, G->inStart[i], i, G->inEnd[i]);
        if (G->inStart[i] > G->nEdge + 1)
            u_errexit("inStart not coherent");
        if (G->inEnd[i] > G->nEdge + 1)
            u_errexit("inEnd not coherent");
    }
    if (G->inStart[1] != 0)
        u_errexit("First index of inStart not right");

    if (G->inStart[G->nVrtx+1] != G->nEdge)
        u_errexit("Last index of inStart not right");
}

void analyzeDGraph(dgraph *G)
{
    idxType i, j, l, w, nlvls, nVrtx = G->nVrtx;
    idxType *toplevels, *histogram;
    ecType *ecLvls; /*we will store from the cost from the level i to level i+1 in ecLvls[i]; for i=0, u < nlvls*/

    idxType wndwLngth = 25 ;

    /*Graph's pointers*/
    idxType *outStart = G->outStart;
    idxType *outEnd = G->outEnd;
    idxType *out = G->out;
    ecType *ecOut = G->ecOut;

    toplevels = (idxType*) malloc(sizeof(idxType) * (nVrtx+1));

    computeToplevels(G, toplevels);

    nlvls = 0;
    for (i=1; i<= nVrtx; i++)
        nlvls = nlvls < toplevels[i] ? toplevels[i] : nlvls;

    nlvls = nlvls+1;

    histogram = (idxType*) calloc(nlvls, sizeof(idxType));
    ecLvls =     (ecType*) calloc(nlvls, sizeof(ecType));/*nlvls is enough for ecLvls but let us keep it parallel*/
    if (nlvls >= 1000)
        wndwLngth = (idxType) ceil(nlvls/25.0);
    else if (nlvls>=500)
        wndwLngth = 12;
    else if (nlvls>=250)
        wndwLngth = 10;
    else if (nlvls>=125)
        wndwLngth = 5;
    else if (nlvls>=65)
        wndwLngth = 2;
    else
        wndwLngth = 1;

    for (i=1; i<= nVrtx; i++) {
      histogram[toplevels[i]]+=G->vw[i];
        for (j=outStart[i]; j<= outEnd[i]; j++) {
            idxType ng = out[j];
            if (toplevels[ng]<= toplevels[i])
                u_errexit("analyzeDGraph: did not like the graph\n");
            if (((int)toplevels[ng]/wndwLngth)!=((int)toplevels[i]/wndwLngth))
		ecLvls[toplevels[i]] += ecOut[j];
        }
    }

    printf("\t\t%ld vertices in %ld levels (S=%ld, T=%ld) \n", (long) nVrtx, (long) nlvls,
           (long) histogram[0], (long)  histogram[nlvls-1]);
    for (w = 0; w < nlvls; w += wndwLngth) {
        idxType nVrtxWnd = 0;
        ecType ecWnd = 0;
        for (l = w; l < nlvls && l < w + wndwLngth; l++) {
            nVrtxWnd += histogram[l];
            ecWnd += ecLvls[l];
        }
        printf("L%ld-%ld [%ld]: ec-> [%.0f]\n", (long) w, (long) l-1, (long) nVrtxWnd, (1.0)*ecWnd);

    }
    free(ecLvls);
    free(histogram);
    free(toplevels);
}

void printNodeInfo(dgraph *G,idxType* part,idxType node)
{
    int j=0;
        printf("Node: %d (%d): [", node, part[node]);
        for (j=G->inStart[node]; j<=G->inEnd[node]; j++){
            idxType father = G->in[j];
            printf("%d (%d), ", father,part[father] );
        }
        printf("]<" );
        for (j=G->outStart[node]; j<=G->outEnd[node]; j++){
            idxType child = G->out[j];
            printf("%d (%d),", child,part[child] );
        }
        printf(">\n");
}

int printPartWeights(dgraph* G, idxType* part)
{
    int nbpart = 0;

    idxType maxpart = part[1];
    idxType i;
    for (i = 1; i<= G->nVrtx; i++) {
        maxpart = part[i] > maxpart  ? part[i] : maxpart;
    }
    nbpart = maxpart+1;

    idxType* partsize = (idxType*) calloc(nbpart, sizeof(idxType));
    idxType minsize = INT_MAX, maxsize = 0;

    for (i = 1; i <= G->nVrtx; i++) {
        partsize[part[i]] += G->vw[i];
    }

    for (i = 0; i < nbpart; i++) {
        minsize = minsize < partsize[i] ? minsize : partsize[i];
        maxsize = maxsize < partsize[i] ? partsize[i] : maxsize;
    }

    free(partsize);
    return (int) maxsize;
}



ecType edgeCut(dgraph* G, idxType* part)
{
    int frmt = G->frmt;
    /*Return the edge cut of the partition*/
    ecType edgecut = 0;
    idxType i,j;
    for (i = 1; i<=G->nVrtx; i++)
        for (j=G->outStart[i]; j<=G->outEnd[i]; j++) {
            if (part[G->out[j]] != part[i]) {
                if (frmt & DG_FRMT_EC)
                    edgecut += G->ecOut[j];
                else
                    edgecut++;
            }
        }
    return edgecut;
}
ecType volume(dgraph* G, idxType* part, idxType nbPart)
{
    ecType vol = 0.0;
    ecType** topart = (ecType**) calloc((G->nVrtx+1),sizeof(ecType*));
    ecType i, j, neigh;
    for (i=0; i <= G->nVrtx; ++i)
        topart[i] = (ecType*) calloc((nbPart), sizeof(ecType));
    for (i=1; i <= G->nVrtx; ++i) {
        for (j = G->outStart[i]; j <= G->outEnd[i]; ++j) {
            neigh = G->out[j];
            if (part[i] != part[neigh]) {
                if (topart[i][part[neigh]] == 0) {
                    if (G->frmt & DG_FRMT_EC)
                        topart[i][part[neigh]] = G->ecOut[j];
                    else
                        ++topart[i][part[neigh]];
                }
            }
        }
    }
    for (i=1; i <= G->nVrtx; ++i) {
        for (j=0; j < nbPart; ++j) {
            vol += topart[i][j];
        }
    }
    return vol;
}

void reverseGraph(dgraph *G)
{

    idxType *inStart = G->inStart;
    idxType *inEnd = G->inEnd;
    idxType *in = G->in;
    ecType *ecIn = G->ecIn;

    idxType *outStart = G->outStart;
    idxType *outEnd = G->outEnd;
    idxType *out = G->out;
    ecType *ecOut = G->ecOut;

    G->inStart = outStart;
    G->inEnd = outEnd;
    G->in = out;
    G->ecIn = ecOut;

    G->outStart = inStart;
    G->outEnd = inEnd;
    G->out = in;
    G->ecOut = ecIn;
}

void randomizeWeights(dgraph *G, vwType vmin, vwType vmax, ecType emin, ecType emax){
    idxType i,j;
    G->totvw = 0;
    G->maxVW = 0;
    G->totec = 0;
    G->maxEC = 0;
    for (i = 1; i <= G->nVrtx; i++) {
        G->vw[i] = vmin + uRandom((int) (vmax - vmin + 1));
        G->totvw += G->vw[i];
        G->maxVW = G->maxVW < G->vw[i] ? G->vw[i] : G->maxVW;
        //printf("Weight %d : %d\n", (int) i, (int) G->vw[i]);
        for (j = G->inStart[i]; j <= G->inEnd[i]; j++) {
            G->ecIn[j] = emin + uRandom((int) (emax - emin + 1));
            G->totec += G->ecIn[j];
            G->maxEC = G->maxEC < G->ecIn[j] ? G->ecIn[j] : G->maxEC;
            //printf("Weight %d -> %d : %d\n", (int) G->in[j], (int) i, (int) G->ecIn[j]);
        }
    }
    fillOutFromIn(G);
}

void applyCCR(dgraph *G, double CCR){
    double mult = CCR * G->totvw / G->totec;
    idxType i,j;
    G->totec = 0;
    for (i = 1; i <= G->nVrtx; i++) {
        for (j = G->inStart[i]; j <= G->inEnd[i]; j++) {
            G->ecIn[j] = (ecType) round(mult * G->ecIn[j]);
            G->totec += G->ecIn[j];
        }
    }
    fillOutFromIn(G);
}

idxType nbPart(dgraph* G, idxType* part, vwType* partsize)
{
    /*Assume partsize in already allocated
      Compute partsize[i] = number of vertex is part i
      and return the number of partition in part*/
    idxType nbpart = 0;
    int frmt = G->frmt;
    idxType i;
    for (i = 1; i <= G->nVrtx; i++)
        if (nbpart < part[i]) {
            nbpart = part[i];
        }
    nbpart++;

    for (i = 0; i<nbpart; i++)
        partsize[i] = 0;

    for (i = 1; i<=G->nVrtx; i++)
        if (frmt & DG_FRMT_VW)
            partsize[part[i]] += G->vw[i];
        else
            partsize[part[i]]++;
    return nbpart;
}

