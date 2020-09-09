#include "dgraphTraversal.h"


void topologic_natural_with_part(dgraph *G, idxType* toporder, idxType *part, idxType *toporderpart, idxType nbpart){
    idxType i,j, cur_part;
    idxType to = 1, nbinput, node;

    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (nbinleft == NULL)
        u_errexit("nbinleft allocation failed\n");
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;

    //Build heap of ready tasks
    idxType** heaps = (idxType**) umalloc(sizeof(idxType*)*nbpart, "heaps");
    idxType** inheaps = (idxType**) umalloc(nbpart * sizeof(idxType*), "inheaps");
    idxType* hsize = (idxType*) calloc(nbpart, sizeof(idxType));
    ecType* priority = (ecType*) umalloc((G->nVrtx+1) * sizeof(ecType), "priority");
    if (hsize == NULL)
        u_errexit("hsize alloc failed\n");
    for (i = 0; i < nbpart; i++) {
        heaps[i] = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "heaps[i]");
        inheaps[i] = (idxType*) calloc(G->nVrtx+1, sizeof(idxType));
        if (inheaps[i] == NULL)
            u_errexit("inheaps[i] alloc failed\n");
    }

    for (i = 1; i <= G->nVrtx; i++){
        nbinput = 0;
        for (j = G->inStart[i]; j <= G->inEnd[i]; j++)
            if (part[G->in[j]] == part[i])
                nbinput++;
        if (nbinput == 0){
            hsize[part[i]]++;
            heaps[part[i]][hsize[part[i]]] = i;
            inheaps[part[i]][i] = hsize[part[i]];
        }
        priority[i] = G->nVrtx - i;
    }
    for (j = 0; j < nbpart; j++) {
        cur_part = toporderpart[j];
        heapBuild(heaps[cur_part], priority, hsize[cur_part], inheaps[cur_part]);
        while (hsize[cur_part] > 0){
            if (to > G->nVrtx)
                u_errexit("Proof topologic_natural_with_part = %d\n", to);
            node = heapExtractMax(heaps[cur_part], priority, &hsize[cur_part], inheaps[cur_part]);
            if (nbinleft [node] != 0){
                u_errexit("not ready node selected %d\n", node);
            }
            toporder[to++] = node;
            for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
                idxType succ = G->out[i];
                nbinleft[succ]--;
                if ((nbinleft[succ] == 0)&&(part[succ] == cur_part))
                    heapInsert(heaps[cur_part], priority, &hsize[cur_part], succ, inheaps[cur_part]);
                else if (nbinleft[succ]<0)
                    u_errexit("topologic_natural_with_part: negative indegree\n");
            }
        }
    }
    free(nbinleft);
    for (i = 0; i < nbpart; i++) {
        free(heaps[i]);
        free(inheaps[i]);
    }
    free(heaps);
    free(inheaps);
    free(hsize);
    free(priority);
    if (to != G->nVrtx+1)
        u_errexit("topologic_natural_with_part : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);

}

/*
void DFStopsort_with_part(dgraph* G, idxType *part, int nbpart, idxType *toporder)
{
    // Assume that toporder is already allocated
    //  Fill toporder to have a topological order of nodes in G
    // with DFS respecting the convexity of the partitioning
    // (every nodes of one part, then another part...)
    idxType i,j;
    idxType* toporderpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "toporderpart");
    topSortOnParts(G, part, toporderpart, nbpart);
    idxType* topsortpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "topsortpart");

    for (i = 0; i<nbpart; i++) {
        topsortpart[toporderpart[i]] = i;
    }

    idxType to = 1;
    int current_part = 0;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType nbready = sourcesList(G, ready);

    idxType* nbinleft;
    nbinleft = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "nbinleft");
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    while (nbready > 0)
    {
        idxType node = -1;
        for (i = 1; i<=nbready; i++) {
            if (part[ready[nbready-i]] == toporderpart[current_part]) {
                node = ready[nbready-i];
                break;
            }
        }
        if (node == -1) {
            current_part++;
            continue;
        }
        for (j = nbready-i; j<nbready-1; j++) {
            ready[j] = ready[j+1]; //BU2JH: a regarder
        }

        nbready--;
        toporder[to++] = node;

        //printf("Loop nbready = %d, node = %d, to = %d, part = %d\n", (int) nbready, (int) node, (int) to, part[node]);
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
        }
    }
    if (to != G->nVrtx +1)
        u_errexit("DFStopsort_with_part: not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(toporderpart);
    free(topsortpart);
    free(ready);
}
*/

void DFStopsort_with_part(dgraph *G, idxType* toporder, idxType *part, idxType *toporderpart, idxType nbpart)
{
    // printf("in new DFStopsortsort_with_part. nbpart = %d\n", nbpart);
    idxType i,j;
    // idxType* toporderpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "toporderpart");
    // topSortOnParts(G, part, toporderpart, nbpart);
    idxType* topsortpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "topsortpart");
    idxType* inpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "inpart");

    for (i = 0; i<nbpart; i++) {
        topsortpart[toporderpart[i]] = i;
    }

    idxType to = 1;
    int current_part = 0;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    // idxType nbready = sourcesList(G, ready);
    idxType nbready = 0;
    idxType* nbinleft;
    nbinleft = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "nbinleft");
    // int c1=0, c2 =0;
    for (i = 1; i <= G->nVrtx; i++) {
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
        // if (part[i] == 0)
        //     c1++;
        // if (part[i]== 1)
        //     c2++;
    }
    // printf("c1c2 %d %d\n", c1, c2);
    idxType succ;
    idxType node;
    while (current_part < nbpart) {
        nbready = 0;
        for (i=1; i<=G->nVrtx; ++i) {
            if (part[i] == toporderpart[current_part]) {
                nbinleft[i] = 0;
                for (j= G->inStart[i]; j<=G->inEnd[i]; ++j) {
                    if (part[G->in[j]] == toporderpart[current_part])
                        ++nbinleft[i];
                }
                if (nbinleft[i] == 0) {
                    ready[nbready++] = i;
                }
            }
        }
        node = -1;
        if (nbready == 0){
            // this is possible. e.g. nbpart = 2; trying merging
            // u_errexit("part %d is empty!", toporderpart[current_part]);
        }
        while (nbready > 0) {
            if (part[ready[nbready-1]] != toporderpart[current_part]) {
                u_errexit("node not in part, in ready array!");
            }
            node = ready[nbready-1];
            -- nbready;
            toporder[to++] = node;
            
            //printf("Loop nbready = %d, node = %d, to = %d, part = %d\n", (int) nbready, (int) node, (int) to, part[node]);
            for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
                succ = G->out[i];
                -- nbinleft[succ];
                if (nbinleft[succ] == 0 && part[succ] == toporderpart[current_part])
                    ready[nbready++] = succ;
            }
        }
        ++ current_part;
    }
    if (to != G->nVrtx +1)
        u_errexit("DFStopsort_with_part: not every node covered! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    // free(toporderpart);
    free(topsortpart);
    free(ready);
}

void BFStopsort_with_part(dgraph *G, idxType* toporder, idxType *part, idxType *toporderpart, idxType nbpart)
{
    // printf("in new DFStopsortsort_with_part. nbpart = %d\n", nbpart);
    idxType i,j;
    // idxType* toporderpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "toporderpart");
    // topSortOnParts_BFS(G, part, toporderpart, nbpart);
    idxType* topsortpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "topsortpart");
    idxType* inpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "inpart");

    for (i = 0; i<nbpart; i++) {
        topsortpart[toporderpart[i]] = i;
    }

    idxType to = 1;
    int current_part = 0;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    // idxType nbready = sourcesList(G, ready);
    idxType nbready = 0;
    idxType* nbinleft;
    nbinleft = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "nbinleft");
    // int c1=0, c2 =0;
    for (i = 1; i <= G->nVrtx; i++) {
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
        // if (part[i] == 0)
        //     c1++;
        // if (part[i]== 1)
        //     c2++;
    }
    // printf("c1c2 %d %d\n", c1, c2);
    idxType succ;
    idxType node;
    idxType beginready = 0;
    while (current_part < nbpart) {
        nbready = 0;
        beginready = 0;
        for (i=1; i<=G->nVrtx; ++i) {
            if (part[i] == toporderpart[current_part]) {
                nbinleft[i] = 0;
                for (j= G->inStart[i]; j<=G->inEnd[i]; ++j) {
                    if (part[G->in[j]] == toporderpart[current_part])
                        ++nbinleft[i];
                }
                if (nbinleft[i] == 0) {
                    ready[nbready++] = i;
                }
            }
        }
        node = -1;
        if (nbready == 0){
            // this is possible. e.g. nbpart = 2; trying merging
            // u_errexit("part %d is empty!", toporderpart[current_part]);
        }
        while (nbready > beginready) {
            if (part[ready[beginready]] != toporderpart[current_part]) {
                u_errexit("node not in part, in ready array!");
            }
            node = ready[beginready];
            ++beginready;
            toporder[to++] = node;
            
            //printf("Loop nbready = %d, node = %d, to = %d, part = %d\n", (int) nbready, (int) node, (int) to, part[node]);
            for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
                succ = G->out[i];
                -- nbinleft[succ];
                if (nbinleft[succ] == 0 && part[succ] == toporderpart[current_part])
                    ready[nbready++] = succ;
            }
        }
        ++ current_part;
    }
    if (to != G->nVrtx +1)
        u_errexit("BFStopsort_with_part: not every node covered! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    // free(toporderpart);
    free(topsortpart);
    free(ready);
}




/*
void BFStopsort_with_part(dgraph* G, idxType *part, int nbpart, idxType *toporder)
{
    // Assume that toporder is already allocated
    //  Fill toporder to have a topological order of nodes in G
    // with DFS respecting the convexity of the partitioning
    // (every nodes of one part, then another part...)
    idxType i,j;
    idxType* toporderpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "toporderpart");
    topSortOnParts(G, part, toporderpart, nbpart);
    idxType* topsortpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "topsortpart");

    for (i = 0; i<nbpart; i++) {
        topsortpart[toporderpart[i]] = i;
    }

    idxType to = 1;
    int current_part = 0;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType nbready = sourcesList(G, ready);

    idxType* nbinleft;
    nbinleft = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "nbinleft");
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;;
    while (nbready > 0)
    {
        idxType node = -1;
        for (i = 1; i<=nbready; i++) {
            if (part[ready[nbready-i]] == toporderpart[current_part]) {
                node = ready[nbready-i];
                break;
            }
        }
        if (node == -1) {
            current_part++;
            continue;
        }
        for (j = nbready-i; j<nbready-1; j++) {
            ready[j] = ready[j+1];  //BU2JH: a regarder
        }

        nbready--;
        toporder[to++] = node;

        //printf("Loop nbready = %d, node = %d, to = %d, part = %d\n", (int) nbready, (int) node, (int) to, part[node]);
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0) {
                //      printf("succ = %d is ready\n",succ);
                for (j = nbready-1; j>=0; j--)
                    ready[j+1] = ready[j];  //BU2JH: a regarder
                ready[0] = succ;
                nbready++;
            }
        }
    }
    if (to != G->nVrtx + 1)
        u_errexit("BFStopsort_with_part, not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(toporderpart);
    free(topsortpart);
    free(ready);
}
*/
void randDFStopsort_with_part(dgraph* G, idxType *part, int nbpart, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a topological order of nodes in G
    with DFS respecting the convexity of the partitioning
    (every nodes of one part, then another part...)*/
    idxType i,j, itmp;
    idxType* toporderpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "toporderpart");
    topSortOnParts(G, part, toporderpart, nbpart);
    idxType* topsortpart = (idxType*) umalloc (sizeof(idxType) * nbpart, "topsortpart");

    for (i = 0; i<nbpart; i++) {
        topsortpart[toporderpart[i]] = i;
    }

    idxType to = 1;
    int current_part = 0;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1),  "ready");
    idxType nbready = sourcesList(G, ready);
    idxType* shuffle = (idxType*) umalloc(sizeof(idxType)*(G->maxoutdegree), "shuffle");

    idxType* nbinleft;
    nbinleft = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "nbinleft");
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    while (nbready > 0)
    {
        idxType node = -1;
        for (i = 1; i<=nbready; i++) {
            if (part[ready[nbready-i]] == toporderpart[current_part]) {
                node = ready[nbready-i];
                break;
            }
        }
        if (node == -1) {
            current_part++;
            continue;
        }
        for (j = nbready-i; j<nbready-1; j++) {
            ready[j] = ready[j+1]; /*BU2JH: a regarder*/
        }

        nbready--;
        toporder[to++] = node;
        shuffleTab(G->outStart[node], G->outEnd[node], shuffle);

        //printf("Loop nbready = %d, node = %d, to = %d, part = %d\n", (int) nbready, (int) node, (int) to, part[node]);
        for (itmp = 0; itmp <= G->outEnd[node]-G->outStart[node]; itmp++)
        {
            i = shuffle[itmp];
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
        }
    }
    if (to != G->nVrtx +1)
        u_errexit("RandDFStopsort_with_part, not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(shuffle);
    free(nbinleft);
    free(toporderpart);
    free(topsortpart);
    free(ready);
}


void checkTopsort(dgraph *G, idxType *toporder)
{
    idxType i, j, neigh, nVrtx = G->nVrtx;
    idxType *outEnd = G->outEnd, *outStart = G->outStart, *out = G->out;
    idxType *topindex = (idxType *) umalloc(sizeof(idxType)* (nVrtx+1), "topindex");

    for (i = 1; i <= nVrtx; i++)
        topindex[toporder[i]] = i;

    for (i = 1; i <= nVrtx; i ++)
    {
        for ( j = outStart[i]; j<= outEnd[i]; j++)
        {
            neigh = out[j];
            if (topindex[i] >= topindex[neigh])
                u_errexit("The topological order is not correct %.0f %.0f\n", 1.0*topindex[i], 1.0*topindex[neigh]);
        }
    }
    free(topindex);
}

void topSortOnParts(dgraph *G, idxType *part, idxType *toporder, idxType nbpart)
{
    // Does DFS approach
    /*Fill toporder assuming it's already allocated*/
    idxType i, j, k;
    idxType** outpart = (idxType**) umalloc(sizeof(idxType*)*nbpart, "outpart");
    idxType* nbout = (idxType*) umalloc(sizeof(idxType)*nbpart, "nbout");
    idxType** inpart = (idxType**) umalloc(sizeof(idxType*)*nbpart, "inpart");
    idxType* nbin = (idxType*) umalloc(sizeof(idxType)*nbpart, "nbin");
    for (i=0; i<nbpart; i++) {
        outpart[i] = (idxType*) umalloc(sizeof(idxType)*nbpart, "outpart[i]");
        nbout[i] = 0;
        inpart[i] = (idxType*) umalloc(sizeof(idxType)*nbpart, "inpart[i]");
        nbin[i] = 0;
    }
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
    // printf("xxtr\n");
    idxType* ready = (idxType*) umalloc(sizeof(idxType)*(nbpart), "ready");
    idxType* nbinleft = (idxType*) umalloc(sizeof(idxType)*(nbpart), "nbinleft");
    int nbready = 0;
    for (i=0; i<nbpart; i++) {
        nbinleft[i] = nbin[i];
        if (nbin[i] == 0)
            ready[nbready++] = i;
    }

    int to = 0;

    while (nbready > 0) {
        idxType pno = ready[nbready-1];
        // printf("pno %d\n", pno);
        nbready--;
        toporder[to++] = pno;
        for (i = 0; i < nbout[pno]; i++) {
            idxType succ = outpart[pno][i];
            nbinleft[succ]--;
            // printf("succ: %d, left %d\n", succ, nbinleft[succ]);
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
            // else if(nbinleft[succ] < 0)
            //     printf("something wrong at nbinleft\n");
        }
    }
    if (to != nbpart) {
        u_errexit("In topsortPart, not all sort: to = %d, nbpart = %d\n", to, nbpart);
    }
    for (i=0; i<nbpart; i++) {
        free(outpart[i]);
        free(inpart[i]);
    }
    free(ready);
    free(nbinleft);
    free(nbin);
    free(nbout);
    free(inpart);
    free(outpart);
}


void topSortOnParts_BFS(dgraph *G, idxType *part, idxType *toporder, idxType nbpart)
{
    /*Fill toporder assuming it's already allocated*/
    idxType i, j, k;
    idxType** outpart = (idxType**) umalloc(sizeof(idxType*)*nbpart, "outpart");
    idxType* nbout = (idxType*) umalloc(sizeof(idxType)*nbpart, "nbout");
    idxType** inpart = (idxType**) umalloc(sizeof(idxType*)*nbpart, "inpart");
    idxType* nbin = (idxType*) umalloc(sizeof(idxType)*nbpart, "nbin");
    for (i=0; i<nbpart; i++) {
        outpart[i] = (idxType*) umalloc(sizeof(idxType)*nbpart, "outpart[i]");
        nbout[i] = 0;
        inpart[i] = (idxType*) umalloc(sizeof(idxType)*nbpart, "inpart[i]");
        nbin[i] = 0;
    }
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
    // printf("xxtr\n");
    idxType* ready = (idxType*) umalloc(sizeof(idxType)*(nbpart), "ready");
    idxType* nbinleft = (idxType*) umalloc(sizeof(idxType)*(nbpart), "nbinleft");
    int nbready = 0;
    for (i=0; i<nbpart; i++) {
        nbinleft[i] = nbin[i];
        if (nbin[i] == 0)
            ready[nbready++] = i;
    }

    int to = 0;
    int beginready = 0;
    while (nbready > beginready) {
        idxType pno = ready[beginready];
        // printf("pno %d\n", pno);
        ++ beginready;
        toporder[to++] = pno;
        for (i = 0; i < nbout[pno]; i++) {
            idxType succ = outpart[pno][i];
            nbinleft[succ]--;
            // printf("succ: %d, left %d\n", succ, nbinleft[succ]);
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
            // else if(nbinleft[succ] < 0)
            //     printf("something wrong at nbinleft\n");
        }
    }
    if (to != nbpart) {
        u_errexit("In topsortPart, not all sort: to = %d, nbpart = %d\n", to, nbpart);
    }
    for (i=0; i<nbpart; i++) {
        free(outpart[i]);
        free(inpart[i]);
    }
    free(ready);
    free(nbinleft);
    free(nbin);
    free(nbout);
    free(inpart);
    free(outpart);
}

void randTopSortOnParts(dgraph *G, idxType *part, idxType *toporder, idxType nbpart)
{
    /*Fill toporder assuming it's already allocated*/
    idxType i, j, k, ip;
    idxType** outpart = (idxType**) umalloc(sizeof(idxType*)*nbpart, "outpart");
    idxType* nbout = (idxType*) umalloc(sizeof(idxType)*nbpart, "nbout");
    idxType** inpart = (idxType**) umalloc(sizeof(idxType*)*nbpart, "inpart");
    idxType* nbin = (idxType*) umalloc(sizeof(idxType)*nbpart, "nbin");
    for (i=0; i<nbpart; i++) {
        outpart[i] = (idxType*) umalloc(sizeof(idxType)*nbpart, "outpart[i]");
        nbout[i] = 0;
        inpart[i] = (idxType*) umalloc(sizeof(idxType)*nbpart, "inpart[i]");
        nbin[i] = 0;
    }
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

    idxType* ready = (idxType*) umalloc(sizeof(idxType)*(nbpart), "ready");
    idxType* nbinleft = (idxType*) umalloc(sizeof(idxType)*(nbpart), "nbinleft");
    int nbready = 0;
    for (i=0; i<nbpart; i++) {
        nbinleft[i] = nbin[i];
        if (nbin[i] == 0)
            ready[nbready++] = i;
    }

    int to = 0;
    idxType* shuffle = (idxType *) umalloc(nbpart * sizeof(idxType), "shuffle");
    while (nbready > 0) {
        idxType pno = ready[nbready-1];
        nbready--;
        toporder[to++] = pno;
        if (nbout[pno]) {
            shuffleTab(0, nbout[pno]-1, shuffle);
            for (ip = 0; ip < nbout[pno]; ip++) {
                i = shuffle[ip];
                idxType succ = outpart[pno][i];
                nbinleft[succ]--;
                if (nbinleft[succ] == 0) {
                    ready[nbready++] = succ;
                }
            }
        }
    }
    for (i=0; i<nbpart; i++) {
        free(outpart[i]);
        free(inpart[i]);
    }
    free(nbinleft);
    free(shuffle);
    free(ready);
    free(nbin);
    free(nbout);
    free(inpart);
    free(outpart);
    if (to != nbpart) {
        u_errexit("In topsortPart, not every nodes are sorted: to = %d, nbpart = %d\n", to, nbpart);
    }
}

void topsortPriorities(dgraph *G, idxType *toporder, ecType *priority)
{
    /*Assume that toporder is already allocated an priorities filled
     Fill toporder to have a topological order of nodes in G
    with DFS*/
    idxType to = 1, i, node;
    // printf("in topsortPriorities G->nVrtx %d\n", G->nVrtx);
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (nbinleft == NULL)
        u_errexit("nbinleft alloc error\n");
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    //Build heap of ready tasks
    idxType* heap = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "heap");
    idxType* inheap = (idxType*) calloc(G->nVrtx+1, sizeof(idxType));
    if (inheap == NULL)
        u_errexit("inheap alloc error\n");
    idxType hsize = 0;
    for (i =0; i < G->nbsources; i++){
        hsize++;
        heap[hsize] = G->sources[i];
        inheap[G->sources[i]] = hsize;
    }
    heapBuild(heap, priority, hsize, inheap);

    while (hsize > 0) {
        if (to > G->nVrtx)
            u_errexit("Proof topsortPrioritiesto = %d\n", to);
        node = heapExtractMax(heap, priority, &hsize, inheap);
        toporder[to++] = node;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0)
                heapInsert(heap, priority, &hsize, succ, inheap);
            else if (nbinleft[succ]<0)
                u_errexit("topsortPriorities: negative indegree\n");
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("topsortPriorities : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(heap);
    free(inheap);
}

void DFStopsort(dgraph* G, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a topological order of nodes in G
    with DFS*/
    idxType to = 1;
    idxType* ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i;
    idxType nbready = G->nbsources;
    for (i = 0; i < nbready; i++)
        ready[i] = G->sources[i];
    sourcesList(G, ready);
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (nbinleft == NULL)
        u_errexit("nbinleft alloc error\n");

    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    while (nbready > 0) {
        if (to > G->nVrtx)
            u_errexit("Proof DFStopsort to = %d\n", to);
        idxType node = ready[nbready-1];
        nbready--;
        toporder[to++] = node;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
            else if (nbinleft[succ]<0)
                u_errexit("DFStopsort: negative indegree\n");
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("DFStopsort : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(ready);
}

void BFStopsort(dgraph* G, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a topological order of nodes in G
    with DFS*/
    idxType to = 1;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i;
    idxType endready = sourcesList(G, ready) - 1;
    idxType beginready = 0;
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (nbinleft == NULL)
        u_errexit("nbinleft alloc error\n");

    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    while (endready >= beginready) {
        if (to > G->nVrtx)
            u_errexit("Proof to = %d\n", to);
        idxType node = ready[beginready];
        beginready++;
        toporder[to++] = node;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0) {
                ready[++endready] = succ;
            }
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("BFStopsort : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(ready);
}

void DFSsort(dgraph* G, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a topological order of nodes in G
    with DFS*/
    idxType to = 1;
    idxType* ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i;
    idxType nbready = sourcesList(G, ready);
    idxType* done = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (done == NULL)
        u_errexit("done alloc error\n");
    while (nbready > 0)
    {
        if (to > G->nVrtx)
            u_errexit("Proof to = %d\n", to);
        idxType node = ready[nbready-1];
        nbready--;
        toporder[to++] = node;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++)
        {
            idxType succ = G->out[i];
            if (done[succ] == 1)
                continue;
            ready[nbready++] = succ;
            done[succ] = 1;
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("DFSsort : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(done);
    free(ready);
}

void BFSsort(dgraph* G, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a topological order of nodes in G
     with DFS*/
    idxType to = 1;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i;
    idxType endready = sourcesList(G, ready) - 1;
    idxType beginready = 0;
    idxType* done = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (done == NULL)
        u_errexit("done alloc error\n");

    while (endready >= beginready)
    {
        if (to > G->nVrtx)
            u_errexit("Proof to = %d\n", to);
        idxType node = ready[beginready];
        beginready++;
        toporder[to++] = node;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            if (done[succ] == 1)
                continue;
            ready[++endready] = succ;
            done[succ] = 1;
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("BFSsort : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(done);
    free(ready);
}

void randDFStopsort(dgraph *G, idxType *toporder)
{
    /*Assume that toporder is already allocated
      Fill toporder to have a random topological order of nodes in G*/
    idxType to = 1;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i,itmp;
    idxType nbready = sourcesList(G, ready);
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (nbinleft == NULL)
        u_errexit("nbinleft alloc error\n");
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    idxType* shuffle = (idxType*) umalloc(sizeof(idxType)*(G->maxoutdegree), "shuffle");
    while (nbready > 0)
    {
        idxType node = ready[nbready-1];
        nbready--;
        toporder[to++] = node;
        shuffleTab(G->outStart[node], G->outEnd[node], shuffle);
        for (itmp = 0; itmp <= G->outEnd[node]-G->outStart[node]; itmp++) {
            i = shuffle[itmp];
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
        }
    }
    free(shuffle);
    free(nbinleft);
    free(ready);
}

void randBFStopsort(dgraph* G, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a topological order of nodes in G
    with DFS*/
    idxType to = 1;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i,itmp;
    idxType endready = sourcesList(G, ready) - 1;
    idxType beginready = 0;
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (nbinleft == NULL)
        u_errexit("nbinleft alloc error\n");
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    idxType* shuffle = (idxType*) umalloc(sizeof(idxType)*(G->maxoutdegree), "shuffle");
    while (endready >= beginready) {
        if (to > G->nVrtx)
            u_errexit("Proof to = %d\n", to);
        idxType node = ready[beginready];
        beginready++;
        toporder[to++] = node;
        shuffleTab(G->outStart[node], G->outEnd[node], shuffle);
        for (itmp = 0; itmp <= G->outEnd[node]-G->outStart[node]; itmp++) {
            i = shuffle[itmp];
            idxType succ = G->out[i];
            assert(succ <= G->nVrtx);
            nbinleft[succ]--;
            if (nbinleft[succ] == 0) {
                ready[++endready] = succ;
            }
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("randBFStopsort : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(shuffle);
    free(ready);
}

void randDFSsort(dgraph *G, idxType *toporder)
{
    /*Assume that toporder is already allocated
      Fill toporder to have a random topological order of nodes in G*/
    idxType to = 1;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i,itmp;
    idxType nbready = sourcesList(G, ready);
    idxType* done = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (done == NULL)
        u_errexit("done alloc error\n");
    idxType* shuffle = (idxType*) umalloc(sizeof(idxType)*(G->maxoutdegree+1), "shuffle");
    while (nbready > 0) {
        idxType node = ready[nbready-1];
        nbready--;
        toporder[to++] = node;
        shuffleTab(G->outStart[node], G->outEnd[node], shuffle);
        for (itmp = 0; itmp <= G->outEnd[node]-G->outStart[node]; itmp++) {
            i = shuffle[itmp];
            idxType succ = G->out[i];
            if (done[succ] == 1)
                continue;
            ready[nbready++] = succ;
            done[succ] = 1;
        }
    }
    free(shuffle);
    free(done);
    free(ready);
}

void randBFSsort(dgraph* G, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a topological order of nodes in G
    with DFS*/
    idxType to = 1;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i,itmp;
    idxType endready = sourcesList(G, ready) - 1;
    idxType beginready = 0;
    idxType* done = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (done == NULL)
        u_errexit("done alloc error\n");
    idxType* shuffle = (idxType*) umalloc(sizeof(idxType)*(G->maxoutdegree), "shuffle");
    while (endready >= beginready)
    {
        if (to > G->nVrtx)
            u_errexit("Proof to = %d\n", to);
        idxType node = ready[beginready];
        beginready++;
        toporder[to++] = node;
        shuffleTab(G->outStart[node], G->outEnd[node], shuffle);
        for (itmp = 0; itmp <= G->outEnd[node]-G->outStart[node]; itmp++) {
            i = shuffle[itmp];
            idxType succ = G->out[i];
            if (done[succ] == 1)
                continue;
            ready[++endready] = succ;
            done[succ] = 1;
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("randBFSsort : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(done);
    free(shuffle);
    free(ready);
}


void randTopsort(dgraph* G, idxType *toporder)
{
    /*Assume that toporder is already allocated
     Fill toporder to have a random topological order of nodes in G
    */
    idxType to = 1;
    idxType* ready;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType i;
    idxType nbready = sourcesList(G, ready);
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (nbinleft == NULL)
        u_errexit("nbinleft alloc error\n");
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    while (nbready > 0) {
        if (to > G->nVrtx)
            u_errexit("Proof to = %d\n", to);
        idxType node_idx = uRandom(nbready);
        idxType node = ready[node_idx];
        ready[node_idx] = ready[nbready-1];
        nbready--;
        toporder[to++] = node;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0) {
                ready[nbready++] = succ;
            }
        }
    }
    if (to != G->nVrtx + 1)
        u_errexit("randTopsort : Not every node concerned! to = %d, nVrtx = %d\n", to, G->nVrtx);
    free(nbinleft);
    free(ready);
}


int gcd ( int a, int b )
{
    int c;
    while ( a != 0 ) {
        c = a; a = b%a;  b = c;
    }
    return b;
}

void mixtopsort(dgraph* G, idxType* toporder, int priority, int first)
{
    //first==0 --> bfs first
    //first==1 --> dfs first
    idxType* ready = (idxType*) umalloc(sizeof(idxType)*(G->nVrtx+1), "ready");
    idxType* mark = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    if (nbinleft == NULL || mark == NULL)
        u_errexit("nbinleft or mark alloc error\n");
    idxType readyend = sourcesList(G,ready) - 1;
    idxType readybegin = 0;
    idxType bfrac=0, dfrac=0;
    idxType biter=0, diter=0;
    idxType node=-1;
    idxType i,j,to=1;
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;

    i=gcd(priority,100-priority);
    if (i==0) {
        bfrac=priority/100;
        dfrac=(100-priority)/100;
    }
    else{
        bfrac=priority/i;
        dfrac=(100-priority)/i;
    }
    if (first)
        diter=dfrac;
    else
        biter=bfrac;

    while (readyend>=readybegin) {
        if (to > G->nVrtx+1) {
            printf("Proof to = %d\n", to);
            break;
        }

        //pick from bfs priority
        if (biter>=0) {
            node = ready[readybegin];
            readybegin++;
            biter--;
            if (biter<0)
                diter=dfrac;
        }
        else { //pick from dfs priority
            node = ready[readyend];
            diter--;
            readyend--;
            if (diter<0)
                biter=bfrac;
        }
        mark[node]=1;
        toporder[to++] = node;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0) {
                readyend++;
                ready[readyend] = succ;
            }
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("Mixtopsort : Not every node concerned!  to = %d, nVrtx = %d\n", to, G->nVrtx);

    free(nbinleft);
    free(mark);
    free(ready);
}
