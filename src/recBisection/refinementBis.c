#include "refinementBis.h"
#include "dgraph.h"


int getBiggerPart(idxType inPartIdx, idxType outPartIdx, idxType *partsize, idxType** heap, idxType* hsize,  double* lb_pw, double* ub_pw, vwType maxvw, ecType *gain, int balNotOk, MLGP_option opt){
    /* This function is called in refinementPostOrder_max
     * It returns the index of the bigger part we wanna send a node from
     * Return -1 if we have to break the while loop
     */
    int biggerPart = -1;
    // this is only for 2 parts
    double partRatioIn,partRatioOut;
    partRatioIn = (double) partsize[inPartIdx] / ub_pw[inPartIdx];
    partRatioOut = (double) partsize[outPartIdx] / ub_pw[outPartIdx];
    ecType gainin, gainout;
    if (balNotOk || opt.refinement == REF_bFM_MAXW) {
        if (partRatioIn > partRatioOut){
            biggerPart = inPartIdx;
            if(hsize[inPartIdx] == 0 )
                return -1;
        }
        else{
            biggerPart = outPartIdx;
            if (hsize[outPartIdx] == 0)
                return -1;
        }
    }
    else {/*opt.refinement == REF_bFM*/
        if (hsize[inPartIdx] == 0){
            if ((partsize[inPartIdx] <= ub_pw[inPartIdx] - maxvw) && (partsize[outPartIdx] >= lb_pw[outPartIdx] + maxvw && hsize[outPartIdx] >= 1))
                biggerPart = outPartIdx;
            else
                return -1;
        }
        else if (hsize[outPartIdx] == 0){
            if ((partsize[outPartIdx] <= ub_pw[outPartIdx] - maxvw) && (partsize[inPartIdx] >= lb_pw[inPartIdx] + maxvw && hsize[inPartIdx] >= 1))
                biggerPart = inPartIdx;
            else
                return -1;
        }
        else {/*both have things in the heap*/
            if (! (partsize[inPartIdx] >= lb_pw[inPartIdx] + maxvw && partsize[outPartIdx] <= ub_pw[outPartIdx] - maxvw))
            {/* P0 cannot send*/
                biggerPart = outPartIdx;
            }
            else if (! (partsize[outPartIdx] >= lb_pw[outPartIdx] + maxvw && partsize[inPartIdx] <= ub_pw[inPartIdx] - maxvw))
            {/*P1 cannot send*/
                biggerPart = inPartIdx;
            }
            else /*both can send*/
            {
                gainin = gain[*(heap[inPartIdx]+1)];
                gainout = gain[*(heap[outPartIdx]+1)];
                if (gainin > gainout)
                    biggerPart = inPartIdx;
                else if (gainout > gainin)
                    biggerPart = outPartIdx;
                else{
                    if (partRatioIn < partRatioOut) biggerPart = outPartIdx;
                    else biggerPart = inPartIdx;
                }
            }
        }
    }
    return biggerPart;
}

static inline void bisInitializeGainAndHeap(dgraph *G, idxType inPartIdx, idxType outPartIdx, inFromOutToCnts *ioCnts, ecType* gain, idxType** heap, idxType* inheap, idxType* hsize, idxType *part, idxType* partsize, idxType* maxpartsize, vwType* maxvw, idxType* partcopy, MLGP_option opt){
    /* This function is called at the begining of refinements for bisection.
     * It initialize the ioCnts structure, the gain table, and build the heap
     * It also initializes partsize, maxpartsize, maxweight and partcopy
     */
    idxType i,j;

    partsize[inPartIdx] = partsize [outPartIdx] = 0;
    for (i=1; i<=G->nVrtx; i++){
        partcopy[i] = part[i];

        ioCnts[i].inFrom[inPartIdx] = ioCnts[i].inFrom[outPartIdx] = ioCnts[i].outTo[inPartIdx] = ioCnts[i].outTo[outPartIdx] = 0;

        for (j=G->inStart[i]; j<=G->inEnd[i]; j++) {
            ioCnts[i].inFrom[part[G->in[j]]] += G->ecIn[j];
        }

        for (j=G->outStart[i]; j<=G->outEnd[i]; j++) {
            ioCnts[i].outTo[part[G->out[j]]] += G->ecOut[j];
        }

        if(opt.debug > 1) {
            if(part[i] == inPartIdx && ioCnts[i].inFrom[outPartIdx] != (ecType) 0)
                u_errexit("bisInitializeGainAndHeap: an incoming edge in the bad direction.");
            if(part[i] == outPartIdx && ioCnts[i].outTo[inPartIdx] != (ecType) 0)
                u_errexit("bisInitializeGainAndHeap: an outgoing edge in the bad direction.");
        }

        if (part[i] == inPartIdx && ioCnts[i].outTo[inPartIdx] == 0){/*can go to outPartIdx*/
            heap[inPartIdx][++hsize[inPartIdx]] = i;
            inheap[i] = hsize[inPartIdx];
            gain[i] = ioCnts[i].outTo[outPartIdx] - ioCnts[i].inFrom[inPartIdx];
        }
        else if (part[i] == outPartIdx && ioCnts[i].inFrom[outPartIdx] == 0){/*can go to 0*/
            heap[outPartIdx][++hsize[outPartIdx]] = i;
            inheap[i] = hsize[outPartIdx];
            gain[i] = -ioCnts[i].outTo[outPartIdx] + ioCnts[i].inFrom[inPartIdx];
        }
        partsize[part[i]] += G->vw[i];
        *maxpartsize = *maxpartsize < partsize[part[i]] ? partsize[part[i]] : *maxpartsize;
    }
    *maxvw = G->maxVW;

    heapBuild(heap[inPartIdx], gain, hsize[inPartIdx], inheap);
    heapBuild(heap[outPartIdx], gain, hsize[outPartIdx], inheap);
}

static inline void bisMoveAndUpdate(dgraph *G, idxType nodeToMove, idxType* moved, idxType* nbmoved, ecType* edgecut, idxType inPartIdx, idxType outPartIdx, inFromOutToCnts* ioCnts, ecType* gain, idxType *part, idxType* partsize, idxType* maxpartsize, idxType** heap, idxType* hsize, idxType* inheap, int* anchored, MLGP_option opt) {
    /* This function is called in the different refinements for bisection.
     * It moves nodeToMove from part[nodeToMove] to 1-part[nodeToMove] in a bisection
     * Update the partsize, the heap, anchored table, moved table, ioCnts table, maxpartsize value and the edgecut
     */
    idxType j;
    idxType origin_part = part[nodeToMove];
    idxType other_part = 1-origin_part;

    part[nodeToMove] = other_part;
    partsize[origin_part] -= G->vw[nodeToMove];
    partsize[other_part] += G->vw[nodeToMove];
    anchored[nodeToMove] = 1;
    moved[*nbmoved] = nodeToMove;
    (*nbmoved)++;
    *edgecut -= gain[nodeToMove];
    *maxpartsize = partsize[inPartIdx] < partsize[outPartIdx] ? partsize[outPartIdx] : partsize[inPartIdx];

    /*nodeToMove moved from origin_part to other_part*/
    for (j=G->inStart[nodeToMove]; j<=G->inEnd[nodeToMove]; j++){
        idxType pred = G->in[j];
        if(anchored[pred]) continue; /*we do not even update their connectivity*/
        ecType ec = G->ecIn[j];
        ioCnts[pred].outTo[other_part] += ec;
        ioCnts[pred].outTo[origin_part] -= ec;

        if (opt.debug > 3)
            if(part[pred] == outPartIdx)
                u_errexit("In bisMoveAndUpdate we thought that this cannot happen (1).\n");

        if(origin_part == inPartIdx) {
            if (ioCnts[pred].outTo[origin_part] == 0){
                gain[pred] = ioCnts[pred].outTo[other_part] - ioCnts[pred].inFrom[origin_part];
                heapInsert(heap[origin_part], gain, &hsize[origin_part], pred, inheap);
            }
        }
        else/*origin_part is outPartIdx*/
        {
            if (inheap[pred] > 0)
                heapDelete(heap[other_part], gain, &hsize[other_part], pred, inheap);
            anchored[pred] = 1;
        }
    }

    for (j=G->outStart[nodeToMove]; j<=G->outEnd[nodeToMove]; j++){
        idxType succ = G->out[j];
        if(anchored[succ]) continue; /*we do not even update their connectivity*/
        ecType ec = G->ecOut[j];
        ioCnts[succ].inFrom[other_part] += ec;
        ioCnts[succ].inFrom[origin_part] -= ec;

        if (opt.debug > 3)
            if(part[succ] == inPartIdx) {
                  u_errexit("In bisMoveAndUpdate nodeToMove %d we thought that this cannot happen (2).\n", nodeToMove);
            }
        if (origin_part == inPartIdx) {
            if (inheap[succ] > 0)
                heapDelete(heap[other_part], gain, &hsize[other_part], succ, inheap);
            anchored[succ] = 1;
        }
        else/*origin_part is outPartIdx*/
        {
            if (ioCnts[succ].inFrom[origin_part] == 0){
                gain[succ] = -ioCnts[succ].outTo[origin_part] + ioCnts[succ].inFrom[other_part];
                heapInsert(heap[origin_part], gain, &hsize[origin_part], succ, inheap);
            }
        }
    }
}

int createCycleSwap(dgraph* G, idxType* part, idxType nbpart, ecType** partmatrix,
                    idxType node0, idxType node1, idxType movableto0,idxType movableto1)
{
    idxType part_save0 = part[node0];
    idxType part_save1 = part[node1];
    int res=0;
    computePartmatrix(G, part, partmatrix, node0, movableto0); // don't care about intermediate matrix since we are swapping, it may be fixed
    part[node0] = movableto0;
    int flag =computePartmatrix(G, part, partmatrix, node1, movableto1);
    part[node1] = movableto1;
    if (flag>0)
        res = thereIsCycle(nbpart, partmatrix);
    computePartmatrix(G, part, partmatrix, node1, part_save1);
    part[node1] = part_save1;
    computePartmatrix(G, part, partmatrix, node0, part_save0);
    part[node0] = part_save0;
    return res;
}

int getFlagger(idxType* node, idxType* secondNodes, idxType inPartIdx, idxType outPartIdx, dgraph *G, idxType *part, idxType* partsize, idxType nbpart, ecType** partmatrix, double* lb_pw, double* ub_pw, idxType** heap, idxType* hsize, idxType* inheap, ecType* gain){
    /* This function is called in refinementPostOrder_Swap
     * Return 1 if we can swap node[inPartIdx] and node[outPartIdx]
     * Otherwise return 2 if we can swap node[inPartIdx] and the second best in outPartIdx
     * Otherwise return 3 if we can swap the second node in inPartIdx and node[outPartIdx]
     * Otherwise return 4 if we can swap the second node in inPartIdx and the second node in outPartIdx
     * Otherwise return 0
     */
    idxType newsize[2]={-1,-1};
    int flagger = 0;
    int isThereCycle = createCycleSwap(G,part,nbpart,partmatrix,node[inPartIdx],node[outPartIdx],1-part[node[inPartIdx]],1-part[node[outPartIdx]]);

    if (!isThereCycle){
        //check new sizes to see if we can move one node from the first part to second one
        newsize[inPartIdx]=partsize[inPartIdx];
        newsize[outPartIdx]=partsize[outPartIdx];
        newsize[inPartIdx] -= G->vw[node[inPartIdx]];
        newsize[outPartIdx] += G->vw[node[inPartIdx]];
        newsize[outPartIdx] -= G->vw[node[outPartIdx]];
        newsize[inPartIdx] += G->vw[node[outPartIdx]];

        if (!((newsize[inPartIdx] < lb_pw[inPartIdx])||(newsize[inPartIdx] > ub_pw[inPartIdx])||(newsize[outPartIdx] < lb_pw[outPartIdx])||(newsize[outPartIdx] > ub_pw[outPartIdx])))
            flagger=1;
    }

    if(flagger==0){ //We can not swap node[inPartIdx] and node[outPartIdx]
        if(hsize[inPartIdx]>0) //Get the second best from heap inPartIdx
            secondNodes[inPartIdx] = heapExtractMax(heap[inPartIdx], gain, &hsize[inPartIdx], inheap);
        else
            secondNodes[inPartIdx] = -1;

        if(hsize[outPartIdx]>0) //Get the second best from heap 1
            secondNodes[outPartIdx] = heapExtractMax(heap[outPartIdx], gain, &hsize[outPartIdx], inheap);
        else
            secondNodes[outPartIdx]=-1;
    }

    // Check cycle and constraint for topFirstHeap and secondSecondHeap
    if(flagger==0 && secondNodes[outPartIdx]!=-1) {
        isThereCycle = createCycleSwap(G,part,nbpart,partmatrix,node[inPartIdx],secondNodes[outPartIdx],1-part[node[inPartIdx]],1-part[secondNodes[outPartIdx]]);
        if (!isThereCycle){
            newsize[inPartIdx]=partsize[inPartIdx];
            newsize[outPartIdx]=partsize[outPartIdx];
            newsize[inPartIdx] -= G->vw[node[inPartIdx]];
            newsize[outPartIdx] += G->vw[node[inPartIdx]];
            newsize[outPartIdx] -= G->vw[secondNodes[outPartIdx]];
            newsize[inPartIdx] += G->vw[secondNodes[outPartIdx]];

            if (!((newsize[inPartIdx] < lb_pw[inPartIdx])||(newsize[inPartIdx] > ub_pw[inPartIdx])||(newsize[outPartIdx] < lb_pw[outPartIdx])||(newsize[outPartIdx] > ub_pw[outPartIdx])))
                flagger=2;
        }
    }

    // Check cycle and constraint for secondFirstHeap and topSecondHeap
    if(flagger==0 && secondNodes[inPartIdx]!=-1) {
        isThereCycle = createCycleSwap(G,part,nbpart,partmatrix,secondNodes[inPartIdx],node[outPartIdx],1-part[secondNodes[inPartIdx]],1-part[node[outPartIdx]]);
        if (!isThereCycle){
            newsize[inPartIdx]=partsize[inPartIdx];
            newsize[outPartIdx]=partsize[outPartIdx];
            newsize[inPartIdx] -= G->vw[secondNodes[inPartIdx]];
            newsize[outPartIdx] += G->vw[secondNodes[inPartIdx]];
            newsize[outPartIdx] -= G->vw[node[outPartIdx]];
            newsize[inPartIdx] += G->vw[node[outPartIdx]];

            if (!((newsize[inPartIdx] < lb_pw[inPartIdx])||(newsize[inPartIdx] > ub_pw[inPartIdx])||(newsize[outPartIdx] < lb_pw[outPartIdx])||(newsize[outPartIdx] > ub_pw[outPartIdx])))
                flagger=3;
        }
    }

    // Check cycle and constraint for secondFirstHeap and secondSecondHeap
    if(flagger==0 && secondNodes[inPartIdx]!=-1 && secondNodes[outPartIdx]!=-1) {
        isThereCycle = createCycleSwap(G,part,nbpart,partmatrix,secondNodes[inPartIdx],secondNodes[outPartIdx],1-part[secondNodes[inPartIdx]],1-part[secondNodes[outPartIdx]]);
        if (!isThereCycle){
            newsize[inPartIdx]=partsize[inPartIdx];
            newsize[outPartIdx]=partsize[outPartIdx];
            newsize[inPartIdx] -= G->vw[secondNodes[inPartIdx]];
            newsize[outPartIdx] += G->vw[secondNodes[inPartIdx]];
            newsize[outPartIdx] -= G->vw[secondNodes[outPartIdx]];
            newsize[inPartIdx] += G->vw[secondNodes[outPartIdx]];

            if (!((newsize[inPartIdx] < lb_pw[inPartIdx])||(newsize[inPartIdx] > ub_pw[inPartIdx])||(newsize[outPartIdx] < lb_pw[outPartIdx])||(newsize[outPartIdx] > ub_pw[outPartIdx])))
                flagger=4;
        }
    }

    return flagger;
}

// SWAPS nodes between parts. Move one from first to second, move one from second to first.
void refinementPostOrder_Swap(dgraph *G, ecType** partmatrix, idxType *toporderpart,
                              double* lb_pw, double* ub_pw,
                              idxType *part, idxType nbpart, ecType* edgecut, MLGP_option opt)
{
    if (opt.nbPart != 2)
        u_errexit("Refinement refinementPostOrder_Swap only works for nbpart = 2\n");

    idxType hsize[2] = {0,0};
    idxType i, maxpartsize = 0, nbmoved = 0, inPartIdx = toporderpart[0], outPartIdx = toporderpart[1];
    vwType maxvw = -1;
    int without_impro = 0, minindex = 0;
    int secondNodes[2]={-1,-1};
    idxType node[2]={-1,-1};
    inFromOutToCnts *ioCnts = (inFromOutToCnts*) malloc(sizeof(inFromOutToCnts) * (G->nVrtx+1));
    idxType* partcopy = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    idxType* moved = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    ecType* gain = (ecType*) malloc(sizeof(ecType)*(G->nVrtx+1));
    idxType* heap[2];
    idxType* inheap = (idxType*) calloc(G->nVrtx+1, sizeof(idxType));
    int* anchored = (int*) calloc(G->nVrtx+1, sizeof(int));
    idxType partsize[2] = {0,0};
    heap[inPartIdx] = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    heap[outPartIdx] = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));

    bisInitializeGainAndHeap(G, inPartIdx, outPartIdx, ioCnts, gain, heap, inheap, hsize, part, partsize, &maxpartsize, &maxvw, partcopy, opt);

    idxType minmaxpartsize = maxpartsize;
    ecType minedgecut = *edgecut;

    while ((hsize[inPartIdx] > 0 && hsize[outPartIdx] > 0)&&(without_impro <= G->nVrtx/4)&&(nbmoved<G->nVrtx)){
        node[inPartIdx] = heapExtractMax(heap[inPartIdx], gain, &hsize[inPartIdx], inheap);
        node[outPartIdx] = heapExtractMax( heap[outPartIdx], gain, &hsize[outPartIdx], inheap);

        // Check cycle and constraint for top pair to get flagger
        int flagger = getFlagger(node, secondNodes, inPartIdx, outPartIdx, G, part, partsize, nbpart, partmatrix, lb_pw, ub_pw, heap, hsize, inheap, gain);

        switch (flagger){
            case 0:
                continue;
                break;

            case 1: //We move node[inPartIdx] and node[outPartIdx]
                break;

            case 2: //We mode node[inPartIdx] and secondNodes[outPartId]
                if(secondNodes[inPartIdx]!=-1)
                    heapInsert(heap[part[secondNodes[inPartIdx]]], gain, &hsize[part[secondNodes[inPartIdx]]], secondNodes[inPartIdx], inheap);
                node[outPartIdx]=secondNodes[outPartIdx];
                break;

            case 3: //We mode secondNodes[inPartIdx] and nodes[outPartIdx]
                if(secondNodes[outPartIdx]!=-1)
                    heapInsert(heap[part[secondNodes[outPartIdx]]], gain, &hsize[part[secondNodes[outPartIdx]]], secondNodes[outPartIdx], inheap);
                node[inPartIdx]=secondNodes[inPartIdx];
                break;

            case 4: //We mode secondNodes[inPartIdx] and secondNodes[outPartIdx]
                node[inPartIdx]=secondNodes[inPartIdx];
                node[outPartIdx]=secondNodes[outPartIdx];
                break;

            default:
                u_errexit("In refinementPostOrder_Swap: wrong value for flagger %d.\n", flagger);
        }

        computePartmatrix(G, part, partmatrix, node[inPartIdx], 1-part[node[inPartIdx]]);
        computePartmatrix(G, part, partmatrix, node[outPartIdx], 1-part[node[outPartIdx]]);

        //Moves the nodes
        bisMoveAndUpdate(G, node[inPartIdx], moved, &nbmoved, edgecut, inPartIdx, outPartIdx, ioCnts, gain, part, partsize, &maxpartsize, heap, hsize, inheap, anchored, opt);
        bisMoveAndUpdate(G, node[outPartIdx], moved, &nbmoved, edgecut, inPartIdx, outPartIdx, ioCnts, gain, part, partsize, &maxpartsize, heap, hsize, inheap, anchored, opt);

        if ((*edgecut < minedgecut)||((*edgecut == minedgecut)&&(maxpartsize <= minmaxpartsize))){
            minedgecut = *edgecut;
            minmaxpartsize = maxpartsize;
            minindex = nbmoved;
            without_impro = -1;
        }
        without_impro ++;
    }

    for (i = nbmoved-1; i >= minindex; i--){
        part[moved[i]] = partcopy[moved[i]];
    }

    *edgecut = minedgecut;
    free(partcopy);
    free(moved);
    free(gain);
    free(heap[0]);
    free(heap[1]);
    free(inheap);
    free(anchored);
}


void refinementPostOrder_Max(dgraph *G, idxType *toporderpart,
                             double* lb_pw, double* ub_pw,
                             idxType *part, ecType* edgecut, MLGP_option opt, idxType *partsize)
{
    if (opt.nbPart != 2)
        u_errexit("RefinementPostOrder_Max only works for nbpart = 2\n");

    idxType i, without_impro = 0, minindex = 0, maxpartsize = 0, nbmoved = 0, newsize1, newsize2, nodeToMove = -1;
    int balNotOk,   biggerPart = -1, other, inPartIdx = toporderpart[0], outPartIdx = toporderpart[1];
    vwType maxvw = -1;
    idxType hsize[2] = {0,0};
    inFromOutToCnts *ioCnts = (inFromOutToCnts*) malloc(sizeof(inFromOutToCnts) * (G->nVrtx+1));
    idxType* partcopy = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    idxType* moved = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    ecType* gain = (ecType*) calloc(G->nVrtx+1, sizeof(ecType));
    idxType* heap[2];
    heap[0] = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    heap[1] = (idxType*) malloc(sizeof(idxType)*(G->nVrtx+1));
    idxType* inheap = (idxType*) calloc(G->nVrtx+1, sizeof(idxType));
    int* anchored = (int*) calloc(G->nVrtx+1, sizeof(int));

    bisInitializeGainAndHeap(G, inPartIdx, outPartIdx, ioCnts, gain, heap, inheap, hsize, part, partsize, &maxpartsize, &maxvw, partcopy, opt);

    balNotOk = partSizeChecker(partsize, ub_pw, 2, maxvw) ;
    if(opt.debug > 10)
    {
        printf("\t\trefintement: initial cut %.0f bal %d\n", (double)*edgecut, balNotOk);

        heapVerify(heap[inPartIdx], gain, hsize[inPartIdx], inheap);
        heapVerify(heap[outPartIdx], gain, hsize[outPartIdx], inheap);
        printf("refinementPostOrder_Max: verified heaps %.0f %.0f at the beginning\n", (double) hsize[inPartIdx], (double) hsize[outPartIdx]);
    }

    idxType minmaxpartsize = maxpartsize;
    ecType minedgecut = *edgecut;

    while ((hsize[inPartIdx] > 0 || hsize[outPartIdx] > 0)&&(without_impro <= G->nVrtx/4)&&(nbmoved<G->nVrtx)){
        biggerPart = getBiggerPart(inPartIdx, outPartIdx, partsize, heap, hsize, lb_pw, ub_pw, maxvw, gain, balNotOk, opt);
        if (biggerPart == -1)
            break;

        other = 1 - biggerPart;

        nodeToMove = heapExtractMax(heap[biggerPart], gain, &hsize[biggerPart], inheap);

        anchored[nodeToMove] = 1;

        //check new sizes to see if we can move one node from the first part to second one
        newsize1 = partsize[biggerPart] - G->vw[nodeToMove];
        newsize2 = partsize[other] + G->vw[nodeToMove];

        if ((newsize1 < lb_pw[biggerPart]-maxvw) || (newsize2 > ub_pw[other]+maxvw))
            continue;

        //Moves the node
        bisMoveAndUpdate(G, nodeToMove, moved, &nbmoved, edgecut, inPartIdx, outPartIdx, ioCnts, gain, part, partsize, &maxpartsize, heap, hsize, inheap, anchored, opt);

        if (balNotOk || ((*edgecut < minedgecut)||((*edgecut == minedgecut)&&(maxpartsize <= minmaxpartsize)))){
            minedgecut = *edgecut;
            minmaxpartsize = maxpartsize;
            minindex = nbmoved;
            without_impro = -1;
            balNotOk = partSizeChecker(partsize, ub_pw, 2, maxvw) ;
        }
        if(balNotOk == 0)
            without_impro ++;
        if (opt.debug > 10)
        {
            printf("refinementPostOrder_Max: verifying heaps at %d %d %d (Pin %.0f Pout %.0f)\n", nbmoved, hsize[inPartIdx], hsize[outPartIdx], (double) partsize[inPartIdx], (double)partsize[outPartIdx]);
            heapVerify(heap[inPartIdx], gain, hsize[inPartIdx], inheap);
            heapVerify(heap[outPartIdx], gain, hsize[outPartIdx], inheap);
        }

    }
    for (i = nbmoved-1; i >= minindex; i--){
        part[moved[i]] = partcopy[moved[i]];
    }

    if (opt.debug > PD_HIGH)
        printf("\t\trefinementPostOrder_Max: Nb moved %d (minindex %d)\n",nbmoved, minindex);

    *edgecut = minedgecut;
    free(ioCnts);
    free(partcopy);
    free(moved);
    free(gain);
    free(heap[inPartIdx]);
    free(heap[outPartIdx]);
    free(inheap);
    free(anchored);
}

void recBisRefinementStep(MLGP_option opt, int level, dgraph* graph, idxType* part, idxType nbpart, MLGP_info* info)
{
    if (opt.refinement == REF_NONE)
        return;

    idxType j, i, maxpartsize = 0, previous_maxpartsize;
    // if (opt.nbPart != nbpart)
    //     u_errexit("nope");
    vwType* partsize = (vwType*) calloc(opt.nbPart+1, sizeof(vwType));
    if (partsize == NULL)
        u_errexit("Error in allocation: partsize.\n");
    int nbbalance = 1;
    double sum = 0;
    double* bound = (double*) umalloc(sizeof(double) * nbpart, "recBisRefinementStep: bound");
    idxType* toporder = (idxType*) umalloc(sizeof(idxType) * nbpart, "recBisRefinementStep: toporder");
    ecType lastobjvalue, objvalue = -1;
    ecType** partmatrix=NULL;

    switch(opt.refinement) {
        case REF_kFM_TB :
        case REF_kFM_POTB :
        case REF_KL :
        case REF_KL_bFM_MAXW :
        case REF_kFM_RANDV:
            partmatrix = (ecType **) umalloc(opt.nbPart * sizeof(ecType * ), "recBisRefinementStep: partmatrix");
            for (j = 0; j < opt.nbPart; j++)
                partmatrix[j] = (ecType *) umalloc(opt.nbPart * sizeof(ecType), "recBisRefinementStep: partmatrix[j]");
            break;
    }

    for (i=1; i<=graph->nVrtx; i++)
        partsize[part[i]] += graph->vw[i];

    for (i = 0; i < nbpart; i++) {
        sum += opt.ub[i];
    }

    for (i = 0; i < nbpart; i++) {
        bound[i] = graph->totvw * opt.ub[i] / sum + 0.9 * (opt.ub[i] - graph->totvw * opt.ub[i] / sum);
        maxpartsize = maxpartsize < partsize[i] ? partsize[i] : maxpartsize;
    }

    info->timing_forced_tab[level] -= u_wseconds();

    randTopSortOnParts(graph, part, toporder, nbpart);

    objvalue = edgeCut(graph, part);

    while ((nbbalance < 2) && partSizeChecker(partsize,bound,nbpart, graph->maxVW)){
        previous_maxpartsize = maxpartsize;

        int oldRef = opt.refinement;
        opt.refinement = REF_bFM_MAXW; //Temporary setting refinement to be sure to use FM bisection from Max
        refinementPostOrder_Max(graph, toporder, opt.lb, bound, part, &objvalue, opt, partsize);
        info->current_edgecut = objvalue;
        opt.refinement  = oldRef;

        nbbalance++;
        randTopSortOnParts(graph, part, toporder, nbpart);
    }

    lastobjvalue = objvalue;

    info->ec_afterforced_tab[level] = info->current_edgecut;
    info->timing_forced_tab[level] += u_wseconds();

    if (opt.debug > 0){
        printf("We ran %d Forced Balance iterations\n", nbbalance);
        for(i=0;i<nbpart;++i)
            printf("%d, ",partsize[i] );
        printf("PARTSIZES\n");
    }

    if(objvalue < 0)
        u_errexit("The objective value has not been set correctly\n");

    info->nbref_tab[level] = 0;

    for (i=1; i<= opt.ref_step; i++) {
        if (opt.debug > 10){
            print_graph(graph, part, opt.file_name_dot, level, i);
        }
        info->timing_refinement_tab_tab[level][i] -= u_wseconds();

        switch(opt.refinement) {
            case REF_bFM:
            case REF_bFM_MAXW :
                refinementPostOrder_Max(graph, toporder, opt.lb, opt.ub, part, &objvalue, opt, partsize);
                break;

            case REF_KL :
                buildPartmatrix(graph, part, opt.nbPart, partmatrix);
                refinementPostOrder_Swap(graph, partmatrix, toporder, opt.lb, opt.ub, part, nbpart, &objvalue, opt);
                break;

            case REF_KL_bFM_MAXW :
                refinementPostOrder_Max(graph, toporder, opt.lb, opt.ub, part, &objvalue, opt, partsize);
                info->current_edgecut = objvalue;
                buildPartmatrix(graph, part, opt.nbPart, partmatrix);
                refinementPostOrder_Swap(graph, partmatrix, toporder, opt.lb, opt.ub, part, nbpart, &objvalue, opt);
                break;

            default:
                u_errexit("Refinement %d is not recognized!\n", opt.refinement);
        }

        info->current_edgecut = objvalue;
        info->timing_refinement_tab_tab[level][i] += u_wseconds();
        info->ec_afterref_tab[level][i] = info->current_edgecut;
        info->nbref_tab[level] = i;

        if(opt.debug>4) {
            printf("\tVerify edge cut for coars level %d refinement step %d (stored ec %.2f)\n", level, i, 1.0*info->current_edgecut);
            ecType cutComputed = edgeCut(graph, part);
            if(cutComputed !=  info->current_edgecut )
                u_errexit("edge cut computer %.2f but stored %.2f\n", 1.0*cutComputed, 1.0*info->current_edgecut );
        }
        if (opt.debug > 0)
            printf("   ref %d\t%d\t%d\t%d\t%d\n", i, graph->nVrtx, graph->nEdge, (int) objvalue, 0);

        if (objvalue >= 0.99*lastobjvalue)
            break;

        lastobjvalue = objvalue;
    }

    switch(opt.refinement) {
        case REF_kFM_TB :
        case REF_kFM_POTB :
        case REF_KL :
        case REF_KL_bFM_MAXW :
        case REF_kFM_RANDV:
            for (j = 0; j < opt.nbPart; j++)
                free(partmatrix[j]);
            free(partmatrix);
            break;
    }

    free(toporder);
    free(bound);
    free(partsize);
}