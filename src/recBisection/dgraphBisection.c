#include "dgraphBisection.h"

void fixAcyclicityBottom(dgraph* G, idxType *part, idxType outpartidx)
{
    /*Fix acyclicity of a bisection where inpartidx -> outpartidx by puting every child of outpartidx to outpartidx*/
    idxType to = 1, i;
    idxType* ready = (idxType*) calloc((G->nVrtx + 1), sizeof(idxType));
    idxType nbready = sourcesList(G, ready);
    idxType* nbinleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    for (i = 1; i <= G->nVrtx; i++)
        nbinleft[i] = G->inEnd[i] - G->inStart[i] + 1;
    while (nbready > 0) {

        idxType node = ready[nbready-1];
        nbready--;
        to = to + 1;
        if (part[node] == outpartidx)
            for (i = G->outStart[node]; i <= G->outEnd[node]; i++)
                part[G->out[i]] = outpartidx;
        for (i = G->outStart[node]; i <= G->outEnd[node]; i++) {
            idxType succ = G->out[i];
            nbinleft[succ]--;
            if (nbinleft[succ] == 0)
                ready[nbready++] = succ;
            else if (nbinleft[succ]<0)
                u_errexit("fixAcyclicityBottom: negative indegree\n");
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("fixAcyclicityBottom Problem, not every node concerned!\n");
    free(nbinleft);
    free(ready);
}


void fixAcyclicityUp(dgraph* G, idxType *part, idxType inpartidx)
{
    /*Fix acyclicity of a bisection where inpartidx -> outpartidx by puting every child of outpartidx to outpartidx*/
    idxType* ready;
    idxType to = 1, i;
    ready = (idxType*) umalloc(sizeof(idxType) * (G->nVrtx+1), "ready");
    idxType nbready = outputList(G, ready);
    idxType* nboutleft = (idxType *) calloc(G->nVrtx+1, sizeof(idxType));
    for (i = 1; i <= G->nVrtx; i++)
        nboutleft[i] = G->outEnd[i] - G->outStart[i] + 1;
    while (nbready > 0) {

        idxType node = ready[nbready-1];
        nbready--;
        to++;
        if (part[node] == inpartidx)
            for (i = G->inStart[node]; i <= G->inEnd[node]; i++)
                part[G->in[i]] = inpartidx;
        for (i = G->inStart[node]; i <= G->inEnd[node]; i++) {
            idxType pred = G->in[i];
            nboutleft[pred]--;
            if (nboutleft[pred] == 0)
                ready[nbready++] = pred;
            else if (nboutleft[pred]<0)
                u_errexit("fixAcyclicityUp: negative indegree\n");
        }
    }
    if (to != G->nVrtx+1)
        u_errexit("fixAcyclicityUp: Problem, not every node concerned! to = %d != %d = nVrtx\n", (int) to, (int) G->nVrtx);
    free(nboutleft);
    free(ready);
}



void fixUndirBisection(ecType* edgecut, coarsen* coars, idxType* part, MLGP_option opt, MLGP_info* info) {
    dgraph* G = coars->graph;
    idxType* part_up = (idxType*) calloc(G->nVrtx+1, sizeof(idxType));
    idxType* part_down = (idxType*) calloc(G->nVrtx+1, sizeof(idxType));
    ecType  ecfix;
    idxType inpartidx, outpartidx;
    idxType i;
    ecType obj_up, obj_down, new_best_obj = ecType_MAX, best_obj = ecType_MAX;
    vwType* partsize = (vwType*) umalloc((opt.nbPart)*sizeof(vwType), "partsize");

    for (inpartidx = 0 ;inpartidx <=  1; inpartidx++) {
        memcpy(part_up, part, sizeof(idxType)*(G->nVrtx+1));
        memcpy(part_down, part, sizeof(idxType)*(G->nVrtx+1));

        outpartidx = 1-inpartidx;
        idxType toporderpart[2];

        toporderpart[0] = inpartidx;
        toporderpart[1] = outpartidx;

        // middown  = ecType_MAX;
        // midup  = ecType_MAX;
        //We try fixAcyclicityBottom with the order inpartidx -> outpartidx
        fixAcyclicityBottom(G, part_down, toporderpart[1]);

        obj_down = edgeCut(G, part_down);
        // middown = obj_down;
        int oldRef = opt.refinement;
        opt.refinement = REF_bFM_MAXW;
        refinementPostOrder_Max(G, toporderpart, opt.lb, opt.ub, part_down, &obj_down, opt, partsize);
        opt.refinement  = oldRef;
        recBisRefinementStep(opt, info->coars_depth, G, part_down, coars->nbpart, info);
        obj_down = edgeCut(G, part_down);
        if ( !(partsize[0] <= opt.ub[0] + G->maxVW && partsize[1] <= opt.ub[1] + G->maxVW))
            obj_down = ecType_MAX;
        //We try fixAcyclicityUp with the order inpartidx -> outpartidx
        fixAcyclicityUp(G, part_up, toporderpart[0]);
        obj_up = edgeCut(G, part_up);
        // midup = obj_up;
        oldRef = opt.refinement;
        opt.refinement = REF_bFM_MAXW;
        refinementPostOrder_Max(G, toporderpart, opt.lb, opt.ub, part_up, &obj_up, opt, partsize);
        opt.refinement  = oldRef;
        recBisRefinementStep(opt, info->coars_depth, G, part_up, coars->nbpart, info);
        obj_up = edgeCut(G, part_up);
        if ( !(partsize[0] <= opt.ub[0] + G->maxVW && partsize[1] <= opt.ub[1] + G->maxVW))
            obj_up = ecType_MAX;

        if (opt.debug > PD_NONE)
            printf("cuts (undir) FixUp: %.0f FixDown: %.0f\n", (double) obj_up,  (double) obj_down);
        new_best_obj = new_best_obj < obj_up ? new_best_obj : obj_up;
        new_best_obj = new_best_obj < obj_down ? new_best_obj : obj_down;

        if (new_best_obj >= best_obj)
            continue;

        if (obj_up > obj_down)
            memcpy(part, part_down, sizeof(idxType)*(G->nVrtx+1));
        else{
            if (obj_up < ecType_MAX)
                memcpy(part, part_up, sizeof(idxType)*(G->nVrtx+1));
            else
                u_errexit("undirAcyclicBisection: This should not have happened\n");
        }
        *edgecut = obj_up < obj_down ? obj_up : obj_down;
        best_obj = *edgecut;
    }
    if (opt.print >=PD_LOW) {
        printf("\t\tInitial+Fix+refinement\t%d\n",(int) *edgecut);
    }
    free(part_up);
    free(part_down);
    free(partsize);
}

void greedyGraphGrowingEqualTopOrder(dgraph *G, double* ub_pw, idxType *part, MLGP_option opt, int fromPart, int otherPart, ecType *edgecut)
{
/* This is a place-holder for the following heuristic for bisection. It seems
 * equal to greedy graph growing. We already have its equivalent for unit weighted
 * graphs. This will be for vertex/edge weighted graphs.
 *
 * We will build the part 0, assuming all vertices in the part 1 at the outset.
 * Heap contains only indegree-0 vertices. A vertex is indegree 0, if there is no
 * incoming edge from a vertex belonging to the part 1.
 * Among those vertices pick the one with the highest
 * (static) key value sum_inedgecost - sum_outedgecost.
 *
 * This key-values are static because we can put a vertex from the part 1 to 0
 * only if all its incoming vertices are in 0 already.
 *
 * So no key update. The gain=cost(incoming edges)-cost(outgoing edges).
 *
 *
 * We will get the best edgecut value from the moment the part0 has enough weights
 * to the moment the part 1 has no more enough weight.
 */
    idxType i, j, vi, neigh, nVrtx = G->nVrtx;
    idxType  *inEnd = G->inEnd,  *inStart = G->inStart;
    idxType *outEnd = G->outEnd, *outStart = G->outStart, *out = G->out;
    ecType  *ecIn = G->ecIn, *ecOut = G->ecOut;

    vwType maxvw;
    idxType *inDegrees;
    idxType heapsz = 0; /*will be from 1 to heapsz*/
    idxType maxind;

    idxType *heapIds, *placedinZero, moved, bestAt, *visitOrder;
    ecType *inMinusOut; /*this will be keyval in the heap*/
    ecType cut, bestCut=ecType_MAX;

    vwType *vw = G->vw;
    double pw[2], bestCutBal;

    inMinusOut = (ecType *) umalloc(sizeof(ecType) * (nVrtx+1), "inMinusOut");
    inDegrees = (idxType*) umalloc(sizeof(idxType) * (nVrtx+1), "inDegrees");
    heapIds = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "heapIds");
    placedinZero = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "placedinZero");
    visitOrder = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "visitOrder");
    int isValidBsctn = 0;
    shuffleTab(1, nVrtx, visitOrder+1);/*this is the randomness */

    maxvw = -1;
    heapIds[0] = -1;
    for (vi = 1; vi <= nVrtx; vi++)
    {
        ecType in_out = 0;
        i = visitOrder[vi];

        part[i] = fromPart;

        if (vw[i] > maxvw) maxvw = vw[i];

        inDegrees[i] = inEnd[i] - inStart[i] + 1;
        if (inDegrees[i] == 0)
            heapIds[++heapsz] = i;

        for (j=inStart[i]; j <= inEnd[i]; j++)
            in_out += ecIn[j];

        for (j=outStart[i]; j <= outEnd[i]; j++)
            in_out -= ecOut[j];

        inMinusOut[i] = in_out;
    }

    maxBuildHeapKeyvals(heapIds, heapsz, inMinusOut, vw);

    pw[otherPart] = 0;
    pw[fromPart] = bestCutBal = G->totvw;
    bestAt = -1;
    moved = 0;
    cut = 0;

    while (pw[otherPart] <= ub_pw[otherPart] + maxvw && pw[fromPart]>=opt.lb[fromPart] - maxvw)
    {
        if (heapsz < 1)
            break;
        maxind = maxHeapExtractKeyvals(heapIds, &heapsz, inMinusOut, vw);
        if (part[maxind] == otherPart)
            u_errexit("greedGraphGrowingEqualTopOrder: already moved there?\n");

        part[maxind] = otherPart;
        pw[otherPart] += vw[maxind];
        pw[fromPart] -= vw[maxind];

        placedinZero[moved++] = maxind;
        cut -= inMinusOut[maxind];
        for (j = outStart[maxind]; j<= outEnd[maxind]; j++)
        {
            neigh = out[j];
            inDegrees[neigh] --;
            if (inDegrees[neigh] == 0)
                maxHeapInsertKeyvals(heapIds, &heapsz, neigh, inMinusOut, vw);
            else if (inDegrees[neigh] == -1)
                u_errexit("wrong index\n");
        }
        double maxpb = mymax(pw[0]/ub_pw[0], pw[1]/ub_pw[1]);
        if (bestAt == -1 && maxpb < bestCutBal) {
            bestCut = cut;
            bestAt = moved - 1;
            bestCutBal = maxpb;
        }
        if (pw[otherPart] <= ub_pw[otherPart]+maxvw && pw[fromPart] <= ub_pw[fromPart]+ maxvw) /*if both are less than their upper bound, this is a valid partition*/
        {
            if ( (isValidBsctn == 0) || (cut < bestCut  || (cut == bestCut && maxpb < bestCutBal)))
            {
                bestCut = cut;
                bestAt = moved - 1;
                bestCutBal = maxpb;
                isValidBsctn = 1;
            }
        }
    }
    if (bestAt == -1)
        u_errexit("greedGraphGrowingEqualTopOrder: cannot find any balanced bipartition\n");
    for (i = bestAt+1; i < moved; i++)
        part[placedinZero[i]] = fromPart;

    *edgecut = bestCut;
    free(visitOrder);
    free(placedinZero);
    free(heapIds);
    free(inDegrees);
    free(inMinusOut);
}

void greedyGraphGrowingTwoPhases(dgraph *G, double* ub_pw, idxType *part, MLGP_option opt, int fromPart, int otherPart, ecType *edgecut)
{
/* This is a place-holder for the following heuristic for bisection.
 */
    idxType i, j, vi, neigh, nVrtx = G->nVrtx;
    idxType  *inEnd = G->inEnd,  *inStart = G->inStart;
    idxType *outEnd = G->outEnd, *outStart = G->outStart, *out = G->out;
    ecType  *ecIn = G->ecIn, *ecOut = G->ecOut;

    vwType maxvw;
    idxType *inDegrees;
    idxType heapsz = 0; /*will be from 1 to heapsz*/
    idxType maxind;
    int isValidBsctn = 0;
    idxType *heapIds, *placedinZero, moved, bestAt, *visitOrder;
    ecType *inMinusOut, *inSum; /*this will be keyval in the heap*/
    ecType cut, bestCut=ecType_MAX;

    vwType *vw = G->vw;
    double pw[2], bestCutBal;

    inMinusOut = (ecType *) umalloc(sizeof(ecType) * (nVrtx+1), "inMinusOut");
    inDegrees = (idxType*) umalloc(sizeof(idxType) * (nVrtx+1), "inDegrees");
    heapIds = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "heapIds");
    placedinZero = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "placedinZero");
    visitOrder = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "visitOrder");
    inSum = (ecType *) umalloc(sizeof(ecType) * (nVrtx+1), "inSum");

    shuffleTab(1, nVrtx, visitOrder+1);/*this is the randomness */

    maxvw = -1;
    heapIds[0] = -1;
    for (vi = 1; vi <= nVrtx; vi++) {
        ecType in_out = 0;
        i = visitOrder[vi];

        // assign everyhing to fromPart
        part[i] = fromPart;

        // find maxvw
        if (vw[i] > maxvw) maxvw = vw[i];

        // if this is a source node, add it to heap
        inDegrees[i] = inEnd[i] - inStart[i] + 1;
        if (inDegrees[i] == 0)
            heapIds[++heapsz] = i;

        // find incoming edge weight and the difference between in and out: -gain
        //      (all outs should be in the succeeding part, all ins should be in the preceding part
        //      if moving from succeeding part to preceding part)
        for (j=inStart[i]; j <= inEnd[i]; j++)
            in_out += ecIn[j];

        inSum[i] = in_out;

        for (j=outStart[i]; j <= outEnd[i]; j++)
            in_out -= ecOut[j];

        inMinusOut[i] = in_out;
    }

    // part weights are assigned wrt the loop above
    pw[otherPart] = 0;
    pw[fromPart] = bestCutBal = G->totvw;
    bestAt = -1;
    moved = 0;
    cut = 0;

    idxType anchorSource = heapIds[1];
    idxType* distFromAnchor = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "idxType");
    computeDistances(G, anchorSource, distFromAnchor);
    for (i = 1; i <= G->nVrtx; i++)
        distFromAnchor[i] = -distFromAnchor[i];

    maxBuildHeapKeyvals(heapIds, heapsz, inSum, distFromAnchor);
    //First Phase where we fill based on in only (not gain) and distance to anchor
    while (pw[otherPart] <= 0.9*(ub_pw[otherPart] + maxvw) && pw[fromPart]>=1.1*(opt.lb[fromPart] - maxvw))
    {
        if (heapsz < 1)
            break;
        maxind = maxHeapExtractKeyvals(heapIds, &heapsz, inSum, distFromAnchor);
        if (part[maxind] == otherPart)
            u_errexit("greedGraphGrowingTwoPhases: already moved there?\n");

        part[maxind] = otherPart;
        pw[otherPart] += vw[maxind];
        pw[fromPart] -= vw[maxind];

        placedinZero[moved++] = maxind;
        cut -= inMinusOut[maxind];
        for (j = outStart[maxind]; j<= outEnd[maxind]; j++)
        {
            neigh = out[j];
            inDegrees[neigh] --;
            if (inDegrees[neigh] == 0)
                maxHeapInsertKeyvals(heapIds, &heapsz, neigh, inSum, distFromAnchor);
            else if (inDegrees[neigh] == -1)
                u_errexit("wrong index\n");
        }
        double maxpb = mymax(pw[0]/ub_pw[0], pw[1]/ub_pw[1]);
        if (bestAt == -1 && maxpb < bestCutBal) {
            bestCut = cut;
            bestAt = moved - 1;
            bestCutBal = maxpb;
        }
        // if both are less than their upper bound, this is a valid partition
        if (pw[otherPart] <= ub_pw[otherPart] + maxvw && pw[fromPart] <= ub_pw[fromPart] + maxvw)
        {
            if ( (isValidBsctn == 0) || (cut < bestCut  || (cut == bestCut && maxpb < bestCutBal)))
            { // if there is a valid bisection, regular comparison. if no valid bisection exists, then update with this.
                bestCut = cut;
                bestAt = moved - 1;
                bestCutBal = maxpb;
                isValidBsctn = 1;
            }
        }
    }

    maxBuildHeapKeyvals(heapIds, heapsz, inMinusOut, vw);

    // printf("second loop\n");

    //Second Phase where we fill based on inMinusOut and vertex size
    while (pw[otherPart] <= ub_pw[otherPart] + maxvw && pw[fromPart]>=opt.lb[fromPart] - maxvw)
    {
        // printf ("pw[fromPart=%d] = %lf , pw[otherPart = %d] = %lf, ub_pw %lf %lf lb %lf %lf maxvw %lf left %lf right %lf\n", fromPart, pw[fromPart], otherPart, pw[otherPart], ub_pw[fromPart], ub_pw[otherPart], opt.lb[fromPart], opt.lb[otherPart], (double) maxvw, 0.9*(ub_pw[otherPart] + maxvw),1.1*(opt.lb[fromPart] - maxvw));
        if (heapsz < 1)
            break;
        maxind = maxHeapExtractKeyvals(heapIds, &heapsz, inMinusOut, vw);
        if (part[maxind] == otherPart)
            u_errexit("greedGraphGrowingTwoPhases: already moved there?\n");

        part[maxind] = otherPart;
        pw[otherPart] += vw[maxind];
        pw[fromPart] -= vw[maxind];

        placedinZero[moved++] = maxind;
        cut -= inMinusOut[maxind];
        for (j = outStart[maxind]; j<= outEnd[maxind]; j++)
        {
            neigh = out[j];
            inDegrees[neigh] --;
            if (inDegrees[neigh] == 0)
                maxHeapInsertKeyvals(heapIds, &heapsz, neigh, inMinusOut, vw);
            else if (inDegrees[neigh] == -1)
                u_errexit("wrong index\n");
        }
        // printf ("cut %lf pw[0] = %lf ub[0] = %lf pw[1] = %lf ub[1] = %lf\t-->\t %lf\n", (double) cut, pw[0], ub_pw[0], pw[1], ub_pw[1], maxpb);
        if (pw[otherPart] <= ub_pw[otherPart]+maxvw && pw[fromPart] <= ub_pw[fromPart]+ maxvw) /*if both are less than their upper bound, this is a valid partition*/
        {
            double maxpb = mymax(pw[0]/ub_pw[0], pw[1]/ub_pw[1]);
            if ( (isValidBsctn == 0)|| (cut < bestCut  || (cut == bestCut && maxpb < bestCutBal )) )
            { // if there is a valid bisection, regular comparison. if no valid bisection exists, then update with this.
                bestCut = cut;
                bestAt = moved - 1;
                bestCutBal = maxpb;
                isValidBsctn = 1;
                // printf ("bestCut %lf, bestCutBal %lf\n", (double) bestCut, (double) bestCutBal);
            }
        }
    }

    if (bestAt == -1)
        u_errexit("greedGraphGrowingTwoPhases: cannot find any balanced bipartition\n");
    for (i = bestAt+1; i < moved; i++)
        part[placedinZero[i]] = fromPart;

    *edgecut = bestCut;
    free(visitOrder);
    free(placedinZero);
    free(heapIds);
    free(inDegrees);
    free(inMinusOut);
}


void greedyGraphGrowingIn(dgraph *G, double* ub_pw, idxType *part, MLGP_option opt, int fromPart, int otherPart, ecType *edgecut)
{
/* This is a place-holder for the following heuristic for bisection.
 */
    idxType i, j, vi, neigh, nVrtx = G->nVrtx;
    idxType  *inEnd = G->inEnd,  *inStart = G->inStart;
    idxType *outEnd = G->outEnd, *outStart = G->outStart, *out = G->out;
    ecType  *ecIn = G->ecIn, *ecOut = G->ecOut;

    vwType maxvw;
    idxType *inDegrees;
    idxType heapsz = 0; /*will be from 1 to heapsz*/
    idxType maxind;

    idxType *heapIds, *placedinZero, moved, bestAt, *visitOrder;
    ecType *inMinusOut, *inSum; /*this will be keyval in the heap*/
    ecType cut, bestCut=ecType_MAX;

    vwType *vw = G->vw;
    double pw[2], bestCutBal;
    int isValidBsctn = 0;
    inMinusOut = (ecType *) umalloc(sizeof(ecType) * (nVrtx+1), "inMinusOut");
    inDegrees = (idxType*) umalloc(sizeof(idxType) * (nVrtx+1), "inDegrees");
    heapIds = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "heapIds");
    placedinZero = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "placedinZero");
    visitOrder = (idxType *) umalloc(sizeof(idxType) * (nVrtx+1), "visitOrder");
    inSum = (ecType *) umalloc(sizeof(ecType) * (nVrtx+1), "inSum");

    shuffleTab(1, nVrtx, visitOrder+1);/*this is the randomness */

    maxvw = -1;
    heapIds[0] = -1;
    for (vi = 1; vi <= nVrtx; vi++)
    {
        ecType in_out = 0;
        i = visitOrder[vi];

        part[i] = fromPart;

        if (vw[i] > maxvw) maxvw = vw[i];

        inDegrees[i] = inEnd[i] - inStart[i] + 1;
        if (inDegrees[i] == 0)
            heapIds[++heapsz] = i;

        for (j=inStart[i]; j <= inEnd[i]; j++)
            in_out += ecIn[j];

        inSum[i] = in_out;

        for (j=outStart[i]; j <= outEnd[i]; j++)
            in_out -= ecOut[j];

        inMinusOut[i] = in_out;
    }

    pw[otherPart] = 0;
    pw[fromPart] = bestCutBal = G->totvw;
    bestAt = -1;
    moved = 0;
    cut = 0;

    idxType anchorSource = heapIds[0];
    idxType* distFromAnchor = (idxType *) malloc(sizeof(idxType) * (nVrtx+1));
    computeDistances(G, anchorSource, distFromAnchor);
    for (i = 1; i <= G->nVrtx; i++)
        distFromAnchor[i] = -distFromAnchor[i];

    maxBuildHeapKeyvals(heapIds, heapsz, inSum, distFromAnchor);

    while (pw[otherPart] <= ub_pw[otherPart] + maxvw && pw[fromPart]>=opt.lb[fromPart] - maxvw)
    {
        if (heapsz < 1)
            break;
        maxind = maxHeapExtractKeyvals(heapIds, &heapsz, inSum, distFromAnchor);
        if (part[maxind] == otherPart)
            u_errexit("greedGraphGrowingIn: already moved there?\n");

        part[maxind] = otherPart;
        pw[otherPart] += vw[maxind];
        pw[fromPart] -= vw[maxind];

        placedinZero[moved++] = maxind;
        cut -= inMinusOut[maxind];
        for (j = outStart[maxind]; j<= outEnd[maxind]; j++)
        {
            neigh = out[j];
            inDegrees[neigh] --;
            if (inDegrees[neigh] == 0)
                maxHeapInsertKeyvals(heapIds, &heapsz, neigh, inSum, distFromAnchor);
            else if (inDegrees[neigh] == -1)
                u_errexit("wrong index\n");
        }
        double maxpb = mymax(pw[0]/ub_pw[0], pw[1]/ub_pw[1]);
        if (bestAt == -1 && maxpb < bestCutBal) {
            bestCut = cut;
            bestAt = moved - 1;
            bestCutBal = maxpb;
        }
        if (pw[otherPart] <= ub_pw[otherPart]+maxvw && pw[fromPart] <= ub_pw[fromPart]+ maxvw) /*if both are less than their upper bound, this is a valid partition*/
        {
            if ( (isValidBsctn == 0) || (cut < bestCut  || (cut == bestCut && maxpb < bestCutBal)))
            {
                bestCut = cut;
                bestAt = moved - 1;
                bestCutBal = maxpb;
                isValidBsctn = 1;
            }
        }
    }

    if (bestAt == -1)
        u_errexit("greedGraphGrowingIn: cannot find any balanced bipartition\n");
    for (i = bestAt+1; i < moved; i++)
        part[placedinZero[i]] = fromPart;

    *edgecut = bestCut;
    free(visitOrder);
    free(placedinZero);
    free(heapIds);
    free(inDegrees);
    free(inMinusOut);
}

