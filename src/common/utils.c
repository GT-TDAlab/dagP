#include "utils.h"
#include <sys/time.h>
#define U_MB     (1024*1024)

#define KEYVAL(i) (G->inEnd[i]-G->inStart[i] - G->outEnd[i] + G->outStart[i] )


/* Random generator: not thread safe */
static int _u_seed=847449649;

int usRandom(int s)
{
    _u_seed = s;
    srand(_u_seed);
    return  _u_seed;
}

int uGetSeed(void)
{
    return _u_seed;
}


int uRandom(int l)
{
    _u_seed = rand();
    return _u_seed % l;
}

void* umalloc(long size, char* msg)
{
    void* ptr = NULL;

    if (size <= 0)
        return NULL;

    ptr = (void*) malloc(size);
    if (!ptr) {
        u_errexit("Memory allocation failed for '%s'.\nRequested size: %ld byte(s) [%.2lf MB].\n", msg, size, (double)size/(double)U_MB);
    }

    return ptr;
}



void shuffleArray(idxType size, idxType* toBeShuffled, idxType cnt){
    idxType i,j,tmp;
    while(cnt>0){
        j = uRandom(size);
        i = uRandom(size);
        tmp = toBeShuffled[i];
        toBeShuffled[i] = toBeShuffled[j];
        toBeShuffled[j] = tmp;
        --cnt;
    }
}

void reverseTab(idxType* tab, idxType size){
    /*Reverse tab (index starting at 1)*/
    idxType tmp;
    int i, half = (int) size / 2;
    for (i = 1; i <= half; i++){
        tmp = tab[i];
        tab[i] = tab[size - i + 1];
        tab[size - i + 1] = tmp;
    }
}


void heapIncreaseKey(idxType *heapIds, ecType* key, idxType sz, idxType id, idxType *inheap, ecType newKey)
{
    UNUSED(sz);
    idxType i = inheap[id]; /*location in heap*/
    key[id] = newKey;

    while ((i>>1)>0 && (key[id] >= key[heapIds[i>>1]]))/*maybe we prefer newly touched guys come close to the top of the heap. */
    {
        heapIds[i] = heapIds[i>>1];
        inheap[heapIds[i]] = i;
        i = i>>1;
    }
    heapIds[i] = id;
    inheap[id] = i;
}

void heapDecreaseKey(idxType *heapIds, ecType *key, idxType sz, idxType id, idxType *inheap, ecType newKey)
{/*this is heapify*/
    key[id] = newKey;
    heapify(heapIds, key, sz, inheap[id], inheap);
}

void heapDelete(idxType *heapIds, ecType *key, idxType *sz, idxType id, idxType *inheap)
{
    idxType i = inheap[id];
    heapIds[i] = heapIds[*sz];
    inheap[heapIds[i] ] = i;
    *sz = *sz - 1;
    inheap[id] = 0;
    if (i <= *sz)
    {
        if (i > 1 && key[heapIds[i]] >= key[heapIds[i>>1]])
            heapIncreaseKey(heapIds, key, *sz, heapIds[i], inheap, key[heapIds[i]]);
        else
            heapDecreaseKey(heapIds, key, *sz, heapIds[i], inheap, key[heapIds[i]]);
    }
}

void heapInsert(idxType *heapIds, ecType* key, idxType *sz, idxType id, idxType *inheap)
{
    idxType j, kv = key[id];
    *sz = *sz+1;
    heapIds[*sz] = id;
    j = *sz;
    if(inheap[id] != 0)
        u_errexit("heapInsert: inserting id %d but already in (index %d)\n", id, inheap[id]);
    while(j > 1 && key[heapIds[j/2]] < kv)
    {
        heapIds[j] = heapIds[j/2] ;
        inheap[heapIds[j/2]] = j;
        j = j/2;
    }
    heapIds[j] = id;
    inheap[id] = j;
}

void heapify(idxType *heapIds, ecType* key, idxType sz, idxType i, idxType *inheap)
{
    idxType largest, j, l,r, tmp;

    largest = j = i;
    while(j<=sz/2)
    {
        l = 2*j;
        r = 2*j + 1;

        if (key[heapIds[l]]>key[heapIds[j]])
            largest = l;
        else
            largest = j;

        if (r<=sz && key[heapIds[r]]>key[heapIds[largest]])
            largest = r;

        if (largest != j)
        {
            tmp = heapIds[largest];
            heapIds[largest] = heapIds[j];
            inheap[heapIds[j]] = largest;

            heapIds[j] = tmp;
            inheap[heapIds[j]] = j;
            j = largest;
        }
        else
            break;
    }
}

void heapUpdateKey(idxType *heapIds, idxType *inheap, ecType *key, idxType sz, idxType id, ecType newKey)
{
    /*PRECOND: assumes that id is already in heap*/
    if(key[id] > newKey)/*this is decrease key*/
        heapDecreaseKey(heapIds, key, sz, id, inheap, newKey);
    else if (key[id] < newKey)
        heapIncreaseKey(heapIds, key, sz, id, inheap, newKey);
}

void heapVerify(idxType *heapIds, ecType* key, idxType sz, idxType *inheap)
{
    idxType i;
    for (i = 1; i <= sz/2; i++)
    {
        idxType l, r;

        l = i << 1 ;
        r = l + 1;
        if (inheap[heapIds[i]] != i)
            u_errexit("heapVerify: location in heap is not correct i.\n");
        if(inheap[heapIds[l]] != l)
            u_errexit("heapVerify: location in heap is not correct l.\n");

        if(key[heapIds[i]] < key[heapIds[l]])
            u_errexit("heapVerify: not a maxheap l. %.0f %.0f (%d %d)\n", (double) key[heapIds[i]], (double) key[heapIds[l]], i, l);

        if (r <= sz)
        {
            if(inheap[heapIds[r]] != r)
                u_errexit("heapVerify: location in heap is not correct r.\n");
            if(key[heapIds[i]] < key[heapIds[r]])
                u_errexit("heapVerify: not a maxheap r.\n");
        }
    }
    for (i = sz/2+1; i <= sz; i++)
         if (inheap[heapIds[i]] != i)
            u_errexit("heapVerify: location in heap is not correct i.\n");

}

idxType heapExtractMax(idxType *heapIds, ecType* key, idxType *sz, idxType *inheap)
{
    idxType maxind ;
    if (*sz < 1)
        u_errexit("heap underflow\n");

    maxind = heapIds[1];
    heapIds[1] = heapIds[*sz];
    inheap[heapIds[1]] = 1;

    *sz = *sz - 1;
    inheap[maxind] = 0;

    heapify(heapIds, key, *sz, 1, inheap);
    return maxind;

}

void heapBuild(idxType *heapIds, ecType* key, idxType sz, idxType *inheap)
{
    idxType i;
    for (i=sz/2; i>=1; i--)
        heapify(heapIds, key, sz, i, inheap);
}

void randHeapInsert(idxType *heapIds, ecType* key, idxType *sz, idxType i, idxType *inheap)
{
    idxType j;
    ecType kv = key[i];
    *sz = *sz+1;
    heapIds[*sz] = i;
    j = *sz;
    if(inheap[i] != 0)
        u_errexit("randHeapInsert: inserting this but already in\n");

    //if the key is equal to the one we are comparing, flip a coin and decide the order.
    while(j > 1 && ((key[heapIds[j/2]] < kv) || ( kv == key[heapIds[j/2]] && uRandom(2)==1) ) )
    {
        heapIds[j] = heapIds[j/2] ;
        inheap[heapIds[j/2]] = j;
        j = j/2;
    }
    heapIds[j] = i;
    inheap[i] = j;
}

void randheapify(idxType *heapIds, ecType* key, idxType sz, idxType i, idxType *inheap)
{
    idxType largest, j, l,r, tmp;

    largest = j = i;
    while(j<=sz/2)
    {
        l = 2*j;
        r = 2*j + 1;

        if ((key[heapIds[l]]>key[heapIds[j]]) || (key[heapIds[l]]==key[heapIds[j]] && uRandom(2)==1))
            largest = l;
        else
            largest = j;

        if(r<=sz && ((key[heapIds[r]]>key[heapIds[largest]]) || (key[heapIds[l]]==key[heapIds[largest]] && uRandom(2)==1)))
            largest = r;

        if(largest != j)
        {
            tmp = heapIds[largest];
            heapIds[largest] = heapIds[j];
            inheap[heapIds[j]] = largest;

            heapIds[j] = tmp;
            inheap[heapIds[j]] = j;
            j = largest;
        }
        else
            break;
    }
}

idxType randHeapExtractMax(idxType *heapIds, ecType* key, idxType *sz, idxType *inheap)
{
    idxType maxind ;
    if( *sz < 1)
        u_errexit("heap underflow\n");

    maxind = heapIds[1];
    heapIds[1] = heapIds[*sz];
    inheap[heapIds[1]] = 1;

    *sz = *sz - 1;
    inheap[maxind] = 0;

    randheapify(heapIds, key, *sz, 1, inheap);
    return maxind;
}

void randHeapBuild(idxType *heapIds, ecType* key, idxType sz, idxType *inheap)
{
    idxType i;
    for (i=sz/2; i>=1; i--)
        randheapify(heapIds, key, sz, i, inheap);
}


void u_errexit(const char  * const f_str,...)
{
    va_list argp;

    fflush(stdout);
    fflush(stderr);
    fprintf(stderr, "\n****** Error ******\n");
    va_start(argp, f_str);
    vfprintf(stderr, f_str, argp);
    va_end(argp);

    fprintf(stderr,"*******************\n");
    fflush(stderr);

    printf("Error in Execution\n");
    exit(12);
}

double u_wseconds(void)
{
    struct timeval tp;

    gettimeofday(&tp, NULL);

    return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
}

void shuffleTab(idxType deb, idxType end, idxType* shuffle)
{
    /*Suppose shuffle is allocated
      Return shuffle filled with a shuffle list deb...end*/
    idxType i, j, tmp;
    for (i=deb; i<=end; i++)
        shuffle[i-deb] = i;
    for (i=deb; i<=end; i++) {
        j = uRandom(end-deb+1);
        tmp = shuffle[i-deb];
        shuffle[i-deb] = shuffle[j];
        shuffle[j] = tmp;
    }
}

void maxHeapInsertKeyvals(idxType *heapIds, idxType *sz, idxType i, ecType *keyvals, vwType *vw)
{
    idxType j;
    ecType kv = keyvals[i];

    *sz = *sz+1;
    heapIds[*sz] = i;
    j = *sz;

    while(j > 1 && (keyvals[heapIds[j/2]] < kv || (keyvals[heapIds[j/2]] == kv && vw[heapIds[j/2]]< vw[heapIds[j]])))
    {
        heapIds[j] = heapIds[j/2] ;
        j = j/2;
    }
    heapIds[j] = i;
}

void maxHeapifyKeyvals(idxType *heapIds, idxType sz, idxType i, ecType *keyvals, vwType *vw)
{
    idxType largest, j, l,r, tmp;

    largest = j = i;
    while(j<=sz/2)
    {
        l = 2*j;
        r = 2*j + 1;

        largest = j;
        if ((keyvals[heapIds[l]]>keyvals[heapIds[j]] || (keyvals[heapIds[l]]==keyvals[heapIds[j]] &&
                                                         vw[heapIds[l]]> vw[heapIds[j]])))
            largest = l;

        if (r<=sz && ( keyvals[heapIds[r]]>keyvals[heapIds[largest]] || (keyvals[heapIds[r]]==keyvals[heapIds[largest]] && vw[heapIds[r]]> vw[heapIds[largest]])))
            largest = r;
        if (largest != j)
        {
            tmp = heapIds[largest];
            heapIds[largest] = heapIds[j];
            heapIds[j] = tmp;
            j = largest;
        }
        else
            break;
    }
}

idxType maxHeapExtractKeyvals(idxType *heapIds, idxType *sz, ecType *keyvals, vwType *vw)
{
    idxType maxind ;
    if ( *sz < 1)
        u_errexit("maxHeapExtractKeyvals: heap underflow\n");

    maxind = heapIds[1];
    heapIds[1] = heapIds[*sz];
    *sz = *sz - 1;
    maxHeapifyKeyvals(heapIds, *sz, 1, keyvals, vw);
    return maxind;

}

void maxBuildHeapKeyvals(idxType *heapIds, idxType sz, ecType *keyvals, vwType *vw)
{
    idxType i;
    for (i=sz/2; i>=1; i--)
        maxHeapifyKeyvals(heapIds, sz, i, keyvals, vw);
}