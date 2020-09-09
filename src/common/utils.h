#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>

#include "dgraph.h"

#define mymax(a, b)	(((a) >= (b)) ? (a) : (b))
#define isPowOfTwo(x)  (((x)!=0) && !((x)&((x)-1)))


#define UNUSED(x) (void)(x)


/* our random generatos */

int usRandom(int s);
int uGetSeed(void);
int uRandom(int l);

void u_errexit(const char  * const f_str,...);
double u_wseconds(void);

void* umalloc(long size, char* msg);

void shuffleArray(idxType size, idxType* toBeShuffled, idxType cnt);

/* heaps */
void heapDelete(idxType *heapIds, ecType* key, idxType *sz, idxType id, idxType *inheap);
void heapUpdateKey(idxType *heapIds, idxType *inheap, ecType *key, idxType sz, idxType id, ecType newKey);

void randHeapRemove(idxType *heapIds, ecType* key, idxType *sz, idxType i,idxType *inheap);
void heapInsert(idxType *heapIds, ecType* key, idxType *sz, idxType id, idxType *inheap);
void randHeapInsert(idxType *heapIds, ecType* key, idxType *sz, idxType i, idxType *inheap);
void heapify(idxType *heapIds, ecType* key, idxType sz, idxType i,  idxType *inheap);
idxType heapExtractMax(idxType *heapIds, ecType* key, idxType *sz, idxType *inheap);
idxType randHeapExtractMax(idxType *heapIds, ecType* key, idxType *sz, idxType *inheap);
void heapBuild(idxType *heapIds, ecType* key, idxType sz, idxType *inheap);
void heapVerify(idxType *heapIds, ecType* key, idxType sz, idxType *inheap);

void randHeapBuild(idxType *heapIds, ecType* key, idxType sz, idxType *inheap);
void shuffleTab(idxType deb, idxType end, idxType* shuffle);

void maxHeapInsertKeyvals(idxType *heapIds, idxType *sz, idxType i, ecType *keyvals, vwType *vw);
void maxHeapifyKeyvals(idxType *heapIds, idxType sz, idxType i, ecType *keyvals, vwType *vw);
idxType maxHeapExtractKeyvals(idxType *heapIds, idxType *sz, ecType *keyvals, vwType *vw);
void maxBuildHeapKeyvals(idxType *heapIds, idxType sz, ecType *keyvals, vwType *vw);



#endif
