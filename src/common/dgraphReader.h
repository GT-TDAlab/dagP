/*
 *
 *   MM reader prt is from
 *   From Matrix Market I/O library for ANSI C
 *
 *   See http://math.nist.gov/MatrixMarket for details.
 *
 *
 */
#ifndef DGRAPHREADER_H
#define DGRAPHREADER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "dgraph.h"
#include "utils.h"


int readDGraphDot(dgraph *G, char* file_name);


int loadFromCSC(dgraph *G, idxType nVrtx, idxType nEdge, idxType *col_ptrs, idxType *row_inds, ecType *vals);
int readDGraph(dgraph *G, char* file_name, int prefer_binary);
int readDGraphDG(dgraph *G, char* file_name);
int readDGraphMtx(dgraph *G,char* file_name);
int readDGraphBinary(dgraph *G, char* file_name);
int readDGraphFab(dgraph *G, char* file_name);
int writePartsToTxt(dgraph *G, char* file_name, idxType* part);
void writePartGraphDot(dgraph* G, char* file_name, idxType* part);
int writeDGraphDot(dgraph *G, char* file_name, idxType* part);
int writeDGraphDotWithProc(dgraph *G, char* file_name, idxType* part, idxType* proc);
int writeDGraphBinary(dgraph *G, char* file_name);
idxType readVector(char* file_name, idxType* vector);

#ifdef __cplusplus
}
#endif


#endif
