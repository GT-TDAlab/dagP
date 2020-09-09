#ifndef DEBUG_H_
#define DEBUG_H_

#include "dgraph.h"
#include "utils.h"

void createDGraph_example(dgraph* G);
void print_dgraph(dgraph G);
void print_graph(dgraph* G, idxType* part, char* path, int arg1, int arg2);
void dgraph_to_dot(dgraph* G, idxType* part, char* file_name);
void dgraph_to_dot_with_matching(dgraph* G, idxType* leader, char* file_name);

#endif
