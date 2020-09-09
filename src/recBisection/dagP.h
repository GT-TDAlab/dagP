#ifndef DAGP_
#define DAGP_


#ifdef __cplusplus
extern "C" {
#endif

#include "dgraph.h"
#include "option.h"


int dagP_read_graph(char* file_name, dgraph *G, const MLGP_option *opt);
int dagP_init_filename(MLGP_option* opt, char* file_name);
int dagP_init_parameters(MLGP_option *opt, const int nbPart);
int dagP_opt_reallocUBLB(MLGP_option *opt, const int nbPart);
ecType dagP_partition_from_dgraph(dgraph *G, const MLGP_option *opt, idxType* parts);
int dagP_free_graph(dgraph* G);
int dagP_free_option(MLGP_option *opt);

#ifdef __cplusplus
}
#endif

#endif
