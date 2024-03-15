#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "igraph.h"
#include "bib/calc.h"
#include "bib/define.h"
#include "bib/BPR.h"
#include "bib/network.h"

int main(){

    struct PARAMETERS BPR_PARAMETERS;
    igraph_vector_int_t fontes;
    igraph_vector_int_init(&fontes, 0);
    igraph_vector_int_t alvos;
    igraph_vector_int_init(&alvos, 0);
   
    igraph_vector_int_t edges;
    
    igraph_vector_int_init(&edges, 0);
    
    int** edge_list = init_parameters(&BPR_PARAMETERS,&edges);


    double** MATRIZ_OD = (double**) malloc(BPR_PARAMETERS.N*sizeof(double*));
    for (int i = 0; i < BPR_PARAMETERS.N; i++) MATRIZ_OD[i] = (double*) calloc(BPR_PARAMETERS.N,sizeof(double));
    load_MATOD(MATRIZ_OD,&fontes,&alvos);
    
    optimize(&BPR_PARAMETERS,edge_list,MATRIZ_OD,&edges,&fontes,&alvos);

    
}