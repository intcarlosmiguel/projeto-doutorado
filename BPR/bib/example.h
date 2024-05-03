#pragma once

#include <stdbool.h>
#include "igraph.h"
#include <network.h>
#include <calc.h>
#include <BPR.h>


int** init_parameters_example(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_int_t* edges){

    char nomeDoArquivo[800];
    sprintf(nomeDoArquivo,"./file/example/dial_edges_algbformat.txt");

    int size,i;

    double** data = lerArquivo(nomeDoArquivo, 4,&size);

    igraph_vector_init(&BPR_PARAMETERS->capacidade, 0);
    igraph_vector_init(&BPR_PARAMETERS->cost_time, 0);

    int** edge_list = (int**)malloc(size*sizeof(int*));
    BPR_PARAMETERS->L = 0;
    BPR_PARAMETERS->N = 0;

    for (i = 0; i < size; i++){

        edge_list[i] = (int*)calloc(2,sizeof(int));
        edge_list[i][0] = data[i][1]-1;
        edge_list[i][1] = data[i][2]-1;
        
        igraph_vector_int_push_back(edges,edge_list[i][0]);
        igraph_vector_int_push_back(edges,edge_list[i][1]);
        igraph_vector_push_back(&BPR_PARAMETERS->capacidade,data[i][4]);
        igraph_vector_push_back(&BPR_PARAMETERS->cost_time,data[i][6]);
        BPR_PARAMETERS->L += 1;
        if(BPR_PARAMETERS->N < edge_list[i][0]+1) BPR_PARAMETERS->N = edge_list[i][0]+1;
            

        free(data[i]);
    }

    free(data);

    return edge_list;
}