#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "igraph.h"
#include "calc.h"
#include "define.h"
#include "BPR.h"
#include "network.h"
//#include "PSO.h"
#include "GM.h"

void remove_nodes(igraph_t *Grafo,struct PARAMETERS* BPR_PARAMETERS){
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_vector_t peso1;
    igraph_vector_init(&peso1,BPR_PARAMETERS->L);
    for (int i = 0; i < BPR_PARAMETERS->L; i++)VECTOR(peso1)[i] = VECTOR(BPR_PARAMETERS->cost_time)[i];
    printf("%ld %ld\n",igraph_ecount(Grafo),igraph_vector_size(&BPR_PARAMETERS->cost_time));
    printf("Erro\n");
    SETEANV(Grafo, "cost_time", &peso1);
    SETEANV(Grafo, "capacity", &BPR_PARAMETERS->capacidade);
    igraph_vector_int_t membership, csize;
    igraph_integer_t no_of_clusters;

    igraph_vector_int_init(&membership, 0);
    igraph_vector_int_init(&csize, 0);

    igraph_connected_components(Grafo, &membership, &csize, &no_of_clusters, IGRAPH_STRONG);
    int index_of_largest = 0;
    int size_of_largest = VECTOR(csize)[0];

    for (int i = 1; i < igraph_vector_int_size(&csize); i++) {
        if (VECTOR(csize)[i] > size_of_largest) {
            size_of_largest = VECTOR(csize)[i];
            index_of_largest = i;
        }
    }
    
    //printf("O maior cluster é o cluster %d com %d vértices.\n", index_of_largest, size_of_largest);
    igraph_vs_t vs;
    igraph_vector_int_t delete_vertices;
    igraph_vector_int_init(&delete_vertices, 0);
    for (int i = 0; i < BPR_PARAMETERS->N; i++){
        if(VECTOR(membership)[i] != index_of_largest) igraph_vector_int_push_back(&delete_vertices, i);
    }
    igraph_vs_vector(&vs, &delete_vertices);

    igraph_delete_vertices(Grafo, vs);
    igraph_vector_int_destroy(&membership);
    igraph_vector_int_destroy(&csize);
    igraph_vector_int_destroy(&delete_vertices);
    igraph_vs_destroy(&vs);
    BPR_PARAMETERS->N = size_of_largest;
    BPR_PARAMETERS->L = igraph_ecount(Grafo);
}

void simulate(){
    struct PARAMETERS BPR_PARAMETERS;
    struct MATRIZ_OD OD;
    //igraph_vector_int_init(&OD.fontes, 0);
    //igraph_vector_int_init(&OD.alvos, 0);
   
    igraph_vector_int_t edges;
    
    igraph_vector_int_init(&edges, 0);
    int** edge_list = init_parameters(&BPR_PARAMETERS,&edges,false);
    igraph_t Grafo;
    igraph_empty(&Grafo, BPR_PARAMETERS.N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);

    remove_nodes(&Grafo,&BPR_PARAMETERS);

    printf("%d %d\n",BPR_PARAMETERS.N,BPR_PARAMETERS.L);
    exit(0);
    int** indexate = (int**) malloc(BPR_PARAMETERS.N*sizeof(int*));

    for (int i = 0; i < BPR_PARAMETERS.N; i++) indexate[i] = (int*) calloc(BPR_PARAMETERS.N,sizeof(int));

    init_OD(&OD,&BPR_PARAMETERS,indexate);
    igraph_vector_t solucao;
    double** matrix_solution = (double**)malloc(BPR_PARAMETERS.L*sizeof(double*));
    for (int i = 0; i < BPR_PARAMETERS.L; i++)matrix_solution[i] = (double*)calloc(OD.N_ALVOS*(OD.N_FONTES-1),sizeof(double));

    optimize(&BPR_PARAMETERS,edge_list,&OD,&Grafo,&solucao,matrix_solution);

    igraph_vector_print(&solucao);
    /*load_MATOD(&OD,true);

    igraph_vector_t solucao;
    double** matrix_solution = (double**)malloc(BPR_PARAMETERS.L*sizeof(double*));
    for (int i = 0; i < BPR_PARAMETERS.L; i++)matrix_solution[i] = (double*)calloc(OD.N_ALVOS*(OD.N_FONTES-1),sizeof(double));
    optimize(&BPR_PARAMETERS,edge_list,&OD,&Grafo,&solucao,matrix_solution);
    //igraph_vector_print(&solucao);
    //caminhos(matrix_solution,&BPR_PARAMETERS,edge_list,0,6,8);
    //GM(&OD,&BPR_PARAMETERS,edge_list,&Grafo,&solucao);
    //igraph_vector_print(&solucao);

    //PSO(50,500,&BPR_PARAMETERS,edge_list,&Grafo,&fontes, &alvos,&solucao);

    igraph_vector_destroy(&BPR_PARAMETERS.cost_time);
    igraph_vector_destroy(&BPR_PARAMETERS.capacidade);
    igraph_vector_int_destroy(&OD.fontes);
    igraph_vector_int_destroy(&OD.alvos);
    igraph_vector_int_destroy(&edges);
    for (int i = 0; i < BPR_PARAMETERS.N; i++) free(OD.MATRIZ[i]);
    free(OD.MATRIZ);
    igraph_destroy(&Grafo);*/
}