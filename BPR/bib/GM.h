#pragma once
#include "bib/network.h"
#include "bib/BPR.h"
#include "bib/calc.h"
#include "mtwister.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void indexate(int** index,struct MATRIZ_OD *OD_i,struct PARAMETERS* BPR_PARAMETERS){
    int i,j,c = 0;
    for ( i = 0; i < BPR_PARAMETERS->N; i++){
        for (int j = 0; j < BPR_PARAMETERS->N; j++){
            if((i!=j) && (OD_i->MATRIZ[i][j] != 0)){
                index[i][j] = c;
                c++;
            }
            else index[i][j] = -1;
        }
    }
    
}

void init_OD(struct MATRIZ_OD *OD_i,struct PARAMETERS* BPR_PARAMETERS){
    igraph_vector_int_init(&OD_i->fontes, 0);
    igraph_vector_int_init(&OD_i->alvos, 0);
    igraph_vector_int_range(&OD_i->fontes,0,BPR_PARAMETERS->N);
    igraph_vector_int_range(&OD_i->alvos,0,BPR_PARAMETERS->N);
    OD_i->MATRIZ = (int**) malloc(BPR_PARAMETERS->N*sizeof(int*));
    for (int i = 0; i < BPR_PARAMETERS->N; i++){
        OD_i->MATRIZ[i] = (int*) calloc(BPR_PARAMETERS->N,sizeof(int));
        for (int j = 0; j < BPR_PARAMETERS->N; j++) if(i !=j) OD_i->MATRIZ[i][j] = 500*genrand64_real1();
    }
}

void GM(struct MATRIZ_OD *OD,struct PARAMETERS* BPR_PARAMETERS,int** edge_list,igraph_t *Grafo,igraph_vector_t *solucao){

    int i,j;
    igraph_vector_t possivel_solucao;
    struct MATRIZ_OD OD_i;
    int* index = (int*) malloc(BPR_PARAMETERS->N*sizeof(int));
    init_OD(&OD_i,BPR_PARAMETERS);

    double** dZ_dOD = (double**)malloc(BPR_PARAMETERS->N*sizeof(double*));
    for (i = 0; i < BPR_PARAMETERS->N; i++)dZ_dOD[i] = (double*)calloc(BPR_PARAMETERS->N,sizeof(double));

    while (true){
        double** matrix_solution = (double**)malloc(BPR_PARAMETERS->L*sizeof(double*));
        for (i = 0; i < BPR_PARAMETERS->L; i++) matrix_solution[i] = (double*)calloc(OD_i.N_ALVOS*OD_i.N_FONTES,sizeof(double));
        optimize(BPR_PARAMETERS,edge_list,&OD_i,Grafo,&possivel_solucao,matrix_solution);
        igraph_vector_sub(&possivel_solucao,solucao);
        for (i = 0; i < BPR_PARAMETERS->N; i++)
            for (j = 0; j < BPR_PARAMETERS->N; j++)
                if((i != j) && (OD_i.MATRIZ[i][j] != 0))
                    //dZ_dOD[i][j] = (1/OD_i.MATRIZ[i][j])*;
        break;
    }
    
}