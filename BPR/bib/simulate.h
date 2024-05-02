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

void obstructed(struct PARAMETERS* BPR_PARAMETERS,int** edge_list,igraph_vector_t* solucao,struct MATRIZ_OD* OD,igraph_vector_int_t* edges){
    int i,j;
    struct PARAMETERS BPR_PARAMETERS_obstructed;

    BPR_PARAMETERS_obstructed.N = BPR_PARAMETERS->N;
    BPR_PARAMETERS_obstructed.L = BPR_PARAMETERS->L;

    igraph_vector_init(&BPR_PARAMETERS_obstructed.capacidade, BPR_PARAMETERS_obstructed.L);
    igraph_vector_init(&BPR_PARAMETERS_obstructed.cost_time, BPR_PARAMETERS_obstructed.L);

    for ( i = 0; i < BPR_PARAMETERS->L; i++){
        VECTOR(BPR_PARAMETERS_obstructed.capacidade)[i] = VECTOR(BPR_PARAMETERS->capacidade)[i];
        VECTOR(BPR_PARAMETERS_obstructed.cost_time)[i] = VECTOR(BPR_PARAMETERS->cost_time)[i];
    }

    double** resultados = malloc(BPR_PARAMETERS->L* sizeof(double*));
    for (i = 0; i < BPR_PARAMETERS->L; i++) resultados[i] = calloc(BPR_PARAMETERS->L, sizeof(double)); // Aloca espaço para cada string

    FILE *file;
    file = fopen("./output/resultados_1.dat","w");

    double** matrix_solution2 = (double**)malloc((BPR_PARAMETERS->L)*sizeof(double*));
    for (int i = 0; i < BPR_PARAMETERS->L; i++)matrix_solution2[i] = (double*)calloc(OD->N_ALVOS*(OD->N_FONTES-1),sizeof(double));

    for (i = 0; i < BPR_PARAMETERS->L; i++){

        VECTOR(BPR_PARAMETERS_obstructed.capacidade)[i] = 1;

        igraph_t Grafo_obstructed;
        igraph_empty(&Grafo_obstructed, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
        igraph_add_edges(&Grafo_obstructed, edges, NULL);
        igraph_vector_t solucao_obstructed;

        optimize(&BPR_PARAMETERS_obstructed,edge_list,OD,&Grafo_obstructed,&solucao_obstructed,matrix_solution2);

        for (j = 0; j < BPR_PARAMETERS->L; j++){
            if(i == j)resultados[i][j] = 0;
            else resultados[i][j] =  VECTOR(*solucao)[j] - VECTOR(solucao_obstructed)[j];
        }

        igraph_vector_destroy(&solucao_obstructed);
        igraph_destroy(&Grafo_obstructed);
        printf("%d/%d\n",i+1,BPR_PARAMETERS->L);
        //if(i == 1) break;
        VECTOR(BPR_PARAMETERS_obstructed.capacidade)[i] = VECTOR(BPR_PARAMETERS->capacidade)[i];
    }
    for (i = 0; i < BPR_PARAMETERS->L; i++) {
        for (j = 0; j < BPR_PARAMETERS->L; j++) {
            fprintf(file, "%.2f ", resultados[i][j]);
        }
        fprintf(file, "\n"); // Quebra de linha entre as linhas da matriz
    }
    fclose(file);
    igraph_vector_destroy(&BPR_PARAMETERS_obstructed.capacidade);
    igraph_vector_destroy(&BPR_PARAMETERS_obstructed.cost_time);
    free(resultados);

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
    //remove_nodes(&Grafo,&BPR_PARAMETERS);

    int** indexate = (int**) malloc(BPR_PARAMETERS.N*sizeof(int*));

    for (int i = 0; i < BPR_PARAMETERS.N; i++) indexate[i] = (int*) calloc(BPR_PARAMETERS.N,sizeof(int));

    init_OD(&OD,&BPR_PARAMETERS,indexate);
    igraph_vector_t solucao;
    double** matrix_solution = (double**)malloc(BPR_PARAMETERS.L*sizeof(double*));
    for (int i = 0; i < BPR_PARAMETERS.L; i++)matrix_solution[i] = (double*)calloc(OD.N_ALVOS*(OD.N_FONTES-1),sizeof(double));

    optimize(&BPR_PARAMETERS,edge_list,&OD,&Grafo,&solucao,matrix_solution);
    printf("Finalizou a solução!\n");
    obstructed(&BPR_PARAMETERS,edge_list,&solucao,&OD,&edges);
    //igraph_vector_print(&solucao);

    
    
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
    igraph_destroy(&Grafo);
}