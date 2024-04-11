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
//#include "bib/PSO.h"
#include "bib/GM.h"

int main(){

    struct PARAMETERS BPR_PARAMETERS;
    struct MATRIZ_OD OD;
    igraph_vector_int_init(&OD.fontes, 0);
    igraph_vector_int_init(&OD.alvos, 0);
   
    igraph_vector_int_t edges;
    
    igraph_vector_int_init(&edges, 0);
    
    int** edge_list = init_parameters(&BPR_PARAMETERS,&edges);

    igraph_t Grafo;
    igraph_empty(&Grafo, BPR_PARAMETERS.N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);

    OD.MATRIZ = (int**) malloc(BPR_PARAMETERS.N*sizeof(int*));
    for (int i = 0; i < BPR_PARAMETERS.N; i++) OD.MATRIZ[i] = (int*) calloc(BPR_PARAMETERS.N,sizeof(int));
    load_MATOD(&OD);

    igraph_vector_t solucao;
    double** matrix_solution = (double**)malloc(BPR_PARAMETERS.L*sizeof(double*));
    for (int i = 0; i < BPR_PARAMETERS.L; i++)matrix_solution[i] = (double*)calloc(OD.N_ALVOS*(OD.N_FONTES-1),sizeof(double));
    optimize(&BPR_PARAMETERS,edge_list,&OD,&Grafo,&solucao,matrix_solution);
    //igraph_vector_print(&solucao);
    int* soma =(int*) calloc(BPR_PARAMETERS.N,sizeof(int));
    double* s =(double*) calloc(BPR_PARAMETERS.N,sizeof(double));
    for (int i = 0; i < BPR_PARAMETERS.L; i++){
        soma[edge_list[i][0]] += VECTOR(solucao)[i];
        soma[edge_list[i][1]] += VECTOR(solucao)[i];
        s[edge_list[i][0]] -= VECTOR(solucao)[i];
        s[edge_list[i][1]] += VECTOR(solucao)[i];
    }
    print_vetor(soma,BPR_PARAMETERS.N,sizeof(int));
    print_vetor(s,BPR_PARAMETERS.N,sizeof(double)); 
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
    igraph_destroy(&Grafo);
}