#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "igraph.h"
#include "bib/calc.h"

#define ALPHA 0.15
#define BETA 4

void BPR(igraph_vector_t* tempo,igraph_vector_int_t* edges,igraph_vector_t* capacidade,igraph_vector_t* cost_time,double* fluxo){
    int L = igraph_vector_size(capacidade);
    for (int i = 0; i < L; i++) VECTOR(*tempo)[i] = (VECTOR(*cost_time)[i])*(1 + ALPHA*pow(fluxo[i]/VECTOR(*capacidade)[i],BETA));
}

void main(){
    igraph_t Grafo;
    const char *nomeDoArquivo = "./file/dial_edges_algbformat.txt"; // Substitua "meuArquivo.txt" pelo nome real do seu arquivo
    int numeroDeColunas = 7; // Substitua 5 pelo número real de colunas que você quer processar

    // Chama a função lerArquivo com o nome do arquivo e o número de colunas especificados
    //printf("%d\n",contarLinhasNoArquivo(nomeDoArquivo));
    int size;
    double** data = lerArquivo(nomeDoArquivo, numeroDeColunas,&size);

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 2*size);
    igraph_vector_t capacidade;
    igraph_vector_init(&capacidade, size);
    igraph_vector_t velocidade;
    igraph_vector_init(&velocidade, size);
    igraph_vector_t tempo;
    igraph_vector_init(&tempo, size);
    igraph_vector_t cost_time;
    igraph_vector_init(&cost_time, size);
    igraph_vector_t direcao;
    igraph_vector_init(&direcao, size);
    
    int N = 0;
    for (int i = 0; i < size; i++){
        VECTOR(edges)[2*i] = data[i][1] - 1;
        VECTOR(edges)[2*i+1] = data[i][2] - 1;
        
        VECTOR(capacidade)[i] = data[i][3];
        VECTOR(velocidade)[i] = data[i][4];
        VECTOR(cost_time)[i] = data[i][4];
        //printf("%f\n",VECTOR(cost_time)[i]);
        if(data[i][1] > N) N = data[i][1];
    }
    igraph_matrix_t distance;
    igraph_matrix_init(&distance, 0, 0);
    while(true){
        //BPR(tempo,edges,capacidade,cost_time,fluxo);
        //double fy = igraph_vector_sum(&fluxo);
        //igraph_t Grafo;
        //igraph_set_attribute_table(&igraph_cattribute_table);
        //igraph_empty(&Grafo, N, IGRAPH_DIRECTED);
        //igraph_add_edges(&Grafo, &edges, NULL);
        //igraph_cattribute_VAN_setv(&Grafo,"tempo",&tempo);
        //igraph_distances_dijkstra(&Grafo, &distance, igraph_vss_all(), igraph_vss_all(),&tempo, IGRAPH_OUT);
        //igraph_vector_scale(fluxo,-1);
        //igraph_vector_add(direcao,fluxo);
        break;
    }
    //BPR(tempo,edges,capacidade,cost_time,fluxo);
    igraph_vector_destroy(&capacidade);
    igraph_vector_destroy(&tempo);
    igraph_vector_destroy(&velocidade);
    igraph_vector_destroy(&cost_time);
}