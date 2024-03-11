#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "igraph.h"
#include "bib/calc.h"

#define ALPHA 0.15
#define BETA 4

struct PARAMETERS{
    igraph_vector_t capacidade;
    igraph_vector_t cost_time;
    igraph_vector_t velocidade;
};


void Dijkstra(igraph_t* Grafo,int* fonte,igraph_vector_int_t* alvos,igraph_vector_t* pesos,igraph_vector_int_t* parents){

    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t inbound;
    igraph_vs_t alvos_vs;
    igraph_vs_vector(&alvos_vs,alvos);
    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    
    igraph_vector_int_init(&inbound, 0);

    igraph_get_shortest_paths_dijkstra(Grafo, &vecs,&evecs,*fonte,alvos_vs,pesos, IGRAPH_IN,parents,&inbound);
    //igraph_vector_int_print(parents);
    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vs_destroy(&alvos_vs);
    
}

void BPR(igraph_vector_t* tempo,igraph_vector_int_t* edges,igraph_vector_t* capacidade,igraph_vector_t* cost_time,double* fluxo){
    int L = igraph_vector_size(capacidade);
    for (int i = 0; i < L; i++) VECTOR(*tempo)[i] = (VECTOR(*cost_time)[i])*(1 + ALPHA*pow(fluxo[i]/VECTOR(*capacidade)[i],BETA));
}

void main(){
    char nomeDoArquivo[800];
    sprintf(nomeDoArquivo,"./file/dial_edges_algbformat.txt"); // Substitua "meuArquivo.txt" pelo nome real do seu arquivo
    int numeroDeColunas = 7,i,j,size; // Substitua 5 pelo número real de colunas que você quer processar
    double** data = lerArquivo(nomeDoArquivo, numeroDeColunas,&size);

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    igraph_vector_t capacidade;
    igraph_vector_init(&capacidade, 0);
    igraph_vector_t velocidade;
    igraph_vector_init(&velocidade, 0);
    igraph_vector_t tempo;
    igraph_vector_init(&tempo, 0);
    igraph_vector_t cost_time;
    igraph_vector_init(&cost_time, 0);

    igraph_vector_int_t fontes;
    igraph_vector_int_init(&fontes, 0);
    igraph_vector_int_t alvos;
    igraph_vector_int_init(&alvos, 0);
    
    int N = 25;
    int total_ligacoes = N*(N-1);
    bool* ligacoes = (bool*) calloc(total_ligacoes,sizeof(bool));
    double* volume = (double*) calloc(total_ligacoes,sizeof(double));
    int site1,site2;
    for (i = 0; i < size; i++){
        site1 = data[i][1]-1;
        site2 = data[i][2]-1;
        if((data[i][5] == 1) && (!ligacoes[site1*(N-1)+site2-1]) && (!ligacoes[site2*(N-1)+site1-1])){
            igraph_vector_int_push_back(&edges,site1);
            igraph_vector_int_push_back(&edges,site2);
            ligacoes[site1*(N-1)+site2-1] = true;
            ligacoes[site2*(N-1)+site1-1] = true;
            igraph_vector_push_back(&capacidade,data[i][3]);
            igraph_vector_push_back(&velocidade,data[i][4]);
            igraph_vector_push_back(&cost_time,data[i][6]);
        }
    }
    free(ligacoes);
    sprintf(nomeDoArquivo,"./file/dial_matod.txt");
    numeroDeColunas = 3;
    data = lerArquivo(nomeDoArquivo, numeroDeColunas,&size);
    bool* sitios = (bool*) calloc(N,sizeof(bool));
    double** MATRIZ_OD = (double**) malloc(N*sizeof(double*));
    for (i = 0; i < N; i++) MATRIZ_OD[i] = (double*) calloc(N,sizeof(double));
    for (i = 0; i < size; i++){
        site1 = data[i][0] - 1;
        site2 = data[i][1] - 1;
        MATRIZ_OD[site1][site2] = data[i][2];
        if(!igraph_vector_int_contains(&fontes,site1)){
            igraph_vector_int_push_back(&fontes,site1 );
            sitios[i] = true;
        }
        if(!igraph_vector_int_contains(&alvos,site2)){
            igraph_vector_int_push_back(&alvos,site2 );
            sitios[i] = true;
        }
        
    }

    igraph_t Grafo;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_empty(&Grafo, N, IGRAPH_UNDIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    int fonte,alvo;
    for ( i = 0; i < igraph_vector_int_size(&fontes); i++){
        fonte = VECTOR(fontes)[i];
        printf("Fonte: %d\n",fonte);
        igraph_vector_int_t parents;
        igraph_vector_int_init(&parents, 0);
        Dijkstra(&Grafo,&fonte,&alvos,&cost_time,&parents);
        for ( j = 0; j < igraph_vector_int_size(&alvos); j++){
            alvo = VECTOR(alvos)[j];
            printf("Alvo: %d\n",alvo);
            while (alvo != fonte){
                volume[fonte*(N-1)+alvo-1] += MATRIZ_OD[fonte][alvo];
                alvo = VECTOR(parents)[alvo];
            }
        }
        

        igraph_vector_int_destroy(&parents);
    }
    
    for ( i = 0; i < total_ligacoes; i++){
        printf("%d %d %f\n",i/25,i%25,volume[i]);
    }
    
    
    printf("# de Ligações: %ld\n",igraph_ecount(&Grafo));
    while(true){
        //BPR(tempo,edges,capacidade,cost_time,fluxo);
        //double fy = igraph_vector_sum(&fluxo);
        //igraph_vector_scale(fluxo,-1);
        //igraph_vector_add(direcao,fluxo);
        break;
    }
    //BPR(tempo,edges,capacidade,cost_time,fluxo);
    igraph_vector_destroy(&capacidade);
    igraph_vector_destroy(&tempo);
    igraph_vector_destroy(&velocidade);
    igraph_vector_destroy(&cost_time);
    free(sitios);
}

/*void main() {

    igraph_t g;
    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t parents, inbound;
    igraph_integer_t i;
    igraph_real_t weights[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_vector_t weights_vec;
    igraph_vs_t vs;

    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    igraph_vector_int_init(&parents, 0);
    igraph_vector_int_init(&inbound, 0);

    igraph_vs_vector_small(&vs, 0, 1, 3, 5, 2, 1, -1);
    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,   1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,         3, 2, 3, 6,
                 4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,         8, 9,
                 5, 2,
                 2, 1,
                 -1);

    igraph_vector_view(&weights_vec, weights, sizeof(weights) / sizeof(weights[0]));
    
    printf("Vertices:\n");
    for (i = 0; i < igraph_vector_int_list_size(&vecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&vecs, i));
    }

    printf("\nEdges:\n");
    for (i = 0; i < igraph_vector_int_list_size(&evecs); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&evecs, i));
    }

    printf("\nParents:\n");
    igraph_vector_int_print(&parents);

    printf("\nInbound:\n");
    igraph_vector_int_print(&inbound);

    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_destroy(&inbound);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

}*/