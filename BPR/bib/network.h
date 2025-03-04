#pragma once
#include <igraph.h>
#include <calc.h>
#include <omp.h>

void print_vector_igraph(igraph_vector_int_t* vetor){
    int N =igraph_vector_int_size(vetor);
    for (int i = 0; i < N; i++){
        if(i != N - 1) printf("%ld,",VECTOR(*vetor)[i]);
        else printf("%ld\n",VECTOR(*vetor)[i]);
    }
    
}

int find_id(int fonte,int alvo,int** edge_list){
    int i = 0;
    while(true){
        if((edge_list[i][0] == fonte)&&(edge_list[i][1] == alvo)) break;
        i++;
    }
    return i;
}

void Dijkstra(igraph_t* Grafo,int fonte,igraph_vector_int_t* alvos,igraph_vector_t* pesos,igraph_vector_int_t* parents){

    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t inbound;
    igraph_vs_t alvos_vs;
    igraph_vs_vector(&alvos_vs,alvos);
    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    
    igraph_vector_int_init(&inbound, 0);
    igraph_get_shortest_paths_dijkstra(Grafo, &vecs,&evecs,fonte,alvos_vs,pesos, IGRAPH_OUT,parents,&inbound);
    //igraph_get_shortest_paths_bellman_ford(Grafo, &vecs,&evecs,*fonte,alvos_vs,pesos, IGRAPH_IN,parents,&inbound);
    //igraph_vector_int_print(parents);
    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vs_destroy(&alvos_vs);
    igraph_vector_int_destroy(&inbound);
}

void atualiza_fluxo(igraph_t *Grafo,struct MATRIZ_OD* OD,int** edge_list,igraph_vector_t* fluxo, igraph_vector_t *pesos){
    int i,j,fonte,c = 0;
    igraph_vector_fill(fluxo,0);
    for ( i = 0; i < OD->N_FONTES; i++){
        fonte = VECTOR(OD->fontes)[i];
        igraph_vector_int_t parents;
        igraph_vector_int_init(&parents, 0);
        Dijkstra(Grafo,fonte,&OD->alvos,pesos,&parents);
        //printf("Fonte: %d\n",fonte);
        //print_vector_igraph(&parents);
        #pragma omp parallel for schedule(dynamic) reduction(+:c)
        for (j = 0; j < OD->N_ALVOS; j++) {
            int alvo = VECTOR(OD->alvos)[j];
            double volume = OD->MATRIZ[fonte][alvo];
            if(VECTOR(parents)[alvo] < 0) continue;
            while (alvo != fonte) {
                int antecessor = VECTOR(parents)[alvo];
                int index = find_id(antecessor, alvo, edge_list);
                #pragma omp atomic
                VECTOR(*fluxo)[index] += volume;
                alvo = antecessor;
            }
            if(VECTOR(OD->alvos)[j] != fonte) c++;
        }
        igraph_vector_int_destroy(&parents);
    }
}