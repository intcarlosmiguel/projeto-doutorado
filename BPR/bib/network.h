#pragma once

#include <calc.h>
#include "igraph.h"

void print_vector_igraph(igraph_vector_int_t* vetor){
    int N =igraph_vector_int_size(vetor);
    for (int i = 0; i < N; i++){
        if(i != N - 1) printf("%ld,",VECTOR(*vetor)[i]);
        else printf("%ld\n",VECTOR(*vetor)[i]);
    }
    
}

int find_id_log(int fonte, int alvo, int** edge_list,int n) {
    int left = 0;
    int right = n - 1;  // assuming Grafo is accessible
    
    while (left <= right) {
        int mid = left + (right - left) / 2;
        
        if (edge_list[mid][0] < fonte) {
            left = mid + 1;
        } else if (edge_list[mid][0] > fonte) {
            right = mid - 1;
        } else {
            // Found matching source, now check target
            if (edge_list[mid][1] < alvo) {
                left = mid + 1;
            } else if (edge_list[mid][1] > alvo) {
                right = mid - 1;
            } else {
                return mid;  // Found exact match
            }
        }
    }
    return -1;  // Not found (error case)
}

int find_id(int fonte,int alvo,int** edge_list){
    int i = 0;
    while(true){
        if((edge_list[i][0] == fonte)&&(edge_list[i][1] == alvo)) break;
        i++;
    }
    return i;
}

void get_bush(igraph_t* Grafo,int fonte,igraph_vector_t* pesos){

    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t inbound,parents;
    igraph_t Bush;
    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    
    igraph_vector_int_init(&inbound, 0);
    igraph_vector_int_init(&parents, 0);

    igraph_get_shortest_paths_dijkstra(
        Grafo, &vecs,&evecs,fonte,igraph_vss_all(),pesos, IGRAPH_OUT,&parents,&inbound
    );

    printf("%ld %ld\n",igraph_vector_int_list_size(&evecs),igraph_ecount(Grafo));
    igraph_vector_int_destroy(&inbound);
    igraph_vector_int_destroy(&parents);

    igraph_es_t edges;
    igraph_vector_int_t edges_vec;
    igraph_vector_int_init(&edges_vec, 0);
    for (long int i = 0; i < igraph_vector_int_list_size(&evecs); i++) {
        igraph_vector_int_t* vec = igraph_vector_int_list_get_ptr(&evecs, i);
        for (long int j = 0; j < igraph_vector_int_size(vec); j++) {
            igraph_vector_int_push_back(&edges_vec, VECTOR(*vec)[j]);
        }
    }


    igraph_es_vector(&edges, &edges_vec);

    igraph_subgraph_from_edges(
        Grafo, &Bush, edges, false
    );
    printf("Grafo Bush:\n");
    printf("Number of edges in Bush: %ld\n", igraph_ecount(&Bush));
    igraph_destroy(&Bush);
    
    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vector_int_destroy(&edges_vec);
    igraph_es_destroy(&edges);
    exit(0);
}

void Dijkstra(igraph_t* Grafo,int fonte,igraph_vector_int_t* alvos,igraph_vector_t* pesos,igraph_vector_int_t* parents){

    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t inbound;
    igraph_vs_t alvos_vs;
    igraph_vs_vector(&alvos_vs,alvos);
    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    
    igraph_vector_int_init(&inbound, 0);
    igraph_get_shortest_paths_dijkstra(Grafo, &vecs,&evecs,fonte,igraph_vss_all(),pesos, IGRAPH_OUT,parents,&inbound);
    
    // Print each vector in the lis
    //igraph_get_shortest_paths_bellman_ford(Grafo, &vecs,&evecs,*fonte,alvos_vs,pesos, IGRAPH_IN,parents,&inbound);
    //igraph_vector_int_print(parents);
    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vs_destroy(&alvos_vs);
    igraph_vector_int_destroy(&inbound);
}

void parallel_atualiza_fluxo(igraph_t *Grafo,struct OD_MATRIX* OD,int** edge_list,igraph_vector_t* fluxo, igraph_vector_t *pesos){
    double** fluxos_paralelo = (double**) malloc(THREADS * sizeof(double*));

    int i,j,fonte,c = 0,thread_id;
    for (i = 0; i < THREADS; i++){
        fluxos_paralelo[i] = (double*) calloc( igraph_vector_size(fluxo), sizeof(double));
    }
    igraph_vector_fill(fluxo,0);
    #pragma omp parallel for schedule(dynamic)
    for ( i = 0; i < OD->size; i++){
        thread_id = omp_get_thread_num();
        fonte = OD->Elementos[i].fonte;
        igraph_vector_int_t parents;
        igraph_vector_int_init(&parents, 0);
        Dijkstra(Grafo,fonte,&OD->Elementos[i].alvos,pesos,&parents);

        //printf("Fonte: %d\n",fonte);
        //print_vector_igraph(&parents);
        for (j = 0; j < igraph_vector_int_size(&OD->Elementos[i].alvos); j++) {
            int alvo = VECTOR(OD->Elementos[i].alvos)[j];
            double volume = VECTOR(OD->Elementos[i].volumes)[j];
            if(VECTOR(parents)[alvo] < 0) continue;
            while (alvo != fonte) {
                int antecessor = VECTOR(parents)[alvo];
                int index = find_id(antecessor, alvo, edge_list);
                fluxos_paralelo[thread_id][index] += volume;
                //VECTOR(*fluxo)[index] += volume;
                alvo = antecessor;
            }
            if(VECTOR(OD->Elementos[i].alvos)[j] != fonte) c++;
        }
        igraph_vector_int_destroy(&parents);
    }
    for (i = 0; i < THREADS; i++){
        for (j = 0; j < igraph_vector_size(fluxo); j++){
            VECTOR(*fluxo)[j] += fluxos_paralelo[i][j];
        }
        free(fluxos_paralelo[i]);
    }
    

}
void atualiza_fluxo(igraph_t *Grafo,struct OD_MATRIX* OD,int** edge_list,igraph_vector_t* fluxo, igraph_vector_t *pesos){
    int i,j,fonte,c = 0;
    igraph_vector_fill(fluxo,0);
    for ( i = 0; i < OD->size; i++){
        fonte = OD->Elementos[i].fonte;
        igraph_vector_int_t parents;
        igraph_vector_int_init(&parents, 0);
        Dijkstra(Grafo,fonte,&OD->Elementos[i].alvos,pesos,&parents);

        //printf("Fonte: %d\n",fonte);
        //print_vector_igraph(&parents);
        for (j = 0; j < igraph_vector_int_size(&OD->Elementos[i].alvos); j++) {
            int alvo = VECTOR(OD->Elementos[i].alvos)[j];
            double volume = VECTOR(OD->Elementos[i].volumes)[j];
            if(VECTOR(parents)[alvo] < 0) continue;
            while (alvo != fonte) {
                int antecessor = VECTOR(parents)[alvo];
                int index = find_id(antecessor, alvo, edge_list);
                VECTOR(*fluxo)[index] += volume;
                alvo = antecessor;
            }
            if(VECTOR(OD->Elementos[i].alvos)[j] != fonte) c++;
        }
        igraph_vector_int_destroy(&parents);
    }
}
