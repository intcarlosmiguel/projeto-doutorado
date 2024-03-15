#pragma once
#include <igraph.h>

int find_id(int fonte,int alvo,int** edge_list,int L){
    int i;
    for (i = 0; i < L; i++)if((edge_list[i][0] == fonte)&&(edge_list[i][1] == alvo)) break;
    return i;
}

void Dijkstra(igraph_t* Grafo,int* fonte,igraph_vector_int_t* alvos,igraph_vector_t* pesos,igraph_vector_int_t* parents){

    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t inbound;
    igraph_vs_t alvos_vs;
    igraph_vs_vector(&alvos_vs,alvos);
    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    
    igraph_vector_int_init(&inbound, 0);
    igraph_get_shortest_paths_dijkstra(Grafo, &vecs,&evecs,*fonte,alvos_vs,pesos, IGRAPH_ALL,parents,&inbound);
    //igraph_vector_int_print(parents);
    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vs_destroy(&alvos_vs);
    igraph_vector_int_destroy(&inbound);
}

void atualiza_fluxo(igraph_t *Grafo,double** MATRIZ_OD,int** edge_list,igraph_vector_t* fluxo,igraph_vector_int_t *fontes,igraph_vector_int_t *alvos, igraph_vector_t *pesos,int L){
    int i,j,fonte,alvo,antecessor,index;
    double volume;
    igraph_vector_fill(fluxo,0);
    for ( i = 0; i < igraph_vector_int_size(fontes); i++){
        fonte = VECTOR(*fontes)[i];
        //printf("=============== FONTE: %d ===============\n",fonte);
        igraph_vector_int_t parents;
        igraph_vector_int_init(&parents, 0);
        
        Dijkstra(Grafo,&fonte,alvos,pesos,&parents);

        for ( j = 0; j < igraph_vector_int_size(alvos); j++){
            alvo = VECTOR(*alvos)[j];
            //printf("Alvo: %d\n",alvo);
            volume = MATRIZ_OD[fonte][alvo];
            while (alvo != fonte){
                antecessor = VECTOR(parents)[alvo];
                index = find_id(antecessor,alvo,edge_list,L);
                //printf("%d : (%d,%d) - %f\n",index,antecessor,alvo,volume);
                VECTOR(*fluxo)[index] += volume;
                alvo = antecessor;
            }
        }
        igraph_vector_int_destroy(&parents);
    }
}