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
    igraph_vector_int_t dual;
};

void load_MATOD(double** MATRIZ_OD,igraph_vector_int_t* fontes,igraph_vector_int_t* alvos){
    char nomeDoArquivo[800];
    int size,i,site1,site2;
    sprintf(nomeDoArquivo,"./file/dial_matod.txt");
    double** data = lerArquivo(nomeDoArquivo, 3,&size);
    
    for (i = 0; i < size; i++){
        site1 = data[i][0] - 1;
        site2 = data[i][1] - 1;
        MATRIZ_OD[site1][site2] = data[i][2];
        if(!igraph_vector_int_contains(fontes,site1)) igraph_vector_int_push_back(fontes,site1 );
        if(!igraph_vector_int_contains(alvos,site2)) igraph_vector_int_push_back(alvos,site2 );
        
    }

}

int find_id(int fonte,int alvo,int** edge_list,int L){
    for (int i = 0; i < L; i++)if((edge_list[i][0] == fonte)&&(edge_list[i][1] == alvo)) return i;
}

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

void BPR(igraph_vector_t* tempo,struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo){
    int L = igraph_vector_size(&BPR_PARAMETERS->capacidade);
    for (int i = 0; i < L; i++) VECTOR(*tempo)[i] = (VECTOR(BPR_PARAMETERS->cost_time)[i])*(1 + ALPHA*pow(VECTOR(*fluxo)[i]/VECTOR(BPR_PARAMETERS->capacidade)[i],BETA));
}

int** init_parameters(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_int_t* edges,int*N){
    char nomeDoArquivo[800];
    sprintf(nomeDoArquivo,"./file/dial_edges_algbformat.txt"); 
    int size,i;
    double** data = lerArquivo(nomeDoArquivo, 7,&size);
    igraph_vector_init(&BPR_PARAMETERS->capacidade, 0);
    igraph_vector_init(&BPR_PARAMETERS->cost_time, 0);
    //igraph_vector_init(&BPR_PARAMETERS->velocidade, 0);
    igraph_vector_int_init(&BPR_PARAMETERS->dual, 0);
    igraph_vector_int_t pointer;
    igraph_vector_int_init(&pointer, size);
    int** edge_list = (int**)malloc(size*sizeof(int*));
    *N = size;
    for (i = 0; i < size; i++){
        edge_list[i] = (int*)calloc(2,sizeof(int));
        edge_list[i][0] = data[i][1]-1;
        edge_list[i][1] = data[i][2]-1;
        
        igraph_vector_int_push_back(edges,edge_list[i][0]);
        igraph_vector_int_push_back(edges,edge_list[i][1]);
        igraph_vector_push_back(&BPR_PARAMETERS->capacidade,data[i][3]);
        igraph_vector_push_back(&BPR_PARAMETERS->cost_time,data[i][6]);
        igraph_vector_int_push_back(&BPR_PARAMETERS->dual,data[i][3]);
        VECTOR(pointer)[i] = -1;
        free(data[i]);
    }
    int index;
    for ( i = 0; i < size; i++){
        //printf("%ld %ld\n",VECTOR(BPR_PARAMETERS->dual)[i],VECTOR(pointer)[i]);
        if((VECTOR(BPR_PARAMETERS->dual)[i] == 1) && (VECTOR(pointer)[i] == -1)){
            index = find_id(edge_list[i][1],edge_list[i][0],edge_list,size);
            VECTOR(pointer)[i] = index;
            VECTOR(pointer)[index] = i;
        }
    }
    igraph_vector_int_init_copy(&BPR_PARAMETERS->dual,&pointer);
    igraph_vector_int_destroy(&pointer);
    free(data);
    return edge_list;
}

void update_x(igraph_t *Grafo,double** MATRIZ_OD,int** edge_list,igraph_vector_int_t* DUAL,igraph_vector_t* fluxo,igraph_vector_int_t *fontes,igraph_vector_int_t *alvos, igraph_vector_t *pesos,int L){
    int i,j,fonte,alvo,antecessor,index;
    double volume;
    for ( i = 0; i < igraph_vector_int_size(fontes); i++){
        fonte = VECTOR(*fontes)[i];
        //printf("=============== FONTE: %d ===============\n",fonte);
        igraph_vector_int_t parents;
        igraph_vector_int_init(&parents, 0);
        Dijkstra(Grafo,&fonte,alvos,pesos,&parents);
        for ( j = 0; j < igraph_vector_int_size(alvos); j++){
            alvo = VECTOR(*alvos)[j];
            volume = MATRIZ_OD[fonte][alvo];
            if(MATRIZ_OD[fonte][alvo] == 0) break;
            while (alvo != fonte){
                antecessor = VECTOR(parents)[alvo];
                index = find_id(antecessor,alvo,edge_list,L);
                VECTOR(*fluxo)[index] += volume;
                if(VECTOR(*DUAL)[index] > 0 ){
                    index = VECTOR(*DUAL)[index];
                    VECTOR(*fluxo)[index] += volume;
                }
                alvo = antecessor;
            }
        }
        igraph_vector_int_destroy(&parents);
    }
}

void main(){

    struct PARAMETERS BPR_PARAMETERS;
    
    int i,j,L,size; // Substitua 5 pelo número real de colunas que você quer processar
   
    igraph_vector_int_t edges;
    
    igraph_vector_int_init(&edges, 0);
    
    int** edge_list = init_parameters(&BPR_PARAMETERS,&edges,&L);

    igraph_vector_int_t fontes;
    igraph_vector_int_init(&fontes, 0);
    igraph_vector_int_t alvos;
    igraph_vector_int_init(&alvos, 0);
    igraph_vector_t tempo;
    igraph_vector_init(&tempo, L);
    
    int N = 25;

    double** MATRIZ_OD = (double**) malloc(N*sizeof(double*));
    for (i = 0; i < N; i++) MATRIZ_OD[i] = (double*) calloc(N,sizeof(double));
    load_MATOD(MATRIZ_OD,&fontes,&alvos);

    igraph_t Grafo;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_empty(&Grafo, N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    
    //for ( i = 0; i < L; i++)printf("Ligação: %d: (%d,%d)\n",i,edge_list[i][0],edge_list[i][1]);
    
    //exit(0);
    int fonte,alvo,antecessor,index;
    igraph_vector_t volumes;
    igraph_vector_init(&volumes,L);
    igraph_vector_t direcao;
    igraph_vector_init(&direcao,L);

    update_x(&Grafo,MATRIZ_OD,edge_list,&BPR_PARAMETERS.dual,&volumes,&fontes,&alvos, &BPR_PARAMETERS.cost_time, L);
    //igraph_vector_int_print(&fontes);
    //igraph_vector_int_print(&alvos);
    //for ( i = 0; i < L; i++) printf("%d - (%d,%d) %f\n",i+1,edge_list[i][0],edge_list[i][1],volumes[i]);
    
    //printf("# de Ligações: %ld\n",igraph_ecount(&Grafo));
    double objetivo;
    
    while(true){
        BPR(&tempo,&BPR_PARAMETERS,&volumes);
        update_x(&Grafo,MATRIZ_OD,edge_list,&BPR_PARAMETERS.dual,&tempo,&fontes,&alvos, &tempo, L);
        objetivo = igraph_vector_sum(&tempo);
        igraph_vector_scale(&volumes,-1);
        
        igraph_vector_add(&direcao,&volumes);
        igraph_vector_print(&direcao);
        //igraph_vector_scale(fluxo,-1);
        //igraph_vector_add(direcao,fluxo);
        break;
    }
    //BPR(tempo,edges,capacidade,cost_time,fluxo);
    igraph_vector_destroy(&BPR_PARAMETERS.capacidade);
    igraph_vector_destroy(&tempo);
    igraph_vector_destroy(&BPR_PARAMETERS.cost_time);
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