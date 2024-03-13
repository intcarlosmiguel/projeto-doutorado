#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "igraph.h"
#include "bib/calc.h"
#include "bib/define.h"

int** init_parameters(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_int_t* edges,int*N){
    char nomeDoArquivo[800];
    sprintf(nomeDoArquivo,"./file/dial_edges_algbformat.txt"); 
    int size,i;
    double** data = lerArquivo(nomeDoArquivo, 7,&size);
    igraph_vector_init(&BPR_PARAMETERS->capacidade, 0);
    igraph_vector_init(&BPR_PARAMETERS->cost_time, 0);
    //igraph_vector_init(&BPR_PARAMETERS->velocidade, 0);
    int** edge_list = (int**)malloc(size*sizeof(int*));
    *N = size;
    for (i = 0; i < size; i++){
        edge_list[i] = (int*)calloc(2,sizeof(int));
        edge_list[i][0] = data[i][1]-1;
        edge_list[i][1] = data[i][2]-1;
        
        igraph_vector_int_push_back(edges,edge_list[i][0]);
        igraph_vector_int_push_back(edges,edge_list[i][1]);
        igraph_vector_push_back(&BPR_PARAMETERS->capacidade,data[i][4]);
        igraph_vector_push_back(&BPR_PARAMETERS->cost_time,data[i][6]);
        BPR_PARAMETERS->L += 1;
        free(data[i]);
    }
    free(data);
    return edge_list;
}

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

void BPR(igraph_vector_t* gradiente,struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo,double *total_time){
    double s;
    *total_time = 0;
    for (int i = 0; i < BPR_PARAMETERS->L; i++){
        s = pow(VECTOR(*fluxo)[i]/VECTOR(BPR_PARAMETERS->capacidade)[i],BETA);
        *total_time += (VECTOR(BPR_PARAMETERS->cost_time)[i] * VECTOR(*fluxo)[i]) * (1.0 + 0.03 * s);
        VECTOR(*gradiente)[i] = (VECTOR(BPR_PARAMETERS->cost_time)[i])*(1 + ALPHA*(1+BETA)*s);
    }
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

double optimize(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo,igraph_vector_t* tempo,igraph_vector_t* direcao,const double* objetivo_incial,igraph_vector_t* y){
    double direcao_gradiente = 0, passo = 1.0;
    int i,iteracoes = 0;
    for ( i = 0; i < BPR_PARAMETERS->L; i++) direcao_gradiente += VECTOR(*direcao)[i]*VECTOR(*tempo)[i];
    if(direcao_gradiente > 0){
        printf("Gradiente não ótimo!\n");
        exit(0);
    }
    double objetivo = 0;
    double dgtest = ftol * direcao_gradiente;
    
    
    while (true){
        for ( i = 0; i < BPR_PARAMETERS->L; i++) VECTOR(*y)[i] = VECTOR(*fluxo)[i] + VECTOR(*direcao)[i]*passo;
        BPR(tempo,BPR_PARAMETERS,y,&objetivo);
        iteracoes++;
        if (objetivo > *objetivo_incial + passo * dgtest) passo *= decresse;
        else break;

        if(passo < min_step) break;
        if(passo > max_step) break;
		if(iteracoes > MAXIMO_ITERACOES) break;

    }
    return objetivo;
    
}

int main(){

    struct PARAMETERS BPR_PARAMETERS;
    
    int i,L,N = 25,iteracoes = 0;
    double objetivo,dx,dy,df,objetivo2;
   
    igraph_vector_int_t edges;
    
    igraph_vector_int_init(&edges, 0);
    
    int** edge_list = init_parameters(&BPR_PARAMETERS,&edges,&L);

    igraph_vector_int_t fontes;
    igraph_vector_int_init(&fontes, 0);
    igraph_vector_int_t alvos;
    igraph_vector_int_init(&alvos, 0);
    igraph_vector_t tempo;
    igraph_vector_init(&tempo, L);

    double** MATRIZ_OD = (double**) malloc(N*sizeof(double*));
    for (i = 0; i < N; i++) MATRIZ_OD[i] = (double*) calloc(N,sizeof(double));
    load_MATOD(MATRIZ_OD,&fontes,&alvos);

    igraph_t Grafo;
    igraph_empty(&Grafo, N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    
    igraph_vector_t volumes;
    igraph_vector_init(&volumes,L);
    igraph_vector_t direcao;
    igraph_vector_init(&direcao,L);
    atualiza_fluxo(&Grafo,MATRIZ_OD,edge_list,&volumes,&fontes,&alvos, &BPR_PARAMETERS.cost_time, L);
    while(true){

        igraph_vector_t novo_fluxo;

        igraph_vector_init(&novo_fluxo,BPR_PARAMETERS.L);

        BPR(&tempo,&BPR_PARAMETERS,&volumes,&objetivo);

        atualiza_fluxo(&Grafo,MATRIZ_OD,edge_list,&direcao,&fontes,&alvos, &tempo, L);

        for ( i = 0; i < L; i++) VECTOR(direcao)[i] -= VECTOR(volumes)[i];

        objetivo2 = optimize(&BPR_PARAMETERS,&volumes,&tempo,&direcao,&objetivo,&novo_fluxo);

        dx = 0.0;
		for(i=0;i<BPR_PARAMETERS.L;i++){
			dy = fabs(VECTOR(volumes)[i]-VECTOR(novo_fluxo)[i]) / (fabs(VECTOR(volumes)[i]) + 1);
			if(dx < dy) dx = dy;
            VECTOR(volumes)[i] = VECTOR(novo_fluxo)[i];
		}
        df = (objetivo - objetivo2) / objetivo;

        igraph_vector_destroy(&novo_fluxo);

        iteracoes++;
        if(iteracoes > MAXIMO_ITERACOES) break;
        if(dx < X_TOLERANCIA) break;
        if(df < 1E-10) break;

    }
    double saldo = 0;
    int index,vizinho;
    for ( i = 0; i < N; i++){
        saldo = 0;
        igraph_vector_int_t vizinhos;
        igraph_vector_int_init(&vizinhos, 0);
        igraph_neighbors(&Grafo, &vizinhos, i,IGRAPH_ALL);
        for (int j = 0; j < igraph_vector_int_size(&vizinhos); j++){
            vizinho = VECTOR(vizinhos)[j];
            index = find_id(vizinho,i,edge_list,L);
            saldo += VECTOR(volumes)[index];
            index = find_id(i,vizinho,edge_list,L);
            saldo -= VECTOR(volumes)[index];
            //printf("%d : (%d,%d) - %f\n",index,antecessor,alvo,volume);
        }
        printf("%d - %f\n",i+1,saldo);
    }
    
    igraph_vector_print(&volumes);
    igraph_vector_destroy(&BPR_PARAMETERS.capacidade);
    igraph_vector_destroy(&tempo);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&direcao);
    igraph_vector_int_destroy(&fontes);
    igraph_vector_int_destroy(&alvos);
    igraph_vector_destroy(&BPR_PARAMETERS.cost_time);
    igraph_destroy(&Grafo);
    for (i = 0; i < N; i++) free(MATRIZ_OD[i]);
    free(MATRIZ_OD);
}