#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "igraph.h"
#include "bib/calc.h"
#include "bib/define.h"

#define ALPHA 0.15
#define BETA 4

struct PARAMETERS{
    igraph_vector_t capacidade;
    igraph_vector_t cost_time;
    igraph_vector_int_t dual;
    int L;
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
    igraph_get_shortest_paths_dijkstra(Grafo, &vecs,&evecs,*fonte,igraph_vss_all(),pesos, IGRAPH_ALL,parents,&inbound);
    //igraph_vector_int_print(parents);
    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vs_destroy(&alvos_vs);
    igraph_vector_int_destroy(&inbound);
}

void BPR(igraph_vector_t* tempo,struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo,double *f){
    double s;
    *f = 0;
    for (int i = 0; i < BPR_PARAMETERS->L; i++){
        s = pow(VECTOR(*fluxo)[i]/VECTOR(BPR_PARAMETERS->capacidade)[i],BETA);
        *f += (VECTOR(BPR_PARAMETERS->cost_time)[i] * VECTOR(*fluxo)[i]) * (1.0 + 0.03 * s);
        VECTOR(*tempo)[i] = (VECTOR(BPR_PARAMETERS->cost_time)[i])*(1 + ALPHA*s);
    }
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
        BPR_PARAMETERS->L += 1;
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

void update_x(igraph_t *Grafo,double** MATRIZ_OD,int** edge_list,igraph_vector_t* fluxo,igraph_vector_int_t *fontes,igraph_vector_int_t *alvos, igraph_vector_t *pesos,int L){
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
        printf("Gradiente não ótimo!");
        exit(0);
    }
    double objetivo = 0;
    double dgtest = ftol * direcao_gradiente;
    
    
    while (true){
        for ( i = 0; i < BPR_PARAMETERS->L; i++) VECTOR(*y)[i] = VECTOR(*fluxo)[i] - VECTOR(*direcao)[i]*passo;
        BPR(tempo,BPR_PARAMETERS,y,&objetivo);
        iteracoes++;
        if (objetivo > *objetivo_incial + passo * dgtest) passo *= decresse;
        else return objetivo;

        if(passo < min_step) break;
        if(passo > max_step) break;
		if(iteracoes > MAXIMO_ITERACOES) break;

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
    igraph_empty(&Grafo, N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    
    int fonte,alvo,antecessor,index;
    igraph_vector_t volumes;
    igraph_vector_init(&volumes,L);
    igraph_vector_t direcao;
    igraph_vector_init(&direcao,L);
    //igraph_vector_int_print(&fontes);
    update_x(&Grafo,MATRIZ_OD,edge_list,&volumes,&fontes,&alvos, &BPR_PARAMETERS.cost_time, L);
    igraph_vector_print(&volumes);
    int iteracoes = 0;
    double objetivo,dx,dy,df,objetivo2;
    while(true){

        igraph_vector_t novo_fluxo;
        igraph_vector_init(&novo_fluxo,BPR_PARAMETERS.L);

        BPR(&tempo,&BPR_PARAMETERS,&volumes,&objetivo);

        update_x(&Grafo,MATRIZ_OD,edge_list,&direcao,&fontes,&alvos, &tempo, L);

        igraph_vector_print(&direcao);

        for ( i = 0; i < L; i++) VECTOR(direcao)[i] -= VECTOR(volumes)[i];

        objetivo2 = optimize(&BPR_PARAMETERS,&volumes,&tempo,&direcao,&objetivo,&novo_fluxo);

        dx = 0.0;
		for(i=0;i<BPR_PARAMETERS.L;i++){
			dy = fabs(VECTOR(volumes)[i]-VECTOR(novo_fluxo)[i]) / (fabs(VECTOR(volumes)[i]) + 1);
			if(dx < dy) dx = dy;
		}
        df = (objetivo - objetivo2) / objetivo;

        igraph_vector_init_copy(&volumes,&novo_fluxo);

        igraph_vector_destroy(&novo_fluxo);

        iteracoes++;
        if(iteracoes > MAXIMO_ITERACOES) break;
        if(dx < X_TOLERANCIA) break;
        if(df < 1E-4) break;

        break;
    }
    igraph_vector_destroy(&BPR_PARAMETERS.capacidade);
    igraph_vector_destroy(&tempo);
    igraph_vector_destroy(&BPR_PARAMETERS.cost_time);
}