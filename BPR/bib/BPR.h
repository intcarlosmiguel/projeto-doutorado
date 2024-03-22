#pragma once
#include <igraph.h>
#include <math.h>
#include <network.h>

struct PARAMETERS{
    igraph_vector_t capacidade;
    igraph_vector_t cost_time;
    int L;
    int N;
};
int** init_parameters(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_int_t* edges){
    char nomeDoArquivo[800];
    sprintf(nomeDoArquivo,"./file/dial_edges_algbformat.txt"); 
    int size,i;
    double** data = lerArquivo(nomeDoArquivo, 7,&size);
    igraph_vector_init(&BPR_PARAMETERS->capacidade, 0);
    igraph_vector_init(&BPR_PARAMETERS->cost_time, 0);
    //igraph_vector_init(&BPR_PARAMETERS->velocidade, 0);
    int** edge_list = (int**)malloc(size*sizeof(int*));
    BPR_PARAMETERS->L = 0;
    BPR_PARAMETERS->N = 0;
    for (i = 0; i < size; i++){
        edge_list[i] = (int*)calloc(2,sizeof(int));
        edge_list[i][0] = data[i][1]-1;
        edge_list[i][1] = data[i][2]-1;
        
        igraph_vector_int_push_back(edges,edge_list[i][0]);
        igraph_vector_int_push_back(edges,edge_list[i][1]);
        igraph_vector_push_back(&BPR_PARAMETERS->capacidade,data[i][4]);
        igraph_vector_push_back(&BPR_PARAMETERS->cost_time,data[i][6]);
        BPR_PARAMETERS->L += 1;
        if(BPR_PARAMETERS->N < edge_list[i][0]+1) BPR_PARAMETERS->N = edge_list[i][0]+1;
        free(data[i]);
    }
    free(data);
    return edge_list;
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

double frank_wolfe(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo,igraph_vector_t* tempo,igraph_vector_t* direcao,const double* objetivo_incial,igraph_vector_t* y){
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
        //printf("%d %f %f\n",iteracoes,objetivo,*objetivo_incial + passo * dgtest);
        if (objetivo > *objetivo_incial + passo * dgtest) passo *= decresse;
        else break;

        if(passo < min_step) break;
        if(passo > max_step) break;
		if(iteracoes > MAXIMO_ITERACOES) break;
    }
    return objetivo;
    
}

void optimize(struct PARAMETERS* BPR_PARAMETERS,int** edge_list,double** MATRIZ_OD,igraph_vector_int_t* edges,igraph_vector_int_t* fontes, igraph_vector_int_t* alvos){

    igraph_vector_t tempo;
    igraph_vector_init(&tempo, BPR_PARAMETERS->L);

    int i,iteracoes = 0;
    igraph_t Grafo;
    igraph_empty(&Grafo, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, edges, NULL);
    
    igraph_vector_t volumes;
    igraph_vector_init(&volumes,BPR_PARAMETERS->L);
    igraph_vector_t gradiente;
    igraph_vector_init(&gradiente,BPR_PARAMETERS->L);
    double objetivo = 0,objetivo2 = 0,dx,dy,df;

    atualiza_fluxo(&Grafo,MATRIZ_OD,edge_list,&volumes,fontes,alvos,&BPR_PARAMETERS->cost_time);

    igraph_vector_t novo_fluxo;
    igraph_vector_init(&novo_fluxo,BPR_PARAMETERS->L);
    while(true){

        BPR(&tempo,BPR_PARAMETERS,&volumes,&objetivo);
        atualiza_fluxo(&Grafo,MATRIZ_OD,edge_list,&gradiente,fontes,alvos, &tempo);
        
        for ( i = 0; i < BPR_PARAMETERS->L; i++) VECTOR(gradiente)[i] -= VECTOR(volumes)[i];

        objetivo2 = frank_wolfe(BPR_PARAMETERS,&volumes,&tempo,&gradiente,&objetivo,&novo_fluxo);

        dx = 0.0;
		for(i=0;i<BPR_PARAMETERS->L;i++){
			dy = fabs(VECTOR(volumes)[i]-VECTOR(novo_fluxo)[i]) / (fabs(VECTOR(volumes)[i]) + 1);
			if(dx < dy) dx = dy;
            VECTOR(volumes)[i] = VECTOR(novo_fluxo)[i];
		}
        df = (objetivo - objetivo2) / objetivo;
        

        iteracoes++;
        if(iteracoes > MAXIMO_ITERACOES) break;
        if(dx < X_TOLERANCIA) break;
        if(df < 1E-4) break;

    }
    igraph_vector_destroy(&novo_fluxo);
    /* double* soma = (double*) calloc(BPR_PARAMETERS->N,sizeof(double));
    for (i = 0; i < BPR_PARAMETERS->L; i++){
        soma[edge_list[i][0]] += VECTOR(volumes)[i];
        soma[edge_list[i][1]] -= VECTOR(volumes)[i];
    }
    print_vetor(soma,BPR_PARAMETERS->N,sizeof(double));

    free(soma); */

    igraph_vector_print(&volumes);
    igraph_vector_destroy(&BPR_PARAMETERS->capacidade);
    igraph_vector_destroy(&tempo);
    igraph_vector_int_destroy(edges);
    igraph_vector_destroy(&gradiente);
    igraph_vector_int_destroy(fontes);
    igraph_vector_int_destroy(alvos);
    igraph_vector_destroy(&BPR_PARAMETERS->cost_time);
    igraph_destroy(&Grafo);
    for (i = 0; i < BPR_PARAMETERS->N; i++) free(MATRIZ_OD[i]);
    free(MATRIZ_OD);
}