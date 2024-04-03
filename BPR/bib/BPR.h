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

double frank_wolfe(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo,igraph_vector_t* tempo,igraph_vector_t* direcao,const double* objetivo_incial,igraph_vector_t* y,double* stp){
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
    *stp = passo;
    return objetivo;
    
}

void optimize(struct PARAMETERS* BPR_PARAMETERS,int** edge_list,struct MATRIZ_OD *OD,igraph_t* Grafo,igraph_vector_t *solucao,double** matrix_solution){

    igraph_vector_t tempo;
    igraph_vector_init(&tempo, BPR_PARAMETERS->L);

    int i,iteracoes = 0,j;
    
    igraph_vector_init(solucao,BPR_PARAMETERS->L);
    igraph_vector_t gradiente;
    igraph_vector_init(&gradiente,BPR_PARAMETERS->L);
    double objetivo = 0,objetivo2 = 0,dx,dy,df;

    atualiza_fluxo(Grafo,OD,edge_list,solucao,&BPR_PARAMETERS->cost_time,matrix_solution);

    igraph_vector_t novo_fluxo;
    igraph_vector_init(&novo_fluxo,BPR_PARAMETERS->L);

    double** matrix_solution2 = (double**)malloc(BPR_PARAMETERS->L*sizeof(double*));
    int n = OD->N_ALVOS*OD->N_FONTES;
    double stp = 0;
    for (i = 0; i < BPR_PARAMETERS->L; i++)matrix_solution2[i] = (double*)calloc(n,sizeof(double));
    while(true){
        for ( i = 0; i < BPR_PARAMETERS->L; i++) for (int j = 0; j < n; j++) matrix_solution2[i][j] = 0;
        BPR(&tempo,BPR_PARAMETERS,solucao,&objetivo);
        atualiza_fluxo(Grafo,OD,edge_list,&gradiente, &tempo,matrix_solution2);
        
        for ( i = 0; i < BPR_PARAMETERS->L; i++){
            VECTOR(gradiente)[i] -= VECTOR(*solucao)[i];
            for (j = 0; j < n; j++) matrix_solution2[i][j] -= matrix_solution[i][j];
        }

        objetivo2 = frank_wolfe(BPR_PARAMETERS,solucao,&tempo,&gradiente,&objetivo,&novo_fluxo,&stp);
        for ( i = 0; i < BPR_PARAMETERS->L; i++) for ( j = 0; j < n; j++) matrix_solution[i][j] += stp*matrix_solution2[i][j];
        dx = 0.0;
		for(i=0;i<BPR_PARAMETERS->L;i++){
			dy = fabs(VECTOR(*solucao)[i]-VECTOR(novo_fluxo)[i]) / (fabs(VECTOR(*solucao)[i]) + 1);
			if(dx < dy) dx = dy;
            VECTOR(*solucao)[i] = VECTOR(novo_fluxo)[i];
		}
        df = (objetivo - objetivo2) / objetivo;
        

        iteracoes++;
        if(iteracoes > MAXIMO_ITERACOES) break;
        if(dx < X_TOLERANCIA) break;
        if(df < 1E-10) break;
        //printf("%d\n",iteracoes);
    }
    for (i = 0; i < BPR_PARAMETERS->L; i++){
        if(matrix_solution[i][0] != 0) printf("%d %d %f\n",edge_list[i][0],edge_list[i][1],matrix_solution[i][0]);
        
        free(matrix_solution2[i]);
    }
    free(matrix_solution2);
    igraph_vector_destroy(&novo_fluxo);
    igraph_vector_destroy(&tempo);
    igraph_vector_destroy(&gradiente);
}