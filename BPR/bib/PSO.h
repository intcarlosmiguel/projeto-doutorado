#pragma once
#include "mtwister.h"

struct posicao{
    int** x;
    int** best;
};
struct particula{
    struct posicao position;
    double** velocidade;
};

void print_matrix(int** m,int N,int L){
    printf("=================================================================================================== \n");
    int i,j;
    for ( i = 0; i < N; i++){
        for ( j = 0; j < L; j++){
            if(j!= L-1) printf("%d ",m[i][j]);
            else printf("%d\n",m[i][j]);
        }
    }
    printf("=================================================================================================== \n");
}

void new_best(int** best,int** x,int N){
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            best[i][j] = x[i][j];
}

void new_best_particle(struct particula* particulas,int particula,int N){
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            particulas[particula].position.best[i][j] = particulas[particula].position.x[i][j];
    
}

void init_particula(struct particula* particulas,int n_particulas,int N,int C){
    int valor;
    for (int i = 0; i < n_particulas; i++){

        particulas[i].velocidade = (double**) malloc(N*sizeof(double*));
        particulas[i].position.x = (int**) malloc(N*sizeof(int*));
        particulas[i].position.best = (int**) malloc(N*sizeof(int*));

        for (int j = 0; j < N; j++){
            particulas[i].velocidade[j] = (double*) malloc(N*sizeof(double));
            particulas[i].position.x[j] = (int*) malloc(N*sizeof(int));
            particulas[i].position.best[j] = (int*) malloc(N*sizeof(int));

            for (int k = 0; k < N; k++){
                valor =  C*genrand64_real1();
                particulas[i].velocidade[j][k] = valor;
                valor = C*genrand64_real1();
                particulas[i].position.x[j][k] = valor;
                particulas[i].position.best[j][k] = (int)particulas[i].position.x[j][k];
                if(j == k){
                    particulas[i].position.x[j][k] = 0;
                    particulas[i].position.best[j][k] = 0;
                    particulas[i].velocidade[j][k] = 0;
                }
            }
            
        }
        
    }
    
}

void atualiza_particula(struct particula* particulas,int particula,int N,int** best){
    double new_velocidade;
    int new_posicao;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if(i != j){
                new_velocidade = 0.5*particulas[particula].velocidade[i][j] + 0.1*genrand64_real1()*(particulas[particula].position.best[i][j] - particulas[particula].position.x[i][j]) + 0.1*genrand64_real1()*(best[i][j] - particulas[particula].position.x[i][j]);
                new_posicao = particulas[particula].position.x[i][j] + new_velocidade;
                particulas[particula].velocidade[i][j] = new_velocidade;

                particulas[particula].position.x[i][j] =(new_posicao > 0)? new_posicao:0;
            }
        }
    }

    
}

int objetivo(igraph_vector_t *solucao,igraph_vector_t *possivel_solucao,int L){
    int soma = 0;
    for (int i = 0; i < L; i++) soma += pow(VECTOR(*solucao)[i] - VECTOR(*possivel_solucao)[i],2);
    return soma/L;
}

void PSO(int n_particulas,int C,struct PARAMETERS* BPR_PARAMETERS,int** edge_list,igraph_t *Grafo,igraph_vector_int_t* fontes, igraph_vector_int_t* alvos,igraph_vector_t *solucao){

    int i;
    init_genrand64(42);
    struct particula* particulas = (struct particula*) malloc(n_particulas*sizeof(struct particula));

    double* objetivo_particula = (double*) malloc(BPR_PARAMETERS->N*sizeof(double));
    int** melhor_posicao = (int**) malloc(BPR_PARAMETERS->N*sizeof(int*));
    for ( i = 0; i < BPR_PARAMETERS->N; i++){
        melhor_posicao[i] = (int*) calloc(BPR_PARAMETERS->N,sizeof(int));
        objetivo_particula[i] = 1e20;
    }

    double melhor_objetivo = 1e20,obj;
    
    init_particula(particulas,n_particulas,BPR_PARAMETERS->N,C);
    double media;
    int iteracoes = 0;
    while (true){
        media = 0;
        for ( i = 0; i < n_particulas; i++){
            
            igraph_vector_t possivel_solucao;

            optimize(BPR_PARAMETERS,edge_list,particulas[i].position.x,Grafo,fontes,alvos,&possivel_solucao);
            //igraph_vector_print(&possivel_solucao);
            obj = objetivo(solucao,&possivel_solucao,BPR_PARAMETERS->L);
            //printf("%f\n",obj);
            if(obj < objetivo_particula[i]){
                new_best_particle(particulas,i,BPR_PARAMETERS->N);
                objetivo_particula[i] = obj;
            }
            if(obj < melhor_objetivo){
                new_best(melhor_posicao,particulas[i].position.x,BPR_PARAMETERS->N);
                melhor_objetivo = obj;
            }
            
            igraph_vector_destroy(&possivel_solucao);
            media += obj;
        }
        for ( i = 0; i < n_particulas; i++){
            atualiza_particula(particulas,i,BPR_PARAMETERS->N,melhor_posicao);
            //
        }
        media /= n_particulas;
        printf("%f %f %d\n",media,melhor_objetivo,iteracoes);
        if(abs(media - melhor_objetivo) <= 1e-4) break;
        //break;
        iteracoes++;
        if(iteracoes == 100) break;
    }
    print_matrix(melhor_posicao,BPR_PARAMETERS->N,BPR_PARAMETERS->N);
    free(objetivo_particula);
}