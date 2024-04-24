#pragma once
#include "network.h"
#include "BPR.h"
#include "calc.h"
#include "mtwister.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void caminhos(double** matrix_solution,struct PARAMETERS* BPR_PARAMETERS,int** edge_list,int par,int fonte,int alvo){
    int i,c = 0;

    igraph_t subgrafo;

    igraph_vector_int_t conexoes,map,indexate,res,h,teste;

    igraph_vector_int_init(&conexoes,0);
    igraph_vector_int_init(&map,0);
    igraph_vector_int_init(&indexate,BPR_PARAMETERS->N);
    igraph_vector_int_fill(&indexate,-1);
    igraph_vector_int_init(&res,0);
    igraph_vector_int_init(&h,0);
    igraph_vector_int_init(&teste,0);


    for (i = 0; i < BPR_PARAMETERS->L; i++){
        //if(fonte == 0) if(alvo == 1) if(matrix_solution[i][par]!= 0)printf("%d %d %.20f\n",edge_list[i][0],edge_list[i][1],matrix_solution[i][par]);
        if(matrix_solution[i][par] != 0) {
            if(!igraph_vector_int_contains(&map,edge_list[i][0])){
                igraph_vector_int_push_back(&map,edge_list[i][0]);
                VECTOR(indexate)[edge_list[i][0]] = c;
                c++;
            }

            igraph_vector_int_push_back(&conexoes,VECTOR(indexate)[edge_list[i][0]]);

            if(!igraph_vector_int_contains(&map,edge_list[i][1])){
                igraph_vector_int_push_back(&map,edge_list[i][1]);
                VECTOR(indexate)[edge_list[i][1]] = c;
                c++;
            }
            igraph_vector_int_push_back(&conexoes,VECTOR(indexate)[edge_list[i][1]]);
        }
        
    }
    igraph_empty(&subgrafo, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
    igraph_add_edges(&subgrafo, &conexoes, NULL);
    igraph_get_all_simple_paths(&subgrafo,&res,VECTOR(indexate)[fonte],igraph_vss_1(VECTOR(indexate)[alvo]),-1,IGRAPH_OUT);
    //igraph_vector_int_print(&res);
    for (i = 0; i < igraph_vector_int_size(&res); i++){
        if(VECTOR(res)[i] != -1)printf("%ld ",VECTOR(map)[VECTOR(res)[i]]);
        else printf("\n");
    }
    igraph_vector_int_destroy(&teste);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_destroy(&conexoes);
    igraph_vector_int_destroy(&indexate);
    igraph_vector_int_destroy(&res);

    igraph_destroy(&subgrafo);
}

void calcula_derivada(double** matrix_solution,struct PARAMETERS* BPR_PARAMETERS,int** edge_list,int par,int fonte,int alvo,double** dZ_dOD,igraph_vector_t *possivel_solucao){
    
    int i,c = 0,index,current;
    double current_valor;

    igraph_t subgrafo;

    igraph_vector_int_t conexoes,map,indexate,res;

    igraph_vector_int_init(&conexoes,0);
    igraph_vector_int_init(&map,0);
    igraph_vector_int_init(&indexate,BPR_PARAMETERS->N);
    igraph_vector_int_fill(&indexate,-1);
    igraph_vector_int_init(&res,0);

    //printf("OD - %d = (%d %d)\n",par,fonte,alvo);

    for (i = 0; i < BPR_PARAMETERS->L; i++){
        //print_vetor(matrix_solution[i],12,sizeof(double));
        if(matrix_solution[i][par] != 0) {
            //printf("%d %d\n",edge_list[i][0],edge_list[i][1]);

            if(!igraph_vector_int_contains(&map,edge_list[i][0])){
                igraph_vector_int_push_back(&map,edge_list[i][0]);
                VECTOR(indexate)[edge_list[i][0]] = c;
                c++;
            }

            igraph_vector_int_push_back(&conexoes,VECTOR(indexate)[edge_list[i][0]]);

            if(!igraph_vector_int_contains(&map,edge_list[i][1])){
                igraph_vector_int_push_back(&map,edge_list[i][1]);
                VECTOR(indexate)[edge_list[i][1]] = c;
                c++;
            }
            igraph_vector_int_push_back(&conexoes,VECTOR(indexate)[edge_list[i][1]]);
        }
        
    }
    igraph_empty(&subgrafo, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
    igraph_add_edges(&subgrafo, &conexoes, NULL);
    //printf("%ld %ld\n",VECTOR(indexate)[fonte],VECTOR(indexate)[alvo]);
    igraph_get_all_simple_paths(&subgrafo,&res,VECTOR(indexate)[fonte],igraph_vss_1(VECTOR(indexate)[alvo]),-1,IGRAPH_OUT);
    /* printf("=================================================================\n");
    for (i = 0; i < igraph_vector_int_size(&res); i++){
        if(VECTOR(res)[i] != -1)printf("%ld ",VECTOR(map)[VECTOR(res)[i]]);
        else printf("\n");
    } */
    //printf("=================================================================\n");
    current_valor = index = 0;
    for (i = 0; i < igraph_vector_int_size(&res); i++){
        //printf("%d %ld %ld\n",i,VECTOR(res)[i],VECTOR(map)[VECTOR(res)[i]]);
        if(VECTOR(res)[i] == -1){
            dZ_dOD[fonte][alvo] += current_valor*VECTOR(*possivel_solucao)[index];
            i++;
            continue;
        }
        current = VECTOR(map)[VECTOR(res)[i]];
        if(current == fonte){
            index = find_id(current,VECTOR(map)[VECTOR(res)[i+1]],edge_list);
            current_valor = matrix_solution[index][par];
        }
        else{
            if(VECTOR(res)[i+1] == -1) continue;
            index = find_id(current,VECTOR(map)[VECTOR(res)[i+1]],edge_list);
            if(current_valor != matrix_solution[index][par]){
                current_valor = matrix_solution[index][par];
                while(VECTOR(res)[i+1] != -1) i++;
            }
        }
        
    }

    //igraph_vector_int_print(&h);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_destroy(&conexoes);
    igraph_vector_int_destroy(&indexate);
    igraph_vector_int_destroy(&res);

    igraph_destroy(&subgrafo);
}


void init_OD(struct MATRIZ_OD *OD_i,struct PARAMETERS* BPR_PARAMETERS,int** indexate){

    igraph_vector_int_init(&OD_i->fontes, 0);
    igraph_vector_int_init(&OD_i->alvos, 0);
    igraph_vector_int_range(&OD_i->fontes,0,BPR_PARAMETERS->N);
    igraph_vector_int_range(&OD_i->alvos,0,BPR_PARAMETERS->N);
    OD_i->N_ALVOS = BPR_PARAMETERS->N;
    OD_i->N_FONTES = BPR_PARAMETERS->N;
    OD_i->MATRIZ = (int**) malloc(BPR_PARAMETERS->N*sizeof(int*));
    OD_i->LIST = (int**) malloc((BPR_PARAMETERS->N)*(BPR_PARAMETERS->N-1)*sizeof(int*));
    int c = 0;
    
    for (int i = 0; i < BPR_PARAMETERS->N; i++)OD_i->MATRIZ[i] = (int*) calloc(BPR_PARAMETERS->N,sizeof(int));
    for (int i = 0; i < BPR_PARAMETERS->N; i++){
        
        for (int j = 0; j < BPR_PARAMETERS->N; j++){
            if(i !=j){
                OD_i->LIST[c] = (int*) malloc(2*sizeof(int));
                OD_i->MATRIZ[i][j] = 500*genrand64_real1();
                OD_i->LIST[c][0] = i;
                OD_i->LIST[c][1] = j;
                indexate[i][j] = c;
                c++;
            }
        }
    }
}

void GM(struct MATRIZ_OD *OD,struct PARAMETERS* BPR_PARAMETERS,int** edge_list,igraph_t *Grafo,igraph_vector_t *solucao){

    int i,j;
    igraph_vector_print(solucao);
    struct MATRIZ_OD OD_i;
    int** indexate = (int**) malloc(BPR_PARAMETERS->N*sizeof(int*));
    /* OD_i.MATRIZ = (int**) malloc(BPR_PARAMETERS->N*sizeof(int*));
    for (int i = 0; i < BPR_PARAMETERS->N; i++) OD_i.MATRIZ[i] = (int*) calloc(BPR_PARAMETERS->N,sizeof(int));
    igraph_vector_int_init(&OD_i.fontes, 0);
    igraph_vector_int_init(&OD_i.alvos, 0);
    load_MATOD(&OD_i); */
    double** dZ_dOD = (double**)malloc(BPR_PARAMETERS->N*sizeof(double*));
    for (i = 0; i < BPR_PARAMETERS->N; i++){
        dZ_dOD[i] = (double*)calloc(BPR_PARAMETERS->N,sizeof(double));
        indexate[i] = (int*) calloc(BPR_PARAMETERS->N,sizeof(int));
    }
    init_OD(&OD_i,BPR_PARAMETERS,indexate);
    
    int n = OD_i.N_ALVOS*(OD_i.N_FONTES-1);
    //int fonte,alvo;
    //int c = 0;
    /* for ( i = 0; i < OD_i.N_FONTES; i++){
        fonte = VECTOR(OD_i.fontes)[i];
        for ( j = 0; j < OD_i.N_ALVOS; j++){
            alvo = VECTOR(OD_i.alvos)[j];
            if(fonte!= alvo) {
                indexate[fonte][alvo] = c;
                c++;
            }
        }
    } */
    double total,lambda,Z;
    while (true){ 
        Z = 0;
        lambda = 0;
        total = 0;
        double** matrix_solution = (double**)malloc(BPR_PARAMETERS->L*sizeof(double*));
        double *derivada = (double*)calloc(BPR_PARAMETERS->L,sizeof(double*));
        igraph_vector_t possivel_solucao;

        for (i = 0; i < BPR_PARAMETERS->L; i++) matrix_solution[i] = (double*)calloc(n,sizeof(double));

        optimize(BPR_PARAMETERS,edge_list,&OD_i,Grafo,&possivel_solucao,matrix_solution);

        //igraph_vector_print(&possivel_solucao);
        for ( i = 0; i < BPR_PARAMETERS->L; i++){
            VECTOR(possivel_solucao)[i] = (VECTOR(possivel_solucao)[i] -VECTOR(*solucao)[i]); 
            Z += VECTOR(possivel_solucao)[i]*VECTOR(possivel_solucao)[i];
        }
        
        for ( i = 0; i < n; i++) if(OD_i.MATRIZ[OD_i.LIST[i][0]][OD_i.LIST[i][1]] != 0)calcula_derivada(matrix_solution,BPR_PARAMETERS,edge_list,indexate[OD_i.LIST[i][0]][OD_i.LIST[i][1]],OD_i.LIST[i][0],OD_i.LIST[i][1],dZ_dOD,&possivel_solucao);
        for ( i = 0; i < BPR_PARAMETERS->N; i++){
            for (j = 0; j < BPR_PARAMETERS->N; j++) if(( i != j) && (OD_i.MATRIZ[i][j] != 0)) dZ_dOD[i][j] /=OD_i.MATRIZ[i][j]; 
            //printf("%f\n",dZ_dOD[0][1]);
            //print_vetor(dZ_dOD[i],BPR_PARAMETERS->N,sizeof(double));
        }

        for ( i = 0; i < BPR_PARAMETERS->L; i++)for ( j = 0; j < n; j++) derivada[i] += dZ_dOD[OD_i.LIST[j][0]][OD_i.LIST[j][1]]*matrix_solution[i][j];
        //print_vetor(derivada,BPR_PARAMETERS->L,sizeof(double));

        for ( i = 0; i < BPR_PARAMETERS->L; i++) total += derivada[i]*derivada[i];
        
        for ( i = 0; i < BPR_PARAMETERS->L; i++) lambda += derivada[i]*VECTOR(possivel_solucao)[i];
        lambda = (total != 0)? lambda/total : 0;
        for ( i = 0; i < BPR_PARAMETERS->N; i++){
            for (j = 0; j < BPR_PARAMETERS->N; j++){
                if(i != j) OD_i.MATRIZ[i][j] =OD_i.MATRIZ[i][j] - lambda*dZ_dOD[i][j];
                if(OD_i.MATRIZ[i][j] < 0) OD_i.MATRIZ[i][j] = 0;
            }
            //print_vetor(OD_i.MATRIZ[i],BPR_PARAMETERS->N,sizeof(int));
        }

        free(derivada);
        for (i = 0; i < BPR_PARAMETERS->L; i++) free(matrix_solution[i]);
        free(matrix_solution);
        //igraph_vector_print(&possivel_solucao);
        igraph_vector_destroy(&possivel_solucao);
        //printf("=====================================%.20f=============================\n",lambda);
        printf("%f\n",Z);
        if(lambda == 0) break;
        //if(df - lambda < 1e-2) break;
    }
    
}