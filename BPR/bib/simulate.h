#pragma once
#include "igraph.h"
#include "define.h"
#include "calc.h"
#include "network.h"
#include "example.h"
#include "BPR.h"
#include "mtwister.h"
#include "dial.h"
//#include "PSO.h"
//#include "GM.h"

void random_OD(struct MATRIZ_OD *OD_i,int N,int seed){
    init_genrand64(seed);
    igraph_vector_int_init(&OD_i->fontes, 0);
    igraph_vector_int_init(&OD_i->alvos, 0);
    igraph_vector_int_range(&OD_i->fontes,0,N);
    igraph_vector_int_range(&OD_i->alvos,0,N);
    OD_i->N_ALVOS = N;
    OD_i->N_FONTES = N;
    OD_i->MATRIZ = (int**) malloc(N*sizeof(int*));
    OD_i->LIST = (int**) malloc((N)*(N-1)*sizeof(int*));
    int c = 0;
    int total = 4000;
    int i,j;
    for ( i = 0; i < N; i++)OD_i->MATRIZ[i] = (int*) calloc(N,sizeof(int));
    for ( i = 0; i < N; i++){
        
        for ( j = 0; j < N; j++){
            if(i !=j){
                OD_i->LIST[c] = (int*) malloc(2*sizeof(int));
                OD_i->LIST[c][0] = i;
                OD_i->LIST[c][1] = j;
                c++;
            }
        }
    }
    while(total > 0 ){
        i = genrand64_real1()*N;
        j = genrand64_real1()*N;
        if(i == j) continue;
        if(OD_i->MATRIZ[i][j] == 0 )OD_i->MATRIZ[i][j] = 5*genrand64_real1();
        else continue;
        if(total - OD_i->MATRIZ[i][j] < 0)OD_i->MATRIZ[i][j] = total;
        total = total - OD_i->MATRIZ[i][j];
    }
}

/* void obstructed(struct PARAMETERS* BPR_PARAMETERS,int** edge_list,igraph_vector_t* solucao,struct MATRIZ_OD* OD,igraph_vector_int_t* edges,char* nomeDoArquivo){
    int i,j;
    struct PARAMETERS BPR_PARAMETERS_obstructed;

    BPR_PARAMETERS_obstructed.N = BPR_PARAMETERS->N;
    BPR_PARAMETERS_obstructed.L = BPR_PARAMETERS->L;

    igraph_vector_init(&BPR_PARAMETERS_obstructed.capacidade, BPR_PARAMETERS_obstructed.L);
    igraph_vector_init(&BPR_PARAMETERS_obstructed.cost_time, BPR_PARAMETERS_obstructed.L);

    for ( i = 0; i < BPR_PARAMETERS->L; i++){
        VECTOR(BPR_PARAMETERS_obstructed.capacidade)[i] = VECTOR(BPR_PARAMETERS->capacidade)[i];
        VECTOR(BPR_PARAMETERS_obstructed.cost_time)[i] = VECTOR(BPR_PARAMETERS->cost_time)[i];
    }

    double** resultados = malloc(BPR_PARAMETERS->L* sizeof(double*));
    for (i = 0; i < BPR_PARAMETERS->L; i++) resultados[i] = calloc(BPR_PARAMETERS->L, sizeof(double)); // Aloca espaço para cada string

    FILE *file;
    file = fopen(nomeDoArquivo,"w");
    for (i = 0; i < BPR_PARAMETERS->L; i++){

        VECTOR(BPR_PARAMETERS_obstructed.capacidade)[i] = 1;

        igraph_t Grafo_obstructed;
        igraph_empty(&Grafo_obstructed, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
        igraph_add_edges(&Grafo_obstructed, edges, NULL);
        igraph_vector_t solucao_obstructed;

        optimize(&BPR_PARAMETERS_obstructed,edge_list,OD,&Grafo_obstructed,&solucao_obstructed);

        for (j = 0; j < BPR_PARAMETERS->L; j++){
            if(i == j)resultados[i][j] = 0;
            else resultados[i][j] =  VECTOR(*solucao)[j] - VECTOR(solucao_obstructed)[j];
            fprintf(file, "%.2f ", resultados[i][j]);
        }
        fprintf(file, "\n");
        igraph_vector_destroy(&solucao_obstructed);
        igraph_destroy(&Grafo_obstructed);
        printf("%d/%d\n",i+1,BPR_PARAMETERS->L);
        //if(i == 1) break;
        VECTOR(BPR_PARAMETERS_obstructed.capacidade)[i] = VECTOR(BPR_PARAMETERS->capacidade)[i];
    }
    fclose(file);
    igraph_vector_destroy(&BPR_PARAMETERS_obstructed.capacidade);
    igraph_vector_destroy(&BPR_PARAMETERS_obstructed.cost_time);
    free(resultados);

} */

void init_simulate(struct PARAMETERS* BPR_PARAMETERS,int** edge_list,struct OD_MATRIX *OD,igraph_vector_t *solucao,igraph_vector_int_t * edges){
    igraph_t Grafo;
    igraph_empty(&Grafo, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, edges, NULL);
    Bushes(&Grafo,OD,edge_list,BPR_PARAMETERS,solucao);
    //optimize(BPR_PARAMETERS,edge_list,OD,&Grafo,solucao);
    igraph_destroy(&Grafo);
}

void load_OD_from_file(const char* filename, struct OD_MATRIX* OD_MATRIX) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        return;
    }
    int current_source = -1;

    int source, target;
    double flow;
    OD_MATRIX->size = 0;
    OD_MATRIX->Elementos = malloc(sizeof(struct ElementOD) * 0);
    while (fscanf(file, "%d %d %lf", &source, &target, &flow) == 3) {
        if (current_source != source - 1) {
            OD_MATRIX->Elementos = realloc(OD_MATRIX->Elementos, sizeof(struct ElementOD) * (OD_MATRIX->size + 1));
            OD_MATRIX->Elementos[OD_MATRIX->size].fonte = source - 1;
            igraph_vector_int_init(&OD_MATRIX->Elementos[OD_MATRIX->size].alvos, 0);
            igraph_vector_int_init(&OD_MATRIX->Elementos[OD_MATRIX->size].volumes, 0);
            igraph_vector_int_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size].alvos, target - 1);
            igraph_vector_int_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size].volumes, flow);
            OD_MATRIX->size++;
            current_source = source - 1;
        }
        else{
            igraph_vector_int_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size - 1].alvos, target - 1);
            igraph_vector_int_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size - 1].volumes, flow);
        }
    }

    fclose(file);
}

void print_OD_matrix(struct OD_MATRIX* OD_MATRIX) {
    printf("Origin-Destination Matrix:\n");
    for (int i = 0; i < OD_MATRIX->size; i++) {
        printf("Origin %d -> Destinations: ", OD_MATRIX->Elementos[i].fonte);
        for (int j = 0; j < igraph_vector_int_size(&OD_MATRIX->Elementos[i].alvos); j++) {
            printf("(D:%ld, V:%ld) ", VECTOR(OD_MATRIX->Elementos[i].alvos)[j], VECTOR(OD_MATRIX->Elementos[i].volumes)[j]);
        }
        printf("\n");
    }
}

void simulate_example(const char* nomeDoArquivo,const char* nomeDoArquivo2){

    struct PARAMETERS BPR_PARAMETERS;
    struct OD_MATRIX OD_MATRIX;


    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    int** edge_list = init_parameters(&BPR_PARAMETERS,&edges,nomeDoArquivo);

    load_OD_from_file(nomeDoArquivo2, &OD_MATRIX);
    print_OD_matrix(&OD_MATRIX);
    igraph_vector_t solucao;
    init_simulate(&BPR_PARAMETERS,edge_list,&OD_MATRIX,&solucao,&edges);
}

/* void simulate(int arquivo,int seed){
    omp_set_num_threads(THREADS);
    char nome_arquivo_edges[100];
    char nome_arquivo_resultados[100];
    char nome_arquivo_solution[100];
    sprintf(nome_arquivo_edges, "./file/edges_%d.txt", arquivo);
    sprintf(nome_arquivo_resultados, "./output/resultados_%d.dat", arquivo);
    sprintf(nome_arquivo_solution, "./output/solution_%d.dat", arquivo);

    struct PARAMETERS BPR_PARAMETERS;
    struct MATRIZ_OD OD;
   
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    int** edge_list = init_parameters(&BPR_PARAMETERS,&edges,nome_arquivo_edges);

    random_OD(&OD,BPR_PARAMETERS.N,seed);

    igraph_vector_t solucao;

    init_simulate(&BPR_PARAMETERS,edge_list,&OD,&solucao,&edges);
    FILE *file;
    file = fopen(nome_arquivo_solution,"w");
    for ( int i = 0; i < BPR_PARAMETERS.L; i++)fprintf(file,"%f ",VECTOR(solucao)[i]);
    fclose(file);
    obstructed(&BPR_PARAMETERS,edge_list,&solucao,&OD,&edges,nome_arquivo_resultados);
    //igraph_vector_print(&solucao);

    
    
    load_MATOD(&OD,true);

    igraph_vector_t solucao;
    double** matrix_solution = (double**)malloc(BPR_PARAMETERS.L*sizeof(double*));
    for (int i = 0; i < BPR_PARAMETERS.L; i++)matrix_solution[i] = (double*)calloc(OD.N_ALVOS*(OD.N_FONTES-1),sizeof(double));
    optimize(&BPR_PARAMETERS,edge_list,&OD,&Grafo,&solucao,matrix_solution);
    //igraph_vector_print(&solucao);
    //caminhos(matrix_solution,&BPR_PARAMETERS,edge_list,0,6,8);
    //GM(&OD,&BPR_PARAMETERS,edge_list,&Grafo,&solucao);
    //igraph_vector_print(&solucao);

    //PSO(50,500,&BPR_PARAMETERS,edge_list,&Grafo,&fontes, &alvos,&solucao);

    igraph_vector_destroy(&BPR_PARAMETERS.cost_time);
    igraph_vector_destroy(&BPR_PARAMETERS.capacidade);
    igraph_vector_int_destroy(&OD.fontes);
    igraph_vector_int_destroy(&OD.alvos);
    igraph_vector_int_destroy(&edges);
    for (int i = 0; i < BPR_PARAMETERS.N; i++) free(OD.MATRIZ[i]);
    free(OD.MATRIZ);
    igraph_destroy(&Grafo);
} */