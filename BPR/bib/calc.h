#pragma once

#include <define.h>
#include "igraph.h"

bool check_inconsistent_edges(const igraph_t *graph, igraph_vector_int_t *edge_id, int **edge_list) {
    igraph_eit_t edge_it;
    bool found_inconsistency = false;
    
    igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &edge_it);
    for (; !IGRAPH_EIT_END(edge_it); IGRAPH_EIT_NEXT(edge_it)) {
        igraph_integer_t eid = IGRAPH_EIT_GET(edge_it);
        igraph_integer_t from, to;
        igraph_edge(graph, eid, &from, &to);
        int index = VECTOR(*edge_id)[eid];
        int a = edge_list[index][0];
        int b = edge_list[index][1];
        
        if (a != from || b != to) {
            printf("Inconsistency found in edge %ld: bush(%ld->%ld) != edge_list(%d->%d)\n", 
                    eid, from, to, a, b);
            found_inconsistency = true;
        }
    }
    igraph_eit_destroy(&edge_it);
    return found_inconsistency;
}

void copy_bush(
    struct BUSH *bush, 
    struct BUSH *anterior_bush
) {
    igraph_copy(&bush->grafo, &anterior_bush->grafo);
    igraph_vector_update(&bush->flow, &anterior_bush->flow);
    igraph_vector_int_update(&bush->edge_id, &anterior_bush->edge_id);
    memcpy(bush->is_ingraph, anterior_bush->is_ingraph, bush->BPR_PARAMETERS.L * sizeof(bool));
    bush->n_alvos = anterior_bush->n_alvos;
    igraph_vector_update(&bush->BPR_PARAMETERS.capacidade, &anterior_bush->BPR_PARAMETERS.capacidade);
    igraph_vector_update(&bush->BPR_PARAMETERS.cost_time, &anterior_bush->BPR_PARAMETERS.cost_time);
    bush->BPR_PARAMETERS.L = anterior_bush->BPR_PARAMETERS.L;
    bush->BPR_PARAMETERS.N = anterior_bush->BPR_PARAMETERS.N;
    //printf("Copia com sucesso!\n");
    bush->steps = (double*) calloc(bush->n_alvos, sizeof(double));
    for (int i = 0; i < bush->n_alvos; i++) {
        bush->steps[i] = anterior_bush->steps[i];
    }
}

void init_path(struct min_max_bush *path, int N){
    igraph_vector_int_init(&path->min_parents, N);
    igraph_vector_int_init(&path->max_parents, N);
    igraph_vector_int_fill(&path->min_parents, -1); // -1 indica nenhum predecessor
    igraph_vector_int_fill(&path->max_parents, -1);
    igraph_vector_init(&path->derivate_dist_shortest_local, N);
    igraph_vector_fill(&path->derivate_dist_shortest_local, 0); // Inicializa com 0
    igraph_vector_init(&path->derivate_dist_longest_local, N);
    igraph_vector_fill(&path->derivate_dist_longest_local, 0); // Inicializa com 0
    igraph_vector_init(&path->dist_shortest_local, N);
    igraph_vector_fill(&path->dist_shortest_local, DBL_MAX);
    igraph_vector_init(&path->dist_longest_local, N);
    igraph_vector_fill(&path->dist_longest_local, -DBL_MAX); // "Infinito negativo"
}

void erase_path(struct min_max_bush *path){
    igraph_vector_int_destroy(&path->min_parents);
    igraph_vector_int_destroy(&path->max_parents);
    igraph_vector_destroy(&path->dist_shortest_local);
    igraph_vector_destroy(&path->dist_longest_local);
    igraph_vector_destroy(&path->derivate_dist_shortest_local);
    igraph_vector_destroy(&path->derivate_dist_longest_local);
    path->derivate_min_cost = 0.0;
    path->derivate_max_cost = 0.0;

}

void reverse_parents(igraph_vector_int_t *parent,int from, int to){
    igraph_vector_int_t reverse_parent;
    igraph_vector_int_init(&reverse_parent, 0);
    int alvo = to;
    while(alvo != from && alvo != -1) {
        igraph_integer_t prev = VECTOR(*parent)[alvo];
        igraph_vector_int_insert(&reverse_parent, 0, alvo); // Insere no início para reverter a ordem
        alvo = prev;
    }
    igraph_vector_int_insert(&reverse_parent, 0, from);
    igraph_vector_int_update(parent, &reverse_parent);
    igraph_vector_int_destroy(&reverse_parent);
}

void print_vetor(void* array,int N,int check){
    if(check == sizeof(int)){
        int* intArray = (int*)array;
        for (int i = 0; i < N; i++){
            if(i!=N-1) printf("%d ",intArray[i]);
            else printf("%d\n",intArray[i]);
        }
    }
    if(check == sizeof(double)){
        double* doubleArray = (double*)array;
        for (int i = 0; i < N; i++){
            if(i!=N-1) printf("%.2f ",doubleArray[i]);
            else printf("%.2f\n",doubleArray[i]);
        }
    }
}

int contarLinhasNoArquivo(const char *nomeArquivo) {
    FILE *arquivo = fopen(nomeArquivo, "r");
    if (!arquivo) {
        perror("Erro ao abrir o arquivo");
        exit(0);
    }
    
    int linhas = 0;
    char c;
    
    while ((c = fgetc(arquivo)) != EOF) {
        if (c == '\n') {
            linhas++;
        }
    }
    
    fclose(arquivo);
    return linhas + 1;
}

double** lerArquivo(const char *nomeArquivo, int nColunas,int* size) {
    int N = contarLinhasNoArquivo(nomeArquivo);
    double** data = (double**) malloc(N*sizeof(double*));
    
    FILE *arquivo;
    char buffer[1024];

    arquivo = fopen(nomeArquivo, "r");
    if (!arquivo) {
        perror("Erro ao abrir o arquivo");
        exit(0);
    }

    // Lê o arquivo linha por linha
    int linha = 0,i;
    while (fgets(buffer, 1024, arquivo) != NULL) {
        // Processa cada linha do arquivo aqui. Neste exemplo, vamos assumir que os dados são números inteiros.
        char *token = strtok(buffer, " "); // Supondo que os dados sejam separados por espaços. Ajuste o delimitador conforme necessário.
        data[linha] =(double*) malloc(nColunas*sizeof(double));
        for (i = 0; i < nColunas && token != NULL; i++) {
            // Converte o token (string) para o tipo de dado desejado, neste caso, int.
            sscanf(token, "%lf", &data[linha][i]);
            
            // Processa o dado da coluna aqui. Neste exemplo, estamos apenas imprimindo.
            //printf("Dado da coluna %d: %f\n", i + 1, data[linha][i]);

            // Avança para o próximo token (próxima coluna)
            token = strtok(NULL, " ");
        }
        linha++;
    }

    // Fecha o arquivo ao terminar de processar.
    fclose(arquivo);
    *size = N;
    return data;
}
int** init_parameters(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_int_t* edges, const char* nomeDoArquivo) {

    int i,e1,e2;
    double capacidade,tempo;

    int N = contarLinhasNoArquivo(nomeDoArquivo);

    igraph_vector_init(&BPR_PARAMETERS->capacidade, 0);
    igraph_vector_init(&BPR_PARAMETERS->cost_time, 0);

    int** edge_list = (int**)malloc(N*sizeof(int*));
    BPR_PARAMETERS->L = N;
    BPR_PARAMETERS->N = 0;
    FILE *arquivo;
    arquivo = fopen(nomeDoArquivo, "r");
    if (!arquivo) {
        perror("Erro ao abrir o arquivo");
        exit(0);
    }
    for (i = 0; i < N; i++){
        if(fscanf(arquivo, "%d %d %lf %lf\n", &e1, &e2, &capacidade, &tempo));
        edge_list[i] = (int*)calloc(2,sizeof(int));
        edge_list[i][0] = e1-1;
        edge_list[i][1] = e2-1;
        igraph_vector_int_push_back(edges,edge_list[i][0]);
        igraph_vector_int_push_back(edges,edge_list[i][1]);
        igraph_vector_push_back(&BPR_PARAMETERS->capacidade,capacidade);
        igraph_vector_push_back(&BPR_PARAMETERS->cost_time,tempo);
        if(BPR_PARAMETERS->N < edge_list[i][0]+1) BPR_PARAMETERS->N = edge_list[i][0]+1;


    }
    fclose(arquivo);

    return edge_list;
}