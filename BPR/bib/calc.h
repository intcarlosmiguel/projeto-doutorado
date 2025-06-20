#pragma once

#include <define.h>
#include "igraph.h"
void print_flow(
    igraph_vector_t *flow,
    igraph_t *Grafo,
    char* filename
) {
    if (igraph_ecount(Grafo) != igraph_vector_size(flow)) {
        printf("Error: Number of edges (%ld) does not match flow vector size (%ld)\n", 
               igraph_ecount(Grafo), igraph_vector_size(flow));
        exit(1);
    }
    if (filename == NULL) {
        printf("Flow values:\n");
        for (long i = 0; i < igraph_ecount(Grafo); i++) {
            igraph_integer_t from, to;
            igraph_edge(Grafo, i, &from, &to);
            printf("%ld %ld %f\n", from, to, VECTOR(*flow)[i]);
        }
        printf("\n");
        return;
    }
    else{
        FILE *file = fopen(filename, "w");
        for (long i = 0; i < igraph_ecount(Grafo); i++) {
            igraph_integer_t from, to;
            igraph_edge(Grafo, i, &from, &to);
            fprintf(file,"%ld %ld %f\n", from, to, VECTOR(*flow)[i]);
        }
        fclose(file);

    }
}
/* void copy_bush(
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
} */

void init_path(struct min_max_bush *path, int N){
    igraph_vector_int_init(&path->min_edges, 0);
    igraph_vector_int_init(&path->max_edges, 0);
    
    igraph_vector_init(&path->dist_shortest_local, N);
    igraph_vector_fill(&path->dist_shortest_local, DBL_MAX);
    igraph_vector_init(&path->dist_longest_local, N);
    igraph_vector_fill(&path->dist_longest_local, -DBL_MAX); // "Infinito negativo"
}

void erase_path(struct min_max_bush *path){
    igraph_vector_int_destroy(&path->min_edges);
    igraph_vector_int_destroy(&path->max_edges);
    igraph_vector_destroy(&path->dist_shortest_local);
    igraph_vector_destroy(&path->dist_longest_local);
    path->derivate = 0.0;

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