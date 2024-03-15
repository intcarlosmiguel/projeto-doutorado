#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
            if(i!=N-1) printf("%f ",doubleArray[i]);
            else printf("%f\n",doubleArray[i]);
        }
    }
}

int contarLinhasNoArquivo(const char *nomeArquivo) {
    FILE *arquivo;
    char buffer[1024]; // Buffer para armazenar cada linha lida do arquivo
    int contadorDeLinhas = 0; // Contador para o número de linhas

    // Abre o arquivo para leitura
    arquivo = fopen(nomeArquivo, "r");
    if (!arquivo) {
        perror("Erro ao abrir o arquivo");
        return -1; // Retorna -1 em caso de erro ao abrir o arquivo
    }

    // Lê o arquivo linha por linha
    while (fgets(buffer, 1024, arquivo) != NULL) {
        contadorDeLinhas++; // Incrementa o contador para cada linha lida
    }

    // Fecha o arquivo
    fclose(arquivo);

    // Retorna o número total de linhas
    return contadorDeLinhas;
}

double** lerArquivo(const char *nomeArquivo, int nColunas,int* size) {
    int N = contarLinhasNoArquivo(nomeArquivo);
    double** data = (double**) malloc(N*sizeof(double*));
    
    FILE *arquivo;
    char buffer[1024]; // Buffer para armazenar cada linha do arquivo. Ajuste o tamanho conforme necessário.

    // Tentativa de abrir o arquivo para leitura. Se falhar, exibe uma mensagem de erro e termina a função.
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