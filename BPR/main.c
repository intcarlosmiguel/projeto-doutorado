#include "bib/simulate.h"
#include "igraph.h"
// 4.077634084490802
// 60
double gera_lei_potencia(double alpha, double xmin) {
    // Gera um número aleatório entre 0 e 1
    double r = (double)rand() / RAND_MAX;

    // Aplica a transformação para gerar um valor que siga a lei de potência
    double valor = xmin * pow(1 - r, -1.0 / (alpha - 1));

    return valor;
}

int main(int argc,char *argv[ ]){
    char *arquivo1 = argv[1];
    char *arquivo2 = argv[2];
    /* igraph_t Grafo;
    igraph_vector_int_t size;
    igraph_vector_int_init(&size,2);
    VECTOR(size)[0] = 10;
    VECTOR(size)[1] = 10;
    igraph_square_lattice(&Grafo , &size, 1,0, 0,0); */
    simulate_example(arquivo1, arquivo2);
    //int arquivo = atoi(argv[1]);
    //int seed = atoi(argv[2]);
    //simulate(arquivo,seed);
}