#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "igraph.h"
#include <limits.h> // For INT_MAX
#include <float.h>  // Para DBL_MAX

#define ALPHA 0.03
#define BETA 4
const double decresse = 0.9, min_step = 1E-8, max_step = 1.0, ftol = 0.2,X_TOLERANCIA = 0.0001;
const int MAXIMO_ITERACOES = 10000;
const int THREADS = 10;
const double EPSILON = 1E-6;

struct MATRIZ_OD{
    int** MATRIZ;
    int N_FONTES;
    int N_ALVOS;
    int** LIST;
    int** indexate;
    igraph_vector_int_t fontes;
    igraph_vector_int_t alvos;
};

struct PARAMETERS{
    igraph_vector_t capacidade;
    igraph_vector_t cost_time;
    int L;
    int N;
};
struct ElementOD{
    int fonte;
    igraph_vector_int_t alvos;
    igraph_vector_int_t volumes;
};
struct OD_MATRIX{
    struct ElementOD* Elementos;
    int size;
};
struct min_max_bush{
    double min_cost;
    double max_cost;
    double derivate;
    igraph_vector_int_t min_edges;
    igraph_vector_int_t max_edges;
    igraph_vector_t dist_shortest_local;
    igraph_vector_t dist_longest_local;
    double mu;
};

struct BUSH{
    igraph_vector_int_t edge_id;
    bool* is_ingraph; // Indica se o nó está no grafo
    igraph_t Grafo; // Grafo da bush
    int n_alvos;
    struct min_max_bush *paths;
    double* steps;
    double* flow_per_alvo; // Fluxo por alvo
};