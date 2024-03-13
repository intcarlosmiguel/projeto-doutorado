#pragma once

#define ALPHA 0.03
#define BETA 4
const double decresse = 0.9, min_step = 1E-8, max_step = 1.0, ftol = 0.2,X_TOLERANCIA = 0.001;
const int MAXIMO_ITERACOES = 1000;
struct PARAMETERS{
    igraph_vector_t capacidade;
    igraph_vector_t cost_time;
    int L;
};