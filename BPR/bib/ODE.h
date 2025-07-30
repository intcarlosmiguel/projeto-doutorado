#pragma once
#include "igraph.h"
#include "network.h"
#include "define.h"
double F2(
    igraph_vector_t *time,
    igraph_vector_t *time_hat,
    int L
) {
    // Compute the objective function value
    double objective_value = 0.0;
    for (int i = 0; i <L; i++) {
        objective_value += pow(VECTOR(*time_hat)[i] - VECTOR(*time)[i], 2)/2;
    }
    printf("Objective Function Value: %.20f\n", objective_value);
    return objective_value;
}
double F1(
    struct OD_MATRIX* OD
) {
    // Compute the objective function value
    double objective_value = 0.0, volume = 0.0;
    int L;
    for (int i = 0; i < OD->size; i++) {
        L = igraph_vector_int_size(&OD->Elementos[i].alvos);
        for (int j = 0; j < L; j++) {
            volume = VECTOR(OD->Elementos[i].volumes)[j];
            if (volume <= 0) continue; // Skip if volume is zero or negative
            objective_value += volume * log10(volume);
        }
    }
    printf("Objective Function Value: %.20f\n", objective_value);
    return objective_value;
}

double ObjectiveFunction(
    igraph_vector_t *time,
    igraph_vector_t *time_hat,
    struct OD_MATRIX* OD,
    int L,
    double lambda
) {
    // Compute the objective function value
    double objective_value = 0.0;
    objective_value += F1(OD);
    objective_value += lambda*F2(time, time_hat, L);
    printf("Objective Function Value: %.20f\n", objective_value);
    return objective_value;
}
double calc_probability(
    struct BUSH* bush,
    int alvo,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *flow,
    igraph_vector_t *flow_hat
){
    double* prob = (double*) calloc(BPR_PARAMETERS->L, sizeof(double));
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);
    igraph_vector_int_t nodes;
    igraph_vector_int_init(&nodes, 0);
    igraph_vector_int_push_back(&nodes, alvo);
    bool * visited_nodes = (bool*) calloc(BPR_PARAMETERS->N, sizeof(bool));
    visited_nodes[alvo] = true;
    int n = 1,node,i,j,edge,from;
    for(i = 0; i < n; i++) {
        node = VECTOR(nodes)[i];
        igraph_incident(&bush->Grafo, &incident_edges, node, IGRAPH_IN);
        for(j = 0; j < igraph_vector_int_size(&incident_edges); j++) {
            edge = VECTOR(incident_edges)[j];
            from = IGRAPH_FROM(&bush->Grafo, edge);
            if(!visited_nodes[from]) {
                igraph_vector_int_push_back(&nodes, from);
                visited_nodes[from] = true;
                n++;
            }
        }
    }
    double* kappa = (double*) calloc(BPR_PARAMETERS->N, sizeof(double));
    kappa[alvo] = 1.0; // Set kappa for the target node
    double total_flow = 0.0;
    for(i = 0; i < n; i++) {
        total_flow = 0.0;
        node = VECTOR(nodes)[i];
        igraph_incident(&bush->Grafo, &incident_edges, node, IGRAPH_IN);
        int tamanho = igraph_vector_int_size(&incident_edges);
        igraph_vector_t phi;
        igraph_vector_init(&phi, tamanho);
        for(j = 0; j < tamanho; j++) {
            edge = VECTOR(incident_edges)[j];
            total_flow += VECTOR(bush->flow_per_origin)[VECTOR(bush->edge_id)[edge]];
        }
        for(j = 0; j < tamanho; j++) {
            edge = VECTOR(incident_edges)[j];
            VECTOR(phi)[j] = VECTOR(bush->flow_per_origin)[VECTOR(bush->edge_id)[edge]] / total_flow;
            if(VECTOR(phi)[j] > 0){
                prob[edge] = kappa[node] * VECTOR(phi)[j];
                from = IGRAPH_FROM(&bush->Grafo, edge);
                kappa[from] += prob[edge]; // Update kappa for the source node
            }
        }
        igraph_vector_destroy(&phi);
    }
    double d = 0;
    for(i = 0; i < BPR_PARAMETERS->L; i++) {
        d += prob[i]*(VECTOR(*flow)[i] - VECTOR(*flow_hat)[i]);
    }
    for(i = 0; i < BPR_PARAMETERS->L; i++) {
        if(prob[i] > 0)printf("(%ld,%ld) - %f\n", 
               IGRAPH_FROM(&bush->Grafo, i) + 1, 
               IGRAPH_TO(&bush->Grafo, i) + 1, 
               prob[i]);
    }
    free(prob);
    free(visited_nodes);
    free(kappa);
    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_int_destroy(&nodes);
    printf("d: %.20f\n", d);
    return d;
    

}
