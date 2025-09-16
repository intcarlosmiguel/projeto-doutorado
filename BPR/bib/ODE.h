#pragma once
#include "igraph.h"
#include "network.h"
#include "define.h"
#include <lbfgs.h>

double F2(
    igraph_vector_t *flow,
    igraph_vector_t *flow_hat,
    int L
) {
    // Compute the objective function value
    double objective_value = 0.0;
    for (int i = 0; i <L; i++) {
        objective_value += pow(VECTOR(*flow_hat)[i] - VECTOR(*flow)[i], 2)/2;
    }
    //printf("Objective Function Value: %.20f\n", objective_value);
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
    //printf("Objective Function Value: %.20f\n", objective_value);
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
    //printf("Objective Function Value: %.20f\n", objective_value);
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
        if(prob[i] > 0)printf("(%ld,%ld) - %f\n", IGRAPH_FROM(&bush->Grafo, i) + 1, IGRAPH_TO(&bush->Grafo, i) + 1, prob[i]);
    }
    free(prob);
    free(visited_nodes);
    free(kappa);
    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_int_destroy(&nodes);
    printf("d: %.20f\n", d);
    return d;
    

}


double OD_gradient(
    struct OD_MATRIX* OD,
    struct BUSH* bush,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *flow,
    igraph_vector_t *flow_hat,
    double lambda,
    lbfgsfloatval_t *gradiente
){

    int size = OD->size;
    double norma = 0.0;
    int value = 0;
    for(int i = 0; i < size; i++) {
        int L = igraph_vector_int_size(&OD->Elementos[i].alvos);
        if(L == 0) continue;
        for(int j = 0; j < L; j++) {
            int alvo = VECTOR(OD->Elementos[i].alvos)[j];
            double P = calc_probability(&bush[i], alvo, BPR_PARAMETERS, flow, flow_hat);
            double volume = VECTOR(OD->Elementos[i].volumes)[j];

            if(volume <= 1e-6) gradiente[value] = 0.0;
            else gradiente[value] = 1.0 + log10(volume) + lambda*P;
            value++;

        }
    }
    norma = sqrt(norma);
    double obj = ObjectiveFunction(flow, flow_hat, OD, BPR_PARAMETERS->L, lambda);
    double RMSE = F2(flow, flow_hat, BPR_PARAMETERS->L);
    printf("Updated Objective Function Value: %.20f\n", obj);
    printf("Norm of the gradient: %.20f\n", norma);
    printf("RMSE: %.20f\n", RMSE);
    printf("OD volumes updated.\n");
    return obj;
}

void read_flow_hat(const char* filename, igraph_vector_t* flow_hat) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    igraph_vector_init(flow_hat, 0);
    double value;
    while (fscanf(file, "%lf", &value) == 1) {
        igraph_vector_push_back(flow_hat, value);
    }

    fclose(file);
}

static lbfgsfloatval_t evaluate(
    void *instance, // Esta é a "caixa selada"
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
) {
    OPTIMIZATION_DATA* DATA = (OPTIMIZATION_DATA*) instance;
    igraph_t Grafo = *(DATA->graph);
    double lambda = DATA->lambda;
    struct OD_MATRIX* OD = DATA->od_matrix;
    struct PARAMETERS* BPR_PARAMETERS = DATA->bpr_params;
    igraph_vector_t* flow_hat = DATA->flow_hat;

    igraph_vector_t flow;
    igraph_vector_init(&flow, 0);
    struct BUSH* bushes = Bushes(&Grafo, &OD, &BPR_PARAMETERS, &flow);
    double obj = OD_gradient(OD, bushes, BPR_PARAMETERS, &flow, flow_hat, lambda, g);

    igraph_vector_destroy(&flow);
    igraph_destroy(&Grafo);

    return obj;
}


void OD_ESTIMATION_LBFGS(struct OD_MATRIX* OD, const char* file_edges, const char* file_od, const char* file_flow_hat, double lambda) {

    const int N = OD->all_elements; // Number of variables (size of x)
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(N);

    int value = 0;
    for(int i = 0; i < OD->size; i++) {
        int L = igraph_vector_int_size(&OD->Elementos[i].alvos);
        for(int j = 0; j < L; j++) {
            x[value] = VECTOR(OD->Elementos[i].volumes)[j];
            value++;
        }
    }

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.m = 10;
    // Exigir uma convergência mais rigorosa, com norma do gradiente menor
    param.epsilon = 1.0e-8;
    
    // Definir um limite máximo de 500 iterações como salvaguarda
    param.max_iterations = 500;

    struct PARAMETERS BPR_PARAMETERS;
    struct OD_MATRIX OD_MATRIX;
    igraph_vector_t flow_hat;
    
    igraph_vector_int_t edges;

    read_flow_hat(file_flow_hat, &flow_hat);
    igraph_vector_int_init(&edges, 0);
    int** edge_list = init_parameters(&BPR_PARAMETERS,&edges,file_edges);

    load_OD_from_file(file_od, &OD_MATRIX);
    print_OD_matrix(&OD_MATRIX);

    igraph_t Grafo;
    igraph_empty(&Grafo, BPR_PARAMETERS.N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    OPTIMIZATION_DATA DATA = {&Grafo, lambda, &OD_MATRIX, &BPR_PARAMETERS, &flow_hat};

    // 4. CHAMADA DO LBFGS
    int ret = 0;
    // Passamos o endereço do nosso rastreador como o argumento 'instance'
    //int ret = lbfgs(N, x, &fx, evaluate, progress, &tracker, &param);

    // Reporta o resultado final
    printf("========================================\n");
    printf("Status de retorno do L-BFGS: %d (%s)\n", ret, lbfgs_strerror(ret));
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);

    // 5. LIMPEZA
    // Libera a memória alocada para o rastreador e para as variáveis
    lbfgs_free(x);

    return 0;
}