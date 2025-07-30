
#pragma once

#include "igraph.h"
#include "network.h"
#include "ODE.h"
void check_solution(
    struct BUSH* bushes,
    igraph_vector_t* solucao,
    struct PARAMETERS* BPR_PARAMETERS,
    int size
){
    for (int i = 0; i <size; i++){
        for (int k = 0; k < BPR_PARAMETERS->L; k++) {
            if (VECTOR(bushes[i].flow_per_origin)[k] > VECTOR(*solucao)[k]) {
                //printf("Error: Bush %d has flow_per_alvo[%d] = %.20f greater than solucao[%d] = %.20f\n",i+1, k, VECTOR(bushes[i].flow_per_origin)[k], k, VECTOR(*solucao)[k]);
                VECTOR(*solucao)[k] = VECTOR(bushes[i].flow_per_origin)[k]; // Ajusta o fluxo da solução para o fluxo da bush
            }
            if(VECTOR(bushes[i].flow_per_origin)[k] < 0){
                //printf("Error: Bush %d has negative flow_per_alvo[%d] = %.20f\n",i+1, k, VECTOR(bushes[i].flow_per_origin)[k]);
                VECTOR(*solucao)[k] = 0.0; // Ajusta o fluxo da solução para 0
            }
        }
    }
}

int get_princial_nodes(
    int *node_M,
    struct min_max_bush *paths,
    igraph_t *Grafo,
    int fonte,
    bool* visited_nodes
) {
    while(VECTOR(paths->max_edges)[*node_M] == VECTOR(paths->min_edges)[*node_M]){
        *node_M = IGRAPH_FROM(Grafo, VECTOR(paths->max_edges)[*node_M]);
        if(*node_M == fonte)return fonte;
    }
    
    int N = IGRAPH_FROM(Grafo, VECTOR(paths->max_edges)[*node_M]);
    //printf("Nó M: %d, Nó N: %d\n", *node_M, N);
    while(!visited_nodes[N]) N = IGRAPH_FROM(Grafo, VECTOR(paths->max_edges)[N]);
    return N;
}

void init_bush(
    struct BUSH *bush, 
    igraph_t *Grafo, 
    int i,
    int** edge_list,
    struct PARAMETERS* BPR_PARAMETERS,
    struct OD_MATRIX* OD,
    igraph_vector_t *flow
) {
    igraph_vector_int_t inbound;
    
    igraph_vector_int_init(&inbound, 0);
    
    int fonte = OD->Elementos[i].fonte,j,k,id;
    //printf("Fonte: %d!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", fonte);
    
    //BPR(&time, BPR_PARAMETERS, flow);
    //igraph_vector_print(&BPR_PARAMETERS->cost_time);
    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL,NULL,fonte,igraph_vss_all(),&BPR_PARAMETERS->cost_time, IGRAPH_OUT,NULL,&inbound
    );
    int alvo,edge_id,from,to;
    double tempo;
    double* total_time = (double*) calloc(BPR_PARAMETERS->N , sizeof(double));
    igraph_vector_int_init(&bush->edge_id, 0);
    bush->is_ingraph = (bool*) calloc(BPR_PARAMETERS->L, sizeof(bool)); // Inicializa o vetor de arestas no grafo
    int index,target,antecessor;
    igraph_vector_int_t edges_vec;
    igraph_vector_int_init(&edges_vec, 0);

    for ( j = 0; j< BPR_PARAMETERS->N; j++){
        if(j == fonte) continue; // Pula a fonte
        alvo = j;
        tempo = 0;
        while(alvo!=fonte){
            edge_id = VECTOR(inbound)[alvo];
            tempo += VECTOR(BPR_PARAMETERS->cost_time)[edge_id];
            alvo = IGRAPH_FROM(Grafo, edge_id);            
        }
        total_time[j] = tempo;
    }
    for (j = 0; j < BPR_PARAMETERS->L; j++) {
        from = edge_list[j][0];
        to = edge_list[j][1];
        if (total_time[from] < total_time[to]){
            igraph_vector_int_push_back(&edges_vec, j);
            bush->is_ingraph[j] = true; // Marca a aresta como parte do grafo da bush
            //printf("Aresta (%ld,%ld) adicionada ao grafo da bush\n", IGRAPH_FROM(Grafo, j)+1, IGRAPH_TO(Grafo, j)+1);
        }
    }
    //exit(0); // Encerra o programa após inicializar a bush
    //printf("Size of bush->flow: %ld\n", igraph_vector_size(&bush->flow));
    //igraph_vector_int_print(&edges_vec);
    igraph_es_t edges;
    igraph_es_vector(&edges,&edges_vec);
    
    igraph_subgraph_from_edges(
        Grafo, &bush->Grafo, edges, false
    );
    //printf("Number of edges in Bush: %ld\n", igraph_ecount(&bush->grafo));
    init_path(&bush->paths, BPR_PARAMETERS->N);
    for(j = 0; j < bush->n_alvos; j++) {
        antecessor = VECTOR(OD->Elementos[i].alvos)[j];
        while(antecessor != fonte) {
            index = VECTOR(inbound)[antecessor];
            VECTOR(*flow)[index] += VECTOR(OD->Elementos[i].volumes)[j];
            //printf("Adicionando fluxo %ld na aresta %d (%ld,%ld)\n", VECTOR(OD->Elementos[i].volumes)[j],index, IGRAPH_FROM(Grafo, index)+1, IGRAPH_TO(Grafo, index)+1);
            VECTOR(bush->flow_per_origin)[index] += VECTOR(OD->Elementos[i].volumes)[j];
            antecessor = IGRAPH_FROM(Grafo, index);

        }
    }

    //exit(0); // Encerra o programa após inicializar a bush

    char filename[100];
    sprintf(filename, "bush_%d.txt", fonte);
    FILE* file = fopen(filename, "w");
    for (long e = 0; e < igraph_ecount(&bush->Grafo); e++) {
        igraph_integer_t from, to;
        igraph_edge(&bush->Grafo, e, &from, &to);
        index = find_id_log(from, to, edge_list, BPR_PARAMETERS->L);
        igraph_vector_int_push_back(&bush->edge_id, index);
        //printf("Aresta %d (%ld,%ld) - %f %f\n",index, from, to, VECTOR(*flow)[index],bush->flow_per_alvo[index]);
        fprintf(file,"%ld %ld %f\n", IGRAPH_FROM(Grafo,index), IGRAPH_TO(Grafo,index), VECTOR(*flow)[index]);
    }
    
    fclose(file);
    igraph_vector_int_destroy(&inbound);
    //exit(0); // Encerra o programa após inicializar a bush
}


void find_dag_shortest_longest_costs_and_parents(
    struct BUSH *bush,
    igraph_t *Grafo,
    igraph_integer_t from,
    struct PARAMETERS *BPR_PARAMETERS,
    igraph_vector_t *total_flow
) {
    //printf("\n=================================================================================================\n");
    //print_flow(total_flow, Grafo, NULL); // Imprime o fluxo total
    //printf("Encontrando caminhos entre %ld ➡️  %ld\n", from, to);
    int N = igraph_vcount(Grafo),j;
    // Garante que os vetores de saída estejam dimensionados e inicializados
    // Vetores locais para armazenar distâncias para todos os nós
    igraph_vector_fill(&bush->paths.dist_shortest_local, DBL_MAX);
    igraph_vector_fill(&bush->paths.dist_longest_local, -DBL_MAX); // "Infinito negativo"
    VECTOR(bush->paths.dist_shortest_local)[from] = 0.0;
    VECTOR(bush->paths.dist_longest_local)[from] = 0.0;
    VECTOR(bush->paths.min_edges)[from] = -1; // Inicializa o vetor de arestas mínimas
    VECTOR(bush->paths.max_edges)[from] = -1; // Inicializa o vetor de arestas máximas
    
    igraph_vector_int_t topological_order;
    igraph_vector_int_init(&topological_order, 0); // Inicializa para evitar problemas no destroy

    int L = igraph_ecount(&bush->Grafo);
    //printf("Número de arestas no grafo da bush: %d\n", L);
    if (igraph_topological_sorting(&bush->Grafo, &topological_order, IGRAPH_OUT) != IGRAPH_SUCCESS) {
        printf("Erro ao ordenar topologicamente o grafo.\n");
        exit(0); // Falha ao ordenar topologicamente, encerra o programa
    }
    
    
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);
    igraph_vector_t flow_bush;
    igraph_vector_init(&flow_bush, 0);
    igraph_vector_t time;
    igraph_vector_init(&time, L);
    igraph_vector_t dtime;
    igraph_vector_init(&dtime, BPR_PARAMETERS->L);
    struct PARAMETERS bush_PARAMETERS;
    igraph_vector_init(&bush_PARAMETERS.capacidade, 0);
    igraph_vector_init(&bush_PARAMETERS.cost_time, 0);
    bush_PARAMETERS.L = L;
    bush_PARAMETERS.N = N;

    int i,index,id,edge_id;
    double fluxo,weight,dfluxo;
    for(i = 0; i < L; i++) {
        index = VECTOR(bush->edge_id)[i];
        igraph_vector_push_back(&flow_bush, VECTOR(*total_flow)[index]);
        //printf("Aresta %d (%ld,%ld) - %f,%f\n", index, IGRAPH_FROM(&bush->Grafo, i), IGRAPH_TO(&bush->Grafo, i), VECTOR(*total_flow)[index], VECTOR(bush->flow_per_alvo)[i]);
        igraph_vector_push_back(&bush_PARAMETERS.capacidade, VECTOR(BPR_PARAMETERS->capacidade)[index]);
        igraph_vector_push_back(&bush_PARAMETERS.cost_time, VECTOR(BPR_PARAMETERS->cost_time)[index]);

    }
    BPR(&time, &bush_PARAMETERS, &flow_bush); // Calcula o tempo de cada aresta
    BPR_derivate(&dtime,BPR_PARAMETERS, total_flow); // Calcula o gradiente de cada aresta
    //printf("Custo do tempo: \n");
    //igraph_vector_print(&time); // Imprime o tempo de cada aresta
    igraph_integer_t eid,u,v_node;

    //printf("Fluxo do tempo: \n");
    //print_flow(&flow_bush, &bush->Grafo, NULL); // Imprime o fluxo do tempo
    for (i = 0; i < igraph_vector_int_size(&topological_order); ++i) {
        u = VECTOR(topological_order)[i];
        // Se 'u' não é alcançável a partir de 'from' (para menor caminho),
        // não pode contribuir para caminhos subsequentes.
        if (VECTOR(bush->paths.dist_shortest_local)[u] == DBL_MAX) continue;
        
        igraph_incident(&bush->Grafo, &incident_edges, u, IGRAPH_OUT);
        
        for (j = 0; j < igraph_vector_int_size(&incident_edges); ++j) {

            eid = VECTOR(incident_edges)[j];
            id = VECTOR(bush->edge_id)[eid]; // Acessa o id da aresta
            v_node = IGRAPH_TO(&bush->Grafo, eid);
            weight = VECTOR(time)[eid]; // Acessa o peso da aresta
            
            
            fluxo = VECTOR(bush->flow_per_origin)[id]; // Acessa o fluxo da aresta
            dfluxo = VECTOR(dtime)[id]; // Acessa o gradiente da aresta

            if (VECTOR(bush->paths.dist_shortest_local)[u] + weight < VECTOR(bush->paths.dist_shortest_local)[v_node]) {

                VECTOR(bush->paths.dist_shortest_local)[v_node] = VECTOR(bush->paths.dist_shortest_local)[u] + weight;
                VECTOR(bush->paths.min_edges)[v_node] = id;
                VECTOR(bush->paths.derivate_shortest)[v_node] = dfluxo + VECTOR(bush->paths.derivate_shortest)[u]; // Atualiza o gradiente do caminho mínimo

            }

            // Relaxamento para o maior caminho
            // Só estender de 'u' se 'u' já tem um caminho mais longo válido (não -DBL_MAX)
            if ((VECTOR(bush->paths.dist_longest_local)[u] != -DBL_MAX) && (fluxo > 0)) {
                 if (VECTOR(bush->paths.dist_longest_local)[u] + weight > VECTOR(bush->paths.dist_longest_local)[v_node]) {
                    VECTOR(bush->paths.dist_longest_local)[v_node] = VECTOR(bush->paths.dist_longest_local)[u] + weight;
                    VECTOR(bush->paths.max_edges)[v_node] = id;
                    VECTOR(bush->paths.derivate_longest)[v_node] = dfluxo + VECTOR(bush->paths.derivate_longest)[u]; // Atualiza o gradiente do caminho máximo
                    
                }
            }
        }
        igraph_vector_int_clear(&incident_edges); // Reutiliza o vetor de arestas
    }
    //igraph_vector_int_print(&bush->paths.min_edges);
    //igraph_vector_int_print(&bush->paths.max_edges);

    //exit(0);
    

    // Libera memória dos vetores locais
    igraph_vector_int_destroy(&topological_order);
    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_destroy(&time);
    igraph_vector_destroy(&flow_bush);
    igraph_vector_destroy(&bush_PARAMETERS.capacidade);
    igraph_vector_destroy(&bush_PARAMETERS.cost_time);
    

}


bool ε_optimum(
    struct BUSH *bush, 
    struct BUSH *anterior_bush, 
    igraph_t *Grafo, 
    int i,
    int** edge_list,
    struct PARAMETERS* BPR_PARAMETERS,
    struct OD_MATRIX* OD,
    igraph_vector_t *flow
){
    igraph_vector_t time;
    igraph_vector_init(&time, BPR_PARAMETERS->L);
    
    int fonte = OD->Elementos[i].fonte;
    
    BPR(&time, BPR_PARAMETERS, flow);
    //print_flow(flow, Grafo, NULL); // Imprime o fluxo atual
    igraph_matrix_t res;
    igraph_matrix_init(&res, 0, BPR_PARAMETERS->N);
    igraph_distances_dijkstra(
        Grafo, &res,igraph_vss_1(fonte),igraph_vss_all(), &time, IGRAPH_OUT
    );
    bool converged = true;
    double tempo = 0;
    int edge_id,alvo,j,antecessor;

    for ( j = 0; j< bush->n_alvos; j++){

        alvo = VECTOR(OD->Elementos[i].alvos)[j];
        if(VECTOR(anterior_bush->paths.dist_longest_local)[alvo] - MATRIX(res,0,alvo) > EPSILON) {
            converged = false;
            //printf("\n\t\tCaminho não convergiu para Bush %d,%d: Maior tempo: %f, Tempo calculado: %f\n\n", fonte,alvo, VECTOR(anterior_bush->paths.dist_longest_local)[alvo], tempo);
            break;
        }
        
    }
    
    if(!converged) {
        //printf("\t\tBush %d não está em equilíbrio, recalculando...🔢🔢🔢🔢🔢🔢🔢🔢🔢\n", i+1);
        bush->is_ingraph = calloc(BPR_PARAMETERS->L,sizeof(bool));
        igraph_vector_init(&bush->flow_per_origin, BPR_PARAMETERS->L);
        igraph_vector_int_init(&bush->edge_id, 0);
        int add = 0;
        igraph_vector_int_t edges_vec;
        igraph_vector_int_init(&edges_vec, 0);
        int from, to;
        bool check;
        
        for (j = 0; j < BPR_PARAMETERS->L; j++) {
            from = edge_list[j][0];
            to = edge_list[j][1];
            // Check if edge should be in new topology
            //if (MATRIX(res,0,from) < MATRIX(res,0,to)){
            if(VECTOR(anterior_bush->paths.dist_shortest_local)[from] < VECTOR(anterior_bush->paths.dist_shortest_local)[to]){
                // Add edge to new bush
                igraph_vector_int_push_back(&edges_vec, j);
                bush->is_ingraph[j] = true;
                
                // Keep existing flow if edge was in previous bush
                //printf("A aresta (%d,%d) está na nova bush", from, to);
                VECTOR(bush->flow_per_origin)[j] = VECTOR(anterior_bush->flow_per_origin)[j]; // Mantém o fluxo existente
                if (!anterior_bush->is_ingraph[j]) {
                    //printf(", mas não estava na anterior %f\n", VECTOR(*flow)[j]);
                    add++;
                    VECTOR(bush->flow_per_origin)[j] = 0;
                }
                //else printf(" e estava na anterior %f\n", VECTOR(*flow)[j]);
            } 
            else {
                // Edge not in new bush, remove any existing flow
                bush->is_ingraph[j] = false;
                //printf("A aresta (%d,%d) não está na nova bush ", from, to);
                if (anterior_bush->is_ingraph[j]) {
                    ///printf("e estava na anterior, removendo fluxo %f\n", VECTOR(anterior_bush->flow_per_alvo)[j]);
                    VECTOR(*flow)[j] -= VECTOR(anterior_bush->flow_per_origin)[j];
                    if (VECTOR(*flow)[j] < 0) VECTOR(*flow)[j] = 0;
                    add++;
                }
                //else printf("e não estava na anterior, mantendo fluxo %f\n", VECTOR(*flow)[j]);
                VECTOR(bush->flow_per_origin)[j] = 0;

            }
        }
        //print_flow(flow, Grafo, NULL); // Imprime o fluxo atual
        //printf("Número de arestas adicionadas: %d\n", add);
        //exit(0); // Encerra o programa para depuração
        int index;
        igraph_es_t edges;
        igraph_es_vector(&edges,&edges_vec);
        igraph_subgraph_from_edges(Grafo, &bush->Grafo, edges, false);
        for (long e = 0; e < igraph_ecount(&bush->Grafo); e++) {
            igraph_integer_t from, to;
            igraph_edge(&bush->Grafo, e, &from, &to);
            
            index = find_id_log(from, to, edge_list, BPR_PARAMETERS->L);
            igraph_vector_int_push_back(&bush->edge_id, index);
            //printf("\n Final Edges = (%ld,%ld) %f\n", from, to, VECTOR(*flow)[index]);
        }

        //exit(0);
        // Copy everything from bush to anterior_bush
        anterior_bush->n_alvos = bush->n_alvos;
        for(j = 0; j < BPR_PARAMETERS->L; j++) {
            //anterior_bush->steps[j] = bush->steps[j];
            VECTOR(anterior_bush->flow_per_origin)[j] = VECTOR(bush->flow_per_origin)[j];
            anterior_bush->is_ingraph[j] = bush->is_ingraph[j];
        }
        igraph_vector_int_update(&anterior_bush->edge_id, &bush->edge_id);
        igraph_destroy(&anterior_bush->Grafo);
        igraph_copy(&anterior_bush->Grafo, &bush->Grafo);

        
        //print_flow(&bush->flow,&bush->grafo,NULL);
        igraph_vector_int_destroy(&edges_vec);
        
        
        // Copy new bush to anterior_bush
        //igraph_vector_int_print(&inbound);
        igraph_vector_destroy(&time);
        igraph_matrix_destroy(&res);
        //exit( 0); // Encerra o programa após atualizar a bush
        return false;

    }
    //igraph_vector_int_print(&inbound);
    igraph_vector_destroy(&time);
    igraph_matrix_destroy(&res);
    return true;
}

bool ε_equilibrium(
    struct BUSH *bush, 
    igraph_t *Grafo, 
    int i,
    struct PARAMETERS* BPR_PARAMETERS,
    struct OD_MATRIX* OD,
    igraph_vector_t *solucao
) {

    int alvo,j;
    int fonte = OD->Elementos[i].fonte;
    // ε-equilibrada
    double min_cost, max_cost;
    bool bush_is_converged;
    bool EQUILIBRED = true;

    while(true){

        find_dag_shortest_longest_costs_and_parents(
            bush,Grafo, fonte,BPR_PARAMETERS, solucao
        );
        bush_is_converged = true;
        
        //print_flow(solucao, Grafo, NULL); // Imprime o fluxo atualizado
        for ( j = 0; j < bush->n_alvos; j++) {

            alvo = VECTOR(OD->Elementos[i].alvos)[j];

            min_cost = VECTOR(bush->paths.dist_shortest_local)[alvo];
            max_cost = VECTOR(bush->paths.dist_longest_local)[alvo];

            if(max_cost == -DBL_MAX){

                VECTOR(bush->paths.dist_longest_local)[alvo] = 0;
                printf("\nCusto máximo -DBL_MAX encontrado para Bush %d, Alvo %d, ajustando para 0\n", fonte, alvo);
                exit(0); // Encerra o programa se o custo máximo for -DBL_MAX
                
                continue;
            }
            
            if(max_cost - min_cost > EPSILON){

                bush_is_converged = false; // A bush não está em equilíbrio
                EQUILIBRED = false;
                //printf("\n🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑 Caminho não equilibrado %d e %d 🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑🛑\n\n", fonte, alvo);
                int k,edge_id;
                /* print_edges(
                    fonte, alvo, &bush->paths.min_edges, Grafo
                );
                print_edges(
                    fonte, alvo, &bush->paths.max_edges, Grafo
                ); */
                igraph_vector_t time;
                igraph_vector_init(&time, BPR_PARAMETERS->L);
                igraph_vector_t dtime;
                igraph_vector_init(&dtime, BPR_PARAMETERS->L);
                BPR_derivate(&dtime,BPR_PARAMETERS, solucao); // Calcula o gradiente de cada aresta
                BPR(&time, BPR_PARAMETERS, solucao); // Calcula o tempo de cada aresta
                bool* visited_nodes = (bool*) calloc(BPR_PARAMETERS->N, sizeof(bool)); // Inicializa o vetor de nós visitados
                
                int node_M = alvo, node_N = -1,antecessor = alvo;
                int size_min = 0, size_max = 0;
                antecessor = alvo;
                while(antecessor != fonte) {
                    visited_nodes[antecessor] = true; // Marca o nó como visitado
                    antecessor = IGRAPH_FROM(Grafo, VECTOR(bush->paths.min_edges)[antecessor]);
                }
                visited_nodes[alvo] = true; // Marca a fonte como visitada
                visited_nodes[fonte] = true; // Marca a fonte como visitada
                double Numerator = 0,derivate = 0,step;
                double mu;
                while(node_N != fonte) {

                    Numerator = 0;
                    derivate = 0;
                    mu = 1e6;
                    node_N = get_princial_nodes( &node_M, &bush->paths, Grafo, fonte, visited_nodes);
                    if(node_M == node_N) break;
                    antecessor = node_M;
                    //printf("Nó M: %d, Nó N: %d\n", node_M, node_N);
                    //exit(0); // Encerra o programa para depuração
                    //exit(0); // Encerra o programa para depuração
                    while(antecessor!=fonte){
                        edge_id = VECTOR(bush->paths.min_edges)[antecessor];
                        //if(edge_id == 45) printf("Aresta %d (%ld,%ld) - %f\n", edge_id, IGRAPH_FROM(Grafo, edge_id), IGRAPH_TO(Grafo, edge_id), VECTOR(bush->flow_per_origin)[edge_id]);
                        Numerator -= VECTOR(time)[edge_id];
                        derivate += VECTOR(dtime)[edge_id];
                        antecessor = IGRAPH_FROM(Grafo, edge_id);
                        //printf("%d\n", antecessor);
                        if(antecessor == node_N) break; // Encerra o programa se antecessor for fonte
                    }
                    antecessor = node_M;
                    while(antecessor!=fonte){
                        edge_id = VECTOR(bush->paths.max_edges)[antecessor];
                        //if(edge_id == 45) printf("Aresta %d (%ld,%ld) - %f\n", edge_id, IGRAPH_FROM(Grafo, edge_id), IGRAPH_TO(Grafo, edge_id), VECTOR(bush->flow_per_origin)[edge_id]);
                        Numerator += VECTOR(time)[edge_id];
                        derivate += VECTOR(dtime)[edge_id];
                        mu = fmin(mu, (double)VECTOR(bush->flow_per_origin)[edge_id]); // Atualiza o valor de mu
                        antecessor = IGRAPH_FROM(Grafo, edge_id);
                        //printf("%d\n", antecessor);
                        if(antecessor == node_N) break; // Encerra o programa se antecessor for fonte
                    }

                    //printf("Numerador: %f, Derivada: %f, Delta: %f, Node_N: %d, Node_M:%d\n", Numerator, derivate,max_cost - min_cost,node_N, node_M);
                    step = Numerator / derivate; // Calcula o passo

                    step = fmin(step, mu); // Limita o passo ao valor de mu
                    if(derivate == 0.0) exit(0); // Encerra o programa se a derivada for zero
                    if(step < 0) step = 0.0; // Ajusta o passo para 0 se for negativo
                    antecessor = node_M;
                    while(antecessor != fonte) {
                        edge_id = VECTOR(bush->paths.min_edges)[antecessor];
                        VECTOR(*solucao)[edge_id] += step;
                        VECTOR(bush->flow_per_origin)[edge_id] += step;
                        antecessor = IGRAPH_FROM(Grafo, edge_id);
                        if(antecessor == node_N) break;
                    }
                    antecessor = node_M;
                    while(antecessor != fonte) {
                        edge_id = VECTOR(bush->paths.max_edges)[antecessor];
                        VECTOR(*solucao)[edge_id] -= step;
                        VECTOR(bush->flow_per_origin)[edge_id] -= step; 
                        if(VECTOR(*solucao)[edge_id] < 0){
                            printf("Fluxo negativo encontrado na aresta %d (%ld,%ld), ajustando para %e\n", edge_id, IGRAPH_FROM(Grafo, edge_id), IGRAPH_TO(Grafo, edge_id), VECTOR(*solucao)[edge_id]);
                            exit(0); // Encerra o programa se o fluxo for negativo
                        }
                        if(VECTOR(bush->flow_per_origin)[edge_id] < 0){
                            printf("Fluxo negativo encontrado na aresta %d (%ld,%ld), ajustando para 0\n", edge_id, IGRAPH_FROM(Grafo, edge_id), IGRAPH_TO(Grafo, edge_id));
                            exit(0); // Encerra o programa se o fluxo for negativo
                        }
                        if(VECTOR(bush->flow_per_origin)[edge_id] > VECTOR(*solucao)[edge_id]){
                            printf("Fluxo da bush maior que o fluxo total na aresta %d (%ld,%ld), ajustando para %f\n", edge_id, IGRAPH_FROM(Grafo, edge_id), IGRAPH_TO(Grafo, edge_id), VECTOR(*solucao)[edge_id]);
                            VECTOR(bush->flow_per_origin)[edge_id] = VECTOR(*solucao)[edge_id]; // Ajusta o fluxo da bush para o fluxo total
                            exit(0); // Encerra o programa se o fluxo da bush for maior que o fluxo total
                        }
                        antecessor = IGRAPH_FROM(Grafo, edge_id);
                        if(antecessor == node_N) break;
                    }
                    node_M = node_N; // Atualiza o nó M para o nó N
                }
                //exit(0); // Encerra o programa após equilibrar a bush
                igraph_vector_destroy(&time);
                igraph_vector_destroy(&dtime);
                free(visited_nodes); // Libera a memória do vetor de nós visitados
                //igraph_vector_print(solucao);
                /* find_dag_shortest_longest_costs_and_parents(
                    bush,Grafo, fonte,BPR_PARAMETERS, solucao
                ); */
                break;
            }
        }
        if(bush_is_converged) break; // Se todas as bushes estão em equilíbrio, sai do loop
    }
    return EQUILIBRED; // Retorna verdadeiro se a bush estiver em equilíbrio
    
}

void Bushes(
    igraph_t *Grafo,
    struct OD_MATRIX* OD,
    int** edge_list,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *solucao
) {

    igraph_vector_init(solucao, BPR_PARAMETERS->L); // Inicializa a solução com zeros
    int i,count,iter = 0;
    struct BUSH* bushes = (struct BUSH*) malloc(OD->size * sizeof(struct BUSH));
    
    for ( i = 0; i < OD->size; i++){

        bushes[i].n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
        igraph_vector_init(&bushes[i].flow_per_origin, BPR_PARAMETERS->L); // Inicializa o vetor de fluxo para cada bush
        init_bush(&bushes[i], Grafo, i,edge_list, BPR_PARAMETERS,OD,solucao);
        
    }

    bool IS_EPSLION_GLOBAL_CONVERGED = false;
    bool IS_EPSLION_EQUILIBRED = false;
    bool IS_EPSLION_OPTIMAL = false;
    bool ALL_BUSHES_CONVERGED = false;
    FILE *file = fopen("gap.txt", "w");

    while(!ALL_BUSHES_CONVERGED) {

        ALL_BUSHES_CONVERGED = true; // Inicializa como verdadeiro, será atualizado se alguma bush não convergir
        //printf("🔢🔢🔢🔢🔢🔢🔢🔢 Iteração %d 🔢🔢🔢🔢🔢🔢🔢🔢\n", iter);

        for ( i = 0; i < OD->size; i++){

            IS_EPSLION_GLOBAL_CONVERGED = false; // Reseta a convergência global para cada bush
            IS_EPSLION_EQUILIBRED = false; // Reseta o equilíbrio
            IS_EPSLION_OPTIMAL = false; // Reseta a otimalidade
            count = 0; // Reseta o contador de iterações para cada bush

            while (!IS_EPSLION_GLOBAL_CONVERGED){
                count++;

                IS_EPSLION_EQUILIBRED = ε_equilibrium(
                    &bushes[i], 
                    Grafo, 
                    i, 
                    BPR_PARAMETERS, 
                    OD, 
                    solucao
                );

                //printf("\n\t\t👌👌👌👌👌👌👌👌👌👌👌👌 Bush %d está em equilíbrio 👌👌👌👌👌👌👌👌👌👌👌👌\n\n", i+1);
                //print_flow(solucao, Grafo, NULL); // Imprime o fluxo final 
                //print_bush(&bushes[i], solucao, BPR_PARAMETERS,NULL ); // Imprime a bush atual
                //printf("Bush %d convergiu com %d alvos\n", i, n_alvos);
                check_solution(bushes,solucao,BPR_PARAMETERS,OD->size);
                //ε-ótima
                
                struct BUSH new_bushes;
                new_bushes.n_alvos = bushes[i].n_alvos;
                IS_EPSLION_OPTIMAL = ε_optimum( &new_bushes,&bushes[i], Grafo, i, edge_list, BPR_PARAMETERS, OD, solucao);
                //printf("\n\n\t\t🎉🎉🎉🎉🎉🎉 Bush %d é ótima! 🎉🎉🎉🎉🎉🎉\n", i+1);
                IS_EPSLION_GLOBAL_CONVERGED = IS_EPSLION_EQUILIBRED && IS_EPSLION_OPTIMAL;
            
                //print_bush(&bushes[i], solucao, BPR_PARAMETERS,NULL ); // Imprime a bush atual
                //igraph_vector_print(&bushes[i].paths.dist_shortest_local);

                //print_flow(solucao, Grafo,&bushes[i], NULL); // Imprime o fluxo final 
                //exit(0); // Encerra o programa após equilibrar a bush
                check_solution(bushes,solucao,BPR_PARAMETERS,OD->size);
                
            }
            //exit(0); // Encerra o programa após equilibrar a bush
            ALL_BUSHES_CONVERGED = ALL_BUSHES_CONVERGED && (count == 1); // Se alguma bush não convergiu, ALL_BUSHES_CONVERGED será falso
            //printf("\n\n\t\t🎉🎉🎉🎉🎉🎉 Bush %d é ótima! 🎉🎉🎉🎉🎉🎉\n", i+1);
            //igraph_vector_print(solucao);
            //exit(0); // Encerra o programa após equilibrar a bush
        }
        double GAP = relative_gap(solucao,Grafo, BPR_PARAMETERS,OD);
        fprintf(file,"%d %e\n", iter, GAP);
        iter++;
    }
    fclose(file);
    print_flow(solucao, Grafo, NULL); // Imprime o fluxo final de todas as bushes
    int j;
    for(i =0; i < OD->size; i++){
        int size = igraph_vector_int_size(&OD->Elementos[i].alvos);
        for (j = 0; j < size; j++) {
            printf("Bush %d, Alvo %ld: ",  OD->Elementos[i].fonte, VECTOR(OD->Elementos[i].alvos)[j]);
            calc_probability(
                &bushes[i], 
                VECTOR(OD->Elementos[i].alvos)[j], 
                BPR_PARAMETERS, 
                solucao,
                solucao
            );
        }
        exit(0); // Encerra o programa após equilibrar todas as bushes
    }
}