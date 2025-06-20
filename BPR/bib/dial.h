
#pragma once

#include "igraph.h"
#include "network.h"

void find_dag_shortest_longest_costs_and_parents(
    struct BUSH *bush,
    igraph_t *Grafo,
    igraph_integer_t from,
    igraph_integer_t to,
    struct min_max_bush *result,
    struct PARAMETERS *BPR_PARAMETERS,
    igraph_vector_t *flow
) {

    //printf("Encontrando caminhos entre %ld ➡️  %ld\n", from, to);
    int N = igraph_vcount(Grafo);
    // Garante que os vetores de saída estejam dimensionados e inicializados
    // Vetores locais para armazenar distâncias para todos os nós
    init_path(result, N);
    result->mu = 1e6;
    VECTOR(result->dist_shortest_local)[from] = 0.0;
    VECTOR(result->dist_longest_local)[from] = 0.0;

    result->min_cost = 0.0;
    result->max_cost = 0.0;

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
    double fluxo;
    for(i = 0; i < L; i++) {
        index = VECTOR(bush->edge_id)[i];
        igraph_vector_push_back(&flow_bush, VECTOR(*flow)[index]);
        igraph_vector_push_back(&bush_PARAMETERS.capacidade, VECTOR(BPR_PARAMETERS->capacidade)[index]);
        igraph_vector_push_back(&bush_PARAMETERS.cost_time, VECTOR(BPR_PARAMETERS->cost_time)[index]);

    }
    BPR(&time, &bush_PARAMETERS, &flow_bush); // Calcula o tempo de cada aresta
    double _ = 0;
    BPR_derivate(
        &dtime, 
        &bush_PARAMETERS,
        flow,
        &_ // Passa um ponteiro para o total_time, mas não é usado aqui
    );
    //printf("Custo do tempo: \n");
    //igraph_vector_print(&time); // Imprime o tempo de cada aresta
    igraph_vector_int_t min_edges, max_edges;
    igraph_vector_int_init(&min_edges, N);
    igraph_vector_int_init(&max_edges, N);
    igraph_vector_int_fill(&min_edges, -1); // -1 indica nenhum predecessor
    igraph_vector_int_fill(&max_edges, -1); // -1 indica nenhum predecessor

    //printf("Fluxo do tempo: \n");
    //print_flow(&flow_bush, &bush->Grafo, NULL); // Imprime o fluxo do tempo
    for (i = 0; i < igraph_vector_int_size(&topological_order); ++i) {
        igraph_integer_t u = VECTOR(topological_order)[i];
        // Se 'u' não é alcançável a partir de 'from' (para menor caminho),
        // não pode contribuir para caminhos subsequentes.
        if (VECTOR(result->dist_shortest_local)[u] == DBL_MAX) continue;
        
        igraph_incident(&bush->Grafo, &incident_edges, u, IGRAPH_OUT);
        
        for (long j = 0; j < igraph_vector_int_size(&incident_edges); ++j) {

            igraph_integer_t eid = VECTOR(incident_edges)[j];
            id = VECTOR(bush->edge_id)[eid]; // Acessa o id da aresta
            igraph_integer_t v_node = IGRAPH_TO(&bush->Grafo, eid);
            igraph_real_t weight,flow,derivate_weight;

            weight = VECTOR(time)[eid]; // Acessa o peso da aresta
            fluxo = bush->flow_per_alvo[id]; // Acessa o fluxo da aresta
            derivate_weight = VECTOR(dtime)[eid]; // Acessa o derivado do tempo da aresta
            if (VECTOR(result->dist_shortest_local)[u] + weight < VECTOR(result->dist_shortest_local)[v_node]) {
                
                VECTOR(result->dist_shortest_local)[v_node] = VECTOR(result->dist_shortest_local)[u] + weight;
                // Atualiza o predecessor para o menor caminho
                VECTOR(min_edges)[v_node] = id;
                if(v_node == to) result->min_cost = VECTOR(result->dist_shortest_local)[v_node];
                // Atualiza o derivado do menor caminho
            }

            // Relaxamento para o maior caminho
            // Só estender de 'u' se 'u' já tem um caminho mais longo válido (não -DBL_MAX)
            if ((VECTOR(result->dist_longest_local)[u] != -DBL_MAX) && (fluxo > 0)) {
                 if (VECTOR(result->dist_longest_local)[u] + weight > VECTOR(result->dist_longest_local)[v_node]) {
                    //printf("Entrou!!!!!!!! %ld %ld\n",u,v_node);
                    VECTOR(result->dist_longest_local)[v_node] = VECTOR(result->dist_longest_local)[u] + weight;
                    VECTOR(max_edges)[v_node] = id;
                    if(v_node == to) result->max_cost = VECTOR(result->dist_longest_local)[v_node];
                }
            }
        }
        igraph_vector_int_clear(&incident_edges); // Reutiliza o vetor de arestas
    }

    //igraph_vector_int_print(&min_edges);
    //igraph_vector_int_print(&max_edges);
    //printf("Caminho mínimo encontrado: ");
    //igraph_vector_int_print(&result->min_edges);
    if(result->max_cost != 0){
        int alvo = to;
        while(VECTOR(min_edges)[alvo] != -1) {
            id = VECTOR(min_edges)[alvo];
            igraph_vector_int_push_back(&result->min_edges, id); // Adiciona a aresta ao caminho mínimo
            //printf("(%ld,%ld, %f)",IGRAPH_FROM(Grafo, id), IGRAPH_TO(Grafo, id),VECTOR(dtime)[id]);
            alvo = IGRAPH_FROM(Grafo, id);
        }
        //printf("\n");
        alvo = to;
        while(VECTOR(max_edges)[alvo] != -1) {
            id = VECTOR(max_edges)[alvo];
            igraph_vector_int_push_back(&result->max_edges, id); // Adiciona a aresta ao caminho mínimo
            fluxo = bush->flow_per_alvo[id]; // Acessa o fluxo da aresta
            //printf("(%ld,%ld, %f)",IGRAPH_FROM(Grafo, id), IGRAPH_TO(Grafo, id),VECTOR(dtime)[id]);
            //printf("%d,%f\n",alvo,fluxo);
            if(fluxo < result->mu) result->mu = fluxo; // Atualiza o valor de mu se o fluxo for menor
            alvo = IGRAPH_FROM(Grafo, id);
        }
        //printf("\n");
        int node_M,node_N,node1,node2,index_m;
        node1 = IGRAPH_TO(Grafo, VECTOR(result->min_edges)[0]);
        node2 = IGRAPH_TO(Grafo, VECTOR(result->max_edges)[0]);
        i = 0;
        while(node1 == node2) {
            node_M = node1;
            index_m = i;
            if(node1 == from) break;
            node1 = IGRAPH_FROM(Grafo, VECTOR(result->min_edges)[i]);
            node2 = IGRAPH_FROM(Grafo, VECTOR(result->max_edges)[i]);
            i++;
        }
        int l_min = igraph_vector_int_size(&result->min_edges);
        int l_max = igraph_vector_int_size(&result->max_edges);
        node1 = IGRAPH_FROM(Grafo, VECTOR(result->min_edges)[l_min-1]);
        node2 = IGRAPH_FROM(Grafo, VECTOR(result->max_edges)[l_max-1]);
        i = 0;
        //printf("%d %d\n", node1, node2);
        while(node1 == node2) {
            node_N = node1;
            if(node1 == to) break;
            node1 = IGRAPH_TO(Grafo, VECTOR(result->min_edges)[l_min-1 - i]);
            node2 = IGRAPH_TO(Grafo, VECTOR(result->max_edges)[l_max-1 - i]);
            i++;
        }
        if((node_M != from) || (node_N != to)) {
            //printf("m = %d,n = %d, i = %d\n", node_M, node_N,index_m);
            //printf("Caminho mínimo e máximo não são iguais, Bush não é ótima!\n");
            result->derivate = 0.0; // Zera o derivado se os caminhos não forem iguais
            for(i = index_m; i < igraph_vector_int_size(&result->min_edges); i++) {
                index = VECTOR(result->min_edges)[i];
                //printf("Removendo aresta %f (%ld,%ld)\n", VECTOR(dtime)[index], IGRAPH_FROM(Grafo, index), IGRAPH_TO(Grafo, index));
                result->derivate += VECTOR(dtime)[index];
                if(IGRAPH_FROM(Grafo, index) == from) break;
            }
            for(i = index_m; i < igraph_vector_int_size(&result->max_edges); i++) {
                index = VECTOR(result->max_edges)[i];
                //printf("Removendo aresta %f (%ld,%ld)\n", VECTOR(dtime)[index], IGRAPH_FROM(Grafo, index), IGRAPH_TO(Grafo, index));
                result->derivate += VECTOR(dtime)[index];
                if(IGRAPH_FROM(Grafo, index) == from) break;
            }
            //printf("Derivada: %f\n", result->derivate);
            //exit(0);
        }

    }
    //exit(0);
    //igraph_vector_int_print(&result->max_parents);
    //igraph_vector_int_print(&inbound_edges);
    

    // Libera memória dos vetores locais
    igraph_vector_int_destroy(&topological_order);
    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_destroy(&time);
    igraph_vector_destroy(&dtime);
    igraph_vector_destroy(&flow_bush);
    igraph_vector_destroy(&bush_PARAMETERS.capacidade);
    igraph_vector_destroy(&bush_PARAMETERS.cost_time);
    igraph_vector_int_destroy(&min_edges);
    igraph_vector_int_destroy(&max_edges);
    

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
    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL,NULL,fonte,igraph_vss_all(),&BPR_PARAMETERS->cost_time, IGRAPH_OUT,NULL,&inbound
    );
    igraph_vector_int_init(&bush->edge_id, 0);
    bush->is_ingraph = (bool*) calloc(BPR_PARAMETERS->L, sizeof(bool)); // Inicializa o vetor de arestas no grafo
    int index,target,antecessor;
    igraph_vector_int_t edges_vec;
    igraph_vector_int_init(&edges_vec, 0);
    for (k = 0; k < BPR_PARAMETERS->N; k++) {
        if(k == fonte) continue; // Pula a fonte
        index = VECTOR(inbound)[k];
        igraph_vector_int_push_back(&edges_vec, index);
        bush->is_ingraph[index] = true; // Marca a aresta como parte do grafo da bush
    }
    //printf("Size of bush->flow: %ld\n", igraph_vector_size(&bush->flow));
    //igraph_vector_int_print(&edges_vec);
    igraph_es_t edges;
    igraph_es_vector(&edges,&edges_vec);
    
    igraph_subgraph_from_edges(
        Grafo, &bush->Grafo, edges, false
    );
    //printf("Number of edges in Bush: %ld\n", igraph_ecount(&bush->grafo));
    for(j = 0; j < bush->n_alvos; j++) {
        antecessor = VECTOR(OD->Elementos[i].alvos)[j];
        bush->steps[j] = 1.; // Inicializa o passo para cada alvo
        while(antecessor != fonte) {
            index = VECTOR(inbound)[antecessor];
            VECTOR(*flow)[index] += VECTOR(OD->Elementos[i].volumes)[j];
            //printf("Adicionando fluxo %ld na aresta %d (%ld,%ld)\n", VECTOR(OD->Elementos[i].volumes)[j],index, IGRAPH_FROM(Grafo, index), IGRAPH_TO(Grafo, index));
            bush->flow_per_alvo[index] += VECTOR(OD->Elementos[i].volumes)[j];
            antecessor = IGRAPH_FROM(Grafo, index);

        }
    }
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
    igraph_vector_int_t inbound;
    igraph_vector_t time;
    igraph_vector_init(&time, BPR_PARAMETERS->L);
    igraph_vector_int_init(&inbound, 0);
    
    int fonte = OD->Elementos[i].fonte;
    
    BPR(&time, BPR_PARAMETERS, flow);
    //print_flow(flow, Grafo, NULL); // Imprime o fluxo atual
    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL,NULL,fonte,igraph_vss_all(),&time, IGRAPH_OUT,NULL,&inbound
    );

    
    bool converged = true;
    double tempo = 0;
    int edge_id,alvo,j;

    for ( j = 0; j< bush->n_alvos; j++){

        alvo = VECTOR(OD->Elementos[i].alvos)[j];
        tempo = 0;

        while(alvo!=fonte){
            edge_id = VECTOR(inbound)[alvo];
            tempo += VECTOR(time)[edge_id];
            //printf("Edge %d: (%ld,%ld) - %f\n", alvo, IGRAPH_FROM(Grafo, edge_id), IGRAPH_TO(Grafo, edge_id), tempo);
            alvo = IGRAPH_FROM(Grafo, edge_id);
        }
        if(fabs(tempo - anterior_bush->paths->max_cost) > EPSILON) {
            converged = false;
            //printf("\n\t\tCaminho não convergiu para Bush %d,%d: Maior tempo: %f, Tempo calculado: %f\n\n", fonte,(int)VECTOR(OD->Elementos[i].alvos)[j], anterior_bush->paths->max_cost, tempo);
            break;
        }
        
    }
    
    if(!converged) {
        //printf("\t\tBush %d não está em equilíbrio, recalculando...🔢🔢🔢🔢🔢🔢🔢🔢🔢\n", i+1);

        double* total_time = (double*) calloc(BPR_PARAMETERS->N , sizeof(double));
        bool* keep_flow = (bool*) calloc(BPR_PARAMETERS->L , sizeof(bool));
        bool* is_ingraph = calloc(BPR_PARAMETERS->L,sizeof(bool));
        bush->is_ingraph = calloc(BPR_PARAMETERS->L,sizeof(bool));
        bush->flow_per_alvo = (double*) calloc(BPR_PARAMETERS->L , sizeof(double));
        igraph_vector_int_init(&bush->edge_id, 0);

        igraph_vector_int_t edges_vec;
        igraph_vector_int_init(&edges_vec, 0);
        int from, to;
        bool check;
        for ( j = 0; j< BPR_PARAMETERS->N; j++){
            if(j == fonte) continue; // Pula a fonte
            alvo = j;
            tempo = 0;
            while(alvo!=fonte){
                edge_id = VECTOR(inbound)[alvo];
                tempo += VECTOR(time)[edge_id];
                
                alvo = IGRAPH_FROM(Grafo, edge_id);
                bush->is_ingraph[edge_id] = true; // Marca a aresta como não parte do grafo da bush
                
            }
            total_time[j] = tempo;
        }
        
        for ( j = 0; j< BPR_PARAMETERS->L; j++){
            from = edge_list[j][0];
            to = edge_list[j][1];
            check = (total_time[from] < total_time[to]);
            
            if(check) {
                
                if(anterior_bush->is_ingraph[j]) {
                    igraph_vector_int_push_back(&edges_vec,j);
                    is_ingraph[j] = true; // Marca a aresta como parte do grafo da bush
                    keep_flow[j] = true; // Marca a aresta para manter o fluxo
                }
                else{
                    igraph_vector_int_push_back(&edges_vec,j);
                    is_ingraph[j] = true; // Marca a aresta como parte do grafo da bush
                    keep_flow[j] = false; 
                    //printf("%d %d %f %d\n",from,to,VECTOR(*flow)[j],keep_flow[j]);
                }
            }
            bush->is_ingraph[j] = false;
            bush->flow_per_alvo[j] = anterior_bush->flow_per_alvo[j]; // Zera o fluxo por alvo na bush
        }
        //print_flow(flow, Grafo, NULL); // Imprime o fluxo atual
        int index;
        igraph_es_t edges;
        igraph_es_vector(&edges,&edges_vec);
        igraph_subgraph_from_edges(Grafo, &bush->Grafo, edges, false);
        for (long e = 0; e < igraph_ecount(&bush->Grafo); e++) {
            igraph_integer_t from, to;
            igraph_edge(&bush->Grafo, e, &from, &to);
            
            index = find_id_log(from, to, edge_list, BPR_PARAMETERS->L);
            igraph_vector_int_push_back(&bush->edge_id, index);
            is_ingraph[index] = true; // Marca a aresta como parte do grafo da bush
            if(!keep_flow[index]) bush->flow_per_alvo[index] = 0; // Zera o fluxo por alvo na bush

            //printf("\n Final Edges = (%ld,%ld) %f\n", from, to, VECTOR(*flow)[index]);
        }

        
        //exit(0);
        free(total_time);
        free(is_ingraph);
        free(keep_flow);
        igraph_vector_int_update(&anterior_bush->edge_id, &bush->edge_id);
        memcpy(bush->is_ingraph, is_ingraph, BPR_PARAMETERS->L * sizeof(bool));
        igraph_copy(&anterior_bush->Grafo,&bush->Grafo);
        anterior_bush->n_alvos = bush->n_alvos; // Atualiza o número de alvos da bush anterior
        memcpy(anterior_bush->flow_per_alvo, bush->flow_per_alvo, BPR_PARAMETERS->L * sizeof(double));

        
        //print_flow(&bush->flow,&bush->grafo,NULL);
        igraph_vector_int_destroy(&edges_vec);


        // Copy new bush to anterior_bush
        //igraph_vector_int_print(&inbound);
        igraph_vector_int_destroy(&inbound);
        igraph_vector_destroy(&time);
        return false;

    }
    printf("Bush %d é ótima! ✅\n", i+1);
    igraph_vector_int_print(&inbound);
    igraph_vector_int_destroy(&inbound);
    igraph_vector_destroy(&time);
    return true;
}

void Bushes(
    igraph_t *Grafo,
    struct OD_MATRIX* OD,
    int** edge_list,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *solucao
) {

    int iter = 0;
    igraph_vector_init(solucao, BPR_PARAMETERS->L); // Inicializa a solução com zeros
    int i,c = 0,alvo,fonte,j,converged_bush;
    struct BUSH* bushes = (struct BUSH*) malloc(OD->size * sizeof(struct BUSH));
    
    for ( i = 0; i < OD->size; i++){

        
        bushes[i].n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
        bushes[i].paths = (struct min_max_bush*) malloc(bushes[i].n_alvos * sizeof(struct min_max_bush));
        bushes[i].flow_per_alvo = (double*) calloc(BPR_PARAMETERS->L , sizeof(double));
        
        //printf("Iniciando Bush %d\n",i+1);
        bushes[i].steps = (double*) calloc(bushes[i].n_alvos , sizeof(double));
        init_bush(&bushes[i], Grafo, i,edge_list, BPR_PARAMETERS,OD,solucao);
        printf("\n\n🎉🎉🎉🎉🎉🎉 Bush Iniciada! %d 🎉🎉🎉🎉🎉🎉\n", i+1);
    }
    bool IS_EPSLION_GLOBAL_CONVERGED = false;

    while (!IS_EPSLION_GLOBAL_CONVERGED){
        
        iter++;
        if(iter > 50) break; // Encerra o programa se não convergir após 3 iterações
        IS_EPSLION_GLOBAL_CONVERGED = true; // Assume que todas as bushes estão convergidas
        printf("🔢🔢🔢🔢🔢🔢🔢🔢 Iteração %d 🔢🔢🔢🔢🔢🔢🔢🔢\n", iter);
        // Verifica se todas as bushes convergiram

        for ( i = 0; i < OD->size; i++){

            converged_bush = 0;
            fonte = OD->Elementos[i].fonte;
            bool* bush_is_converged = (bool*) calloc(bushes[i].n_alvos , sizeof(bool));
            // ε-equilibrada
            
            while(converged_bush < bushes[i].n_alvos){
                converged_bush = 0; // Reseta o contador de bushes convergidas
                for ( j = 0; j < bushes[i].n_alvos; j++) {
                    //printf("\n Equilibrando Bush %d, Alvo %d:\n", i+1, j+1);
                    alvo = VECTOR(OD->Elementos[i].alvos)[j];
                    if(bush_is_converged[j]){
                        converged_bush++;
                        continue; // Se já está convergida, pula para o próximo alvo
                    }
                    find_dag_shortest_longest_costs_and_parents(
                        &bushes[i],Grafo, fonte, alvo, &bushes[i].paths[j],BPR_PARAMETERS, solucao
                    );
                    //printf("\nCusto mínimo: %f(%f), Custo máximo: %f(%f)\n", bushes[i].paths[j].min_cost,bushes[i].paths[j].derivate_min_cost, bushes[i].paths[j].max_cost,bushes[i].paths[j].derivate_max_cost);
                    bush_is_converged[j] = (bushes[i].paths[j].max_cost - bushes[i].paths[j].min_cost < EPSILON);
                    //printf("%f\n",bushes[i].paths[j].max_cost);
                    if(!bush_is_converged[j]){
                        //printf("\nCaminho com custo diferente encontrado entre %d e %d\n", fonte, alvo);
                        //printf("\nCusto mínimo: %f, Custo máximo: %f Derivada: %f\n", bushes[i].paths[j].min_cost, bushes[i].paths[j].max_cost, bushes[i].paths[j].derivate);
                        //printf("Pais mínimos: ");
                        // Remove first and last elements from min_parents vector
                        int l_min = igraph_vector_int_size(&bushes[i].paths[j].min_edges);
                        int l_max = igraph_vector_int_size(&bushes[i].paths[j].max_edges);
                        //igraph_vector_int_print(&bushes[i].paths[j].min_parents);
                        //printf("Pais máximos: ");
                        // Remove first and last elements from max_parents vector
                        //igraph_vector_int_print(&bushes[i].paths[j].max_parents);
                        
                        
                        int edge_id,l,k = 0;
                        
                        double Numerator = bushes[i].paths[j].max_cost - bushes[i].paths[j].min_cost;
                        double step = bushes[i].steps[j] + Numerator / bushes[i].paths[j].derivate; // Calcula o passo
                        step = fmin(step, bushes[i].paths[j].mu); // Limita o passo ao valor de mu
                        //exit(0); // Encerra o programa para depuração
                        if(bushes[i].paths[j].derivate == 0.0) {
                            //print_flow(solucao, Grafo, NULL); // Imprime o fluxo atual
                            erase_path(&bushes[i].paths[j]); // Libera a memória da estrutura de caminho
                            continue;
                        }
                        if(step < 0) {
                            printf("Passo negativo encontrado: %f, ajustando para 0\n", step);
                            exit(0); // Encerra o programa se o passo for negativo
                            step = 0.0; // Ajusta o passo para 0 se for negativo
                        }
                        int to = alvo;
                        for(k = 0; k < l_min; k++){
                            edge_id = VECTOR(bushes[i].paths[j].min_edges)[k];
                            //printf("Atualizando aresta (%d,%d) com passo %f\n", node1, node2, step);
                            VECTOR(*solucao)[edge_id] += step;
                            bushes[i].flow_per_alvo[edge_id] += step; // Atualiza o fluxo por alvo
                        }
                        for(k = 0; k < l_max; k++){
                            edge_id = VECTOR(bushes[i].paths[j].max_edges)[k];
                            VECTOR(*solucao)[edge_id] -= step;
                            bushes[i].flow_per_alvo[edge_id] -= step; // Atualiza o fluxo por alvo
                        }
                        
                        //igraph_vector_print(&bushes[i].flow);
                        bushes[i].steps[j] = step;
                        //printf("\nDenominator: %e, Numerator: %e, step: %e %f\n", Denominator, Numerator, bushes[i].steps[j],bushes[i].paths[j].mu);
                        //print_flow(solucao, Grafo, NULL); // Imprime o fluxo atualizado
                        //exit(0); // Encerra o programa após encontrar o passo

                    }
                    else converged_bush += bush_is_converged[j];
                    erase_path(&bushes[i].paths[j]); // Libera a memória da estrutura de caminho
                }
            }
            free(bush_is_converged);
            //printf("\n\t\t👌👌👌👌👌👌👌👌👌👌👌👌 Bush %d está em equilíbrio 👌👌👌👌👌👌👌👌👌👌👌👌\n\n", i+1);
            //print_flow(solucao, Grafo, NULL); // Imprime o fluxo final 
            
            //exit(0); // Se não encontrar caminho, encerra o programa
            //printf("Bush %d convergiu com %d alvos\n", i, n_alvos);

            //ε-ótima
            //igraph_vector_print(&bushes[i].flow);
            
            struct BUSH new_bushes;
            new_bushes.n_alvos = bushes[i].n_alvos;
            new_bushes.paths = (struct min_max_bush*) malloc(new_bushes.n_alvos * sizeof(struct min_max_bush));

            //if(i == 0) print_flow(solucao, Grafo);
            //if(i == 0) print_flow(&bushes[i].flow, &bushes[i].grafo);

            bool s = ε_optimum( &new_bushes,&bushes[i], Grafo, i, edge_list, BPR_PARAMETERS, OD, solucao);
            IS_EPSLION_GLOBAL_CONVERGED = IS_EPSLION_GLOBAL_CONVERGED && s;
            //exit(0); // Encerra o programa se não convergir após 3 iterações
            //printf("\n\n\t\t🎉🎉🎉🎉🎉🎉 Bush %d é ótima! 🎉🎉🎉🎉🎉🎉\n", i+1);
            //print_flow(solucao, Grafo, NULL); // Imprime o fluxo final 
            //exit(0); // Encerra o programa após equilibrar a bush
            
        }
    }
    igraph_vector_print(solucao);
}