
#pragma once

#include "igraph.h"
#include "network.h"

void find_dag_shortest_longest_costs_and_parents(
    struct BUSH *bush,
    igraph_integer_t from,
    igraph_integer_t to,
    struct min_max_bush *result
) {

    printf("Encontrando caminhos entre %ld ➡️  %ld\n", from, to);
    igraph_integer_t N = igraph_vcount(&bush->grafo);

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
    if (igraph_topological_sorting(&bush->grafo, &topological_order, IGRAPH_OUT) != IGRAPH_SUCCESS) {
        printf("Erro ao ordenar topologicamente o grafo.\n");
        exit(0); // Falha ao ordenar topologicamente, encerra o programa
    }

    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);
    igraph_vector_t time;
    igraph_vector_init(&time, bush->BPR_PARAMETERS.L);
    igraph_vector_t dtime;
    igraph_vector_init(&dtime, bush->BPR_PARAMETERS.L);

    BPR(&time, &bush->BPR_PARAMETERS, &bush->flow); // Calcula o tempo de cada aresta
    double _ = 0;
    BPR_derivate(
        &dtime, 
        &bush->BPR_PARAMETERS,
        &bush->flow,
        &_ // Passa um ponteiro para o total_time, mas não é usado aqui
    );
    printf("Fluxo do tempo: ");
    igraph_vector_print(&bush->flow);
    for (int i = 0; i < igraph_vector_int_size(&topological_order); ++i) {
        igraph_integer_t u = VECTOR(topological_order)[i];
        // Se 'u' não é alcançável a partir de 'from' (para menor caminho),
        // não pode contribuir para caminhos subsequentes.
        if (VECTOR(result->dist_shortest_local)[u] == DBL_MAX) continue;

        igraph_incident(&bush->grafo, &incident_edges, u, IGRAPH_OUT);
        for (long j = 0; j < igraph_vector_int_size(&incident_edges); ++j) {
            igraph_integer_t eid = VECTOR(incident_edges)[j];
            igraph_integer_t v_node = IGRAPH_TO(&bush->grafo, eid);
            igraph_real_t weight,flow,derivate_weight;

            weight = VECTOR(time)[eid]; // Acessa o peso da aresta
            flow = VECTOR(bush->flow)[eid]; // Acessa o fluxo da aresta
            derivate_weight = VECTOR(dtime)[eid]; // Acessa o derivado do tempo da aresta
            // Relaxamento para o menor caminho
            if (VECTOR(result->dist_shortest_local)[u] + weight < VECTOR(result->dist_shortest_local)[v_node]) {
                
                VECTOR(result->dist_shortest_local)[v_node] = VECTOR(result->dist_shortest_local)[u] + weight;
                // Atualiza o predecessor para o menor caminho
                VECTOR(result->min_parents)[v_node] = u;
                if(v_node == to) result->min_cost = VECTOR(result->dist_shortest_local)[v_node];
                // Atualiza o derivado do menor caminho
                VECTOR(result->derivate_dist_shortest_local)[v_node] = VECTOR(result->derivate_dist_shortest_local)[u] + derivate_weight;
                if(v_node == to)result->derivate_min_cost = VECTOR(result->derivate_dist_shortest_local)[v_node];
            }

            // Relaxamento para o maior caminho
            // Só estender de 'u' se 'u' já tem um caminho mais longo válido (não -DBL_MAX)
            if ((VECTOR(result->dist_longest_local)[u] != -DBL_MAX) && (flow > 0)) {
                 if (VECTOR(result->dist_longest_local)[u] + weight > VECTOR(result->dist_longest_local)[v_node]) {
                    printf("Entrou!!!!!!!! %ld %ld\n",u,v_node);
                    VECTOR(result->dist_longest_local)[v_node] = VECTOR(result->dist_longest_local)[u] + weight;
                    VECTOR(result->max_parents)[v_node] = u;
                    if(v_node == to) result->max_cost = VECTOR(result->dist_longest_local)[v_node];

                    VECTOR(result->derivate_dist_longest_local)[v_node] = VECTOR(result->derivate_dist_longest_local)[u] + derivate_weight;
                    if(v_node == to) result->derivate_max_cost = VECTOR(result->derivate_dist_longest_local)[v_node];

                    if(flow < result->mu) result->mu = flow; // Atualiza o valor de mu se necessário


                }
            }
        }
        igraph_vector_int_clear(&incident_edges); // Reutiliza o vetor de arestas
    }
    reverse_parents(
        &result->min_parents, from, to
    );
    //igraph_vector_int_print(&result->min_parents);
    reverse_parents(
        &result->max_parents, from, to
    );
    //igraph_vector_int_print(&result->max_parents);
    //exit(0); // Encerra o programa após encontrar o caminho mínimo
    // Reverse min_parents and max_parents vectors for convenience
    //igraph_vector_int_print(&result->min_parents);
    //igraph_vector_int_print(&result->max_parents);
    // Print path from 'to' to 'from' using max_parents
    //printf("Caminho mais longo de %ld até %ld: ", from, to);
    /* igraph_integer_t current = to;
    while (current != from && current != -1) {
        igraph_integer_t prev = VECTOR(result->min_parents)[current];
        if (prev != -1) {
            // Find edge id between prev and current
            igraph_integer_t eid;
            igraph_get_eid(&bush->grafo, &eid, prev, current, IGRAPH_DIRECTED, 0);
            printf("%ld <- %ld (flow: %.2f)\n", prev,current, VECTOR(bush->time)[eid]);
        }
        current = prev;
    }
    if (current == from) {
        printf("%ld\n", from);
    } else {
        printf("(caminho não encontrado)\n");
    } */
    //igraph_vector_int_print(&topological_order);
    // Define os custos finais para o nó 'to' na estrutura de resultado

    // Libera memória dos vetores locais
    igraph_vector_int_destroy(&topological_order);
    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_destroy(&time);

}

void init_bush(
    struct BUSH *bush, 
    igraph_t *Grafo, 
    int i,
    int** edge_list,
    struct PARAMETERS* BPR_PARAMETERS,
    struct OD_MATRIX* OD,
    igraph_vector_t *flow,
    bool with_assignments
) {
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_vector_int_list_t evecs;
    igraph_vector_int_t inbound,parents;
    igraph_vector_t time;
    
    igraph_vector_init(&time, BPR_PARAMETERS->L);
    igraph_vector_int_list_init(&evecs, 0);
    igraph_vector_int_init(&inbound, 0);
    igraph_vector_int_init(&parents, 0);
    
    int fonte = OD->Elementos[i].fonte,j,k,id;
    //printf("Fonte: %d!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", fonte);
    
    BPR(&time, BPR_PARAMETERS, flow);
    
    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL,&evecs,fonte,igraph_vss_all(),&time, IGRAPH_OUT,&parents,&inbound
    );

    //printf("%ld %ld\n",igraph_vector_int_list_size(&evecs),igraph_ecount(Grafo));
    igraph_vector_int_destroy(&inbound);
    
    igraph_es_t edges;
    igraph_vector_int_t edges_vec;
    igraph_vector_int_init(&edges_vec, 0);
    
    bush->is_ingraph = (bool*) calloc(BPR_PARAMETERS->L, sizeof(bool));
    
    for (k = 0; k < igraph_vector_int_list_size(&evecs); k++) {

        igraph_vector_int_t* vec = igraph_vector_int_list_get_ptr(&evecs, k);
        for (j = 0; j < igraph_vector_int_size(vec); j++) {
            id = VECTOR(*vec)[j];
            if(bush->is_ingraph[id]) continue; // Se já foi adicionado, pula
            
            igraph_vector_int_push_back(&edges_vec, id);
            //printf("Adicionando aresta %d: (%d,%d) - %f, %f\n", id, edge_list[id][0], edge_list[id][1],VECTOR(*flow)[id], VECTOR(time)[id]);
            bush->is_ingraph[id] = true; // Marca como adicionado
        }
    }
    //printf("Size of bush->flow: %ld\n", igraph_vector_size(&bush->flow));
    //igraph_vector_int_print(&edges_vec);

    igraph_vector_int_list_destroy(&evecs);
    igraph_es_vector(&edges, &edges_vec);
    
    igraph_subgraph_from_edges(
        Grafo, &bush->grafo, edges, false
    );
    //printf("Grafo Bush:\n");
    //printf("Number of edges in Bush: %ld\n", igraph_ecount(&bush->grafo));
    igraph_vector_int_destroy(&edges_vec);
    igraph_es_destroy(&edges);

    int L = igraph_ecount(&bush->grafo);

    igraph_vector_init(&bush->flow, L);
    igraph_vector_init(&bush->BPR_PARAMETERS.capacidade, L);
    igraph_vector_init(&bush->BPR_PARAMETERS.cost_time, L);
    igraph_vector_int_init(&bush->edge_id, L);

    int** edge_list_bush = (int**) malloc(L * sizeof(int*));
    for (int j = 0; j < L; j++) edge_list_bush[j] = (int*) malloc(2 * sizeof(int));

    bush->BPR_PARAMETERS.L = L;
    bush->BPR_PARAMETERS.N = BPR_PARAMETERS->N;


    int target,antecessor,index;
    igraph_eit_t edge_it;
    j = 0;
    igraph_eit_create(&bush->grafo, igraph_ess_all(IGRAPH_EDGEORDER_ID), &edge_it);
    for (; !IGRAPH_EIT_END(edge_it); IGRAPH_EIT_NEXT(edge_it)) {
        igraph_integer_t eid = IGRAPH_EIT_GET(edge_it);
        igraph_integer_t from, to;
        igraph_edge(&bush->grafo, eid, &from, &to);
        //printf("Edge %ld: %ld -> %ld\n", eid, from, to);
        edge_list_bush[eid][0] = from;
        edge_list_bush[eid][1] = to;
        index = find_id_log(from, to, edge_list,BPR_PARAMETERS->L);

        VECTOR(bush->edge_id)[j] = index; // Armazena o id da aresta na bush
        VECTOR(bush->flow)[j] = 0.0; // Inicializa o fluxo da arest
        VECTOR(bush->BPR_PARAMETERS.capacidade)[j] = VECTOR(BPR_PARAMETERS->capacidade)[index];
        VECTOR(bush->BPR_PARAMETERS.cost_time)[j] = VECTOR(BPR_PARAMETERS->cost_time)[index];

        j++;
    }
    //igraph_vector_int_print(&parents);
    igraph_eit_destroy(&edge_it);
    if(with_assignments){
        for(j = 0; j < bush->n_alvos; j++) {
            target = VECTOR(OD->Elementos[i].alvos)[j];
            //printf("Alvo %d: %d\n", j, target);
            antecessor = target;
            while(target != fonte) {
                antecessor = VECTOR(parents)[target];
                //printf("(%d , %d)\n", antecessor,target);
                index = find_id_log(antecessor, target, edge_list_bush,L);
                VECTOR(bush->flow)[index] += VECTOR(OD->Elementos[i].volumes)[j];
                target = antecessor;
            }
        }
        char filename[100];
        sprintf(filename, "bush_%d.txt", fonte);
        FILE* file = fopen(filename, "w");
        for(j = 0; j < L; j++) {
            fprintf(file,"%d %d %f\n", edge_list_bush[j][0], edge_list_bush[j][1], VECTOR(bush->flow)[j]);
        }
        fclose(file);
    } 
    igraph_vector_destroy(&time);
    igraph_vector_int_destroy(&parents);
    for(int j = 0; j < L; j++) {
       free(edge_list_bush[j]); // Libera cada linha do edge_list
    }
    free(edge_list_bush); // Libera o array de ponteiros

}

bool equilibrium(
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
    
    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL,NULL,fonte,igraph_vss_all(),&time, IGRAPH_OUT,NULL,&inbound
    );
    
    bool converged = true;
    double tempo = 0;
    int edge_id,alvo,j;

    for ( j = 0; j< bush->n_alvos; j++){
        alvo = VECTOR(OD->Elementos[i].alvos)[j];
        while(alvo!=fonte){
            edge_id = VECTOR(inbound)[alvo];
            tempo += VECTOR(time)[edge_id];
            //printf("Edge %d: (%ld,%ld) - %f\n", alvo, IGRAPH_FROM(Grafo, edge_id), IGRAPH_TO(Grafo, edge_id), tempo);
            alvo = IGRAPH_FROM(Grafo, edge_id);
        }
        if(fabs(tempo - anterior_bush->paths->max_cost) > EPSILON) {
            converged = false;
            printf("\n Caminho não convergiu para Bush %d,%d: Maior tempo: %f, Tempo calculado: %f\n\n", fonte,(int)VECTOR(OD->Elementos[i].alvos)[j], anterior_bush->paths->max_cost, tempo);
            converged = false;
            //exit(0); // Encerra o programa se não convergir
            break;
        }
        tempo = 0;
    }
    
    if(!converged) {
        printf("Bush %d não está em equilíbrio, recalculando...🔢🔢🔢🔢🔢🔢🔢🔢🔢\n", i+1);

        double* total_time = (double*) calloc(bush->BPR_PARAMETERS.N , sizeof(double));
        bush->is_ingraph = (bool*) calloc(BPR_PARAMETERS->L, sizeof(bool));
        bool* keep_flow = (bool*) calloc(BPR_PARAMETERS->L , sizeof(bool));
        igraph_vector_int_t edges_vec;
        igraph_vector_int_init(&edges_vec, 0);
        int from, to;
        bool check;
        
        for ( j = 0; j< bush->BPR_PARAMETERS.N; j++){
            if(j == fonte) continue; // Pula a fonte
            alvo = j;
            tempo = 0;
            while(alvo!=fonte){
                edge_id = VECTOR(inbound)[alvo];
                tempo += VECTOR(time)[edge_id];
                alvo = IGRAPH_FROM(Grafo, edge_id);
                bush->is_ingraph[edge_id] = true; // Marca a aresta como parte do grafo da bush
            }
            total_time[j] = tempo;
        }
        
        for ( j = 0; j< BPR_PARAMETERS->L; j++){
            from = edge_list[j][0];
            to = edge_list[j][1];
            edge_id = VECTOR(inbound)[from];
            check = (total_time[from] < total_time[to]);
            //printf("%d/%d\n", j,BPR_PARAMETERS->L);
            
            if(check) {
                if(anterior_bush->is_ingraph[edge_id]) {
                    
                    if(bush->is_ingraph[edge_id]){
                        //printf("Aresta %d (%d,%d) está na bush!!!\n", edge_id, from, to);
                        igraph_vector_int_push_back(&edges_vec, edge_id);
                        //printf("Aresta %d (%d,%d) está na bush!!!\n", edge_id, from, to);
                        bush->is_ingraph[edge_id] = true; // Marca a aresta como parte do grafo da bush
                        keep_flow[edge_id] = 1; // Marca a aresta para manter o fluxo
                    }
                    
                }
                else{
                    if(bush->is_ingraph[edge_id]){
                        //printf("Aresta %d (%d,%d) está na bush!!!\n", edge_id, from, to);
                        igraph_vector_int_push_back(&edges_vec, edge_id);
                        bush->is_ingraph[edge_id] = true; // Marca a aresta como parte do grafo da bush
                        keep_flow[edge_id] = 0; // Marca a aresta para manter o fluxo
                    }
                }
            }
        }
        igraph_es_t edges;
        igraph_es_vector(&edges, &edges_vec);
    
        igraph_subgraph_from_edges(
            Grafo, &bush->grafo, edges, false
        );
        //printf("Grafo Bush:\n");
        //printf("Number of edges in Bush: %ld\n", igraph_ecount(&bush->grafo));
        bush->BPR_PARAMETERS.L = igraph_ecount(&bush->grafo);
        int L = bush->BPR_PARAMETERS.L;
        igraph_vector_init(&bush->flow, L);
        igraph_vector_print(&bush->flow);
        printf("\n");
        igraph_vector_init(&bush->BPR_PARAMETERS.capacidade, L);
        igraph_vector_init(&bush->BPR_PARAMETERS.cost_time, L);
        igraph_vector_int_init(&bush->edge_id, L);

        igraph_eit_t edge_it;
        int id;
        j = 0;
        igraph_eit_create(&bush->grafo, igraph_ess_all(IGRAPH_EDGEORDER_ID), &edge_it);

        for (; !IGRAPH_EIT_END(edge_it); IGRAPH_EIT_NEXT(edge_it)) {

            igraph_integer_t eid = IGRAPH_EIT_GET(edge_it);
            igraph_integer_t from, to;
            igraph_edge(&bush->grafo, eid, &from, &to);

            edge_id = find_id_log(from, to, edge_list,BPR_PARAMETERS->L);

            VECTOR(bush->edge_id)[j] = edge_id; // Armazena o id da aresta na bush

            id = VECTOR(anterior_bush->edge_id)[j];

            
            igraph_integer_t new_edge_id;
            igraph_get_eid(&anterior_bush->grafo, &new_edge_id, from, to, IGRAPH_DIRECTED, 0);
            VECTOR(bush->flow)[j] = keep_flow[edge_id] ? VECTOR(anterior_bush->flow)[new_edge_id] : 0.0;
            //printf("Aresta %d: (%ld,%ld) - %f\n", j, from, to, VECTOR(anterior_bush->flow)[new_edge_id]);
            VECTOR(bush->BPR_PARAMETERS.capacidade)[j] = VECTOR(BPR_PARAMETERS->capacidade)[edge_id];
            VECTOR(bush->BPR_PARAMETERS.cost_time)[j] = VECTOR(BPR_PARAMETERS->cost_time)[edge_id];
            //printf("%ld,%ld %f %f\n",from,to,VECTOR(BPR_PARAMETERS->cost_time)[edge_id],VECTOR(anterior_bush->BPR_PARAMETERS.cost_time)[new_edge_id]);
            j++;

        }
        
        printf("\nFluxo anterior: ");
        igraph_vector_print(&anterior_bush->flow);
        printf("Fluxo novo: ");
        igraph_vector_print(&bush->flow);
        printf("\n");
        free(total_time);

        igraph_vector_int_destroy(&edges_vec);
        igraph_es_destroy(&edges);
        igraph_eit_destroy(&edge_it);
        // Copy new bush to anterior_bush
        copy_bush(
            anterior_bush,
            bush
        );
        //exit(0); // Encerra o programa se não convergir
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

        //printf("Iniciando Bush %d\n",i+1);

        bushes[i].n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
        bushes[i].paths = (struct min_max_bush*) malloc(bushes[i].n_alvos * sizeof(struct min_max_bush));
        init_bush(&bushes[i], Grafo, i,edge_list, BPR_PARAMETERS,OD,solucao,true);
        bushes[i].steps = (double*) calloc(bushes[i].n_alvos , sizeof(double));

        if (check_inconsistent_edges(&bushes[i].grafo, &bushes[i].edge_id, edge_list)) {
            printf("Bush %d has inconsistent edges!\n", i+1);
            exit(1);
        }
    }

    bool IS_EPSLION_GLOBAL_CONVERGED = false;

    while (!IS_EPSLION_GLOBAL_CONVERGED){
        
        iter++;
        if(iter > 2) exit(0); // Encerra o programa se não convergir após 3 iterações
        IS_EPSLION_GLOBAL_CONVERGED = true; // Assume que todas as bushes estão convergidas
        printf("\n\n🔢🔢🔢🔢🔢🔢🔢🔢 Iteração %d 🔢🔢🔢🔢🔢🔢🔢🔢\n", iter);
        // Verifica se todas as bushes convergiram

        for ( i = 0; i < OD->size; i++){

            igraph_vector_fill(solucao, 0.0);
            converged_bush = 0;
            fonte = OD->Elementos[i].fonte;
            bool* bush_is_converged = (bool*) calloc(bushes[i].n_alvos , sizeof(bool));
            // ε-equilibrada
            
            while(converged_bush < bushes[i].n_alvos){
                converged_bush = 0; // Reseta o contador de bushes convergidas
                for ( j = 0; j < bushes[i].n_alvos; j++) {
                    
                    alvo = VECTOR(OD->Elementos[i].alvos)[j];
                    if(bush_is_converged[j]){
                        converged_bush++;
                        continue; // Se já está convergida, pula para o próximo alvo
                    }
                    find_dag_shortest_longest_costs_and_parents(
                        &bushes[i], fonte, alvo, &bushes[i].paths[j]
                    );
                    //exit(0);
                    printf("\nCusto mínimo: %f(%f), Custo máximo: %f(%f)\n", bushes[i].paths[j].min_cost,bushes[i].paths[j].derivate_min_cost, bushes[i].paths[j].max_cost,bushes[i].paths[j].derivate_max_cost);
                    bush_is_converged[j] = (bushes[i].paths[j].max_cost - bushes[i].paths[j].min_cost < EPSILON);
                    printf("%f\n",bushes[i].paths[j].max_cost);
                    if(!bush_is_converged[j]){
                        printf("\nCaminho com custo diferente encontrado entre %d e %d\n", fonte, alvo);
                        printf("\nCusto mínimo: %f(%f), Custo máximo: %f(%f)\n", bushes[i].paths[j].min_cost,bushes[i].paths[j].derivate_min_cost, bushes[i].paths[j].max_cost,bushes[i].paths[j].derivate_max_cost);
                        printf("Pais mínimos: ");
                        igraph_vector_int_print(&bushes[i].paths[j].min_parents);
                        printf("Pais máximos: ");
                        igraph_vector_int_print(&bushes[i].paths[j].max_parents);
                        exit(0); // Encerra o programa se não convergir
                        if(bushes[i].paths[j].max_cost,bushes[i].paths[j].derivate_max_cost ==0){
                            printf("Custo máximo é zero, encerrando o programa.\n");
                            exit(0); // Encerra o programa se o custo máximo for zero
                        }

                        int node1 = 0,node2 = 0,l,k = 0;
                        node1 = VECTOR(bushes[i].paths[j].min_parents)[k];
                        node2 = VECTOR(bushes[i].paths[j].max_parents)[k];
                        while(VECTOR(bushes[i].paths[j].min_parents)[k+1] == VECTOR(bushes[i].paths[j].max_parents)[k+1]){
                            node1 = VECTOR(bushes[i].paths[j].min_parents)[k];
                            node2 = VECTOR(bushes[i].paths[j].max_parents)[k];
                            k++;
                        }

                        //printf("Nó divergente: %d %d - %d\n", node1,node2,k);
                        double Denominator = (bushes[i].paths[j].derivate_max_cost + bushes[i].paths[j].derivate_min_cost);
                        double Numerator = bushes[i].paths[j].max_cost - bushes[i].paths[j].min_cost;
                        double step = bushes[i].steps[j] + Numerator / Denominator;
                        step = fmin(step, bushes[i].paths[j].mu); // Limita o passo ao valor de mu
                        if(step < 0) {
                            printf("Passo negativo encontrado: %f, ajustando para 0\n", step);
                            exit(0); // Encerra o programa se o passo for negativo
                            step = 0.0; // Ajusta o passo para 0 se for negativo
                        }
                        igraph_integer_t new_edge_id;

                        for(l = k+1; l < igraph_vector_int_size(&bushes[i].paths[j].min_parents); l++) {
                            node1 = VECTOR(bushes[i].paths[j].min_parents)[l-1];
                            node2 = VECTOR(bushes[i].paths[j].min_parents)[l];
                            igraph_get_eid(&bushes[i].grafo, &new_edge_id, node1, node2, IGRAPH_DIRECTED, 0);
                            VECTOR(bushes[i].flow)[new_edge_id] += step;
                        }

                        for(l = k+1; l < igraph_vector_int_size(&bushes[i].paths[j].max_parents); l++) {
                            node1 = VECTOR(bushes[i].paths[j].max_parents)[l-1];
                            node2 = VECTOR(bushes[i].paths[j].max_parents)[l];
                            igraph_get_eid(&bushes[i].grafo, &new_edge_id, node1, node2, IGRAPH_DIRECTED, 0);
                            VECTOR(bushes[i].flow)[new_edge_id] -= step;
                        }
                        
                        //igraph_vector_print(&bushes[i].flow);
                        bushes[i].steps[j] = step;
                        if(i == 0) printf("\nDenominator: %e, Numerator: %e, step: %e %f\n", Denominator, Numerator, bushes[i].steps[j],bushes[i].paths[j].mu);
                        if(i == 0) igraph_vector_print(&bushes[i].flow);

                    }
                    else converged_bush += bush_is_converged[j];
                    erase_path(&bushes[i].paths[j]); // Libera a memória da estrutura de caminho
                }
            }
            free(bush_is_converged);
            //printf("\n👌👌👌👌👌👌👌👌👌👌👌👌 Bush %d está em equilíbrio 👌👌👌👌👌👌👌👌👌👌👌👌\n\n", i+1);
            
            
            //exit(0); // Se não encontrar caminho, encerra o programa
            //printf("Bush %d convergiu com %d alvos\n", i, n_alvos);

            //ε-ótima
            //igraph_vector_print(&bushes[i].flow);
            int id;
            for (int k = 0; k < i+1; k++){
                for ( j = 0; j < bushes[k].BPR_PARAMETERS.L; j++){
                    id = VECTOR(bushes[k].edge_id)[j];
                    VECTOR(*solucao)[id] += VECTOR(bushes[k].flow)[j];
                    if(VECTOR(*solucao)[id] < 0) {
                        printf("Fluxo negativo encontrado na aresta %d: %f\n", id, VECTOR(*solucao)[id]);
                        exit(0); // Encerra o programa se encontrar fluxo negativo
                    }

                }
            }
            
            struct BUSH new_bushes;
            new_bushes.n_alvos = bushes[i].n_alvos;
            new_bushes.BPR_PARAMETERS.N = BPR_PARAMETERS->N;
            new_bushes.paths = (struct min_max_bush*) malloc(new_bushes.n_alvos * sizeof(struct min_max_bush));

            printf("Solução: \n");
            igraph_vector_print(solucao);

            bool s = equilibrium( &new_bushes,&bushes[i], Grafo, i, edge_list, BPR_PARAMETERS, OD, solucao);
            IS_EPSLION_GLOBAL_CONVERGED = IS_EPSLION_GLOBAL_CONVERGED && s;
            if (check_inconsistent_edges(&bushes[i].grafo, &bushes[i].edge_id, edge_list)) {
                printf("Bush %d has inconsistent edges!\n", i+1);
                exit(1);
            }
            //exit(0); // Encerra o programa se não convergir após 3 iterações
        }
    }
}