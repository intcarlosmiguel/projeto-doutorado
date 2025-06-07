#pragma once

#include <igraph.h>
#include "network.h"

void find_dag_shortest_longest_costs_and_parents(
    struct BUSH *bush,
    igraph_integer_t from,
    igraph_integer_t to,
    struct min_max_bush *result
) {
    igraph_integer_t N = igraph_vcount(&bush->grafo);

    // Garante que os vetores de saída estejam dimensionados e inicializados
    igraph_vector_int_init(&result->min_parents, N);
    igraph_vector_int_init(&result->max_parents, N);
    igraph_vector_int_fill(&result->min_parents, -1); // -1 indica nenhum predecessor
    igraph_vector_int_fill(&result->max_parents, -1);

    // Vetores locais para armazenar distâncias para todos os nós
    igraph_vector_t dist_shortest_local;
    igraph_vector_t dist_longest_local;
    igraph_vector_init(&dist_shortest_local, N);
    igraph_vector_fill(&dist_shortest_local, DBL_MAX);
    igraph_vector_init(&dist_longest_local, N);
    igraph_vector_fill(&dist_longest_local, -DBL_MAX); // "Infinito negativo"
    VECTOR(dist_shortest_local)[from] = 0.0;
    VECTOR(dist_longest_local)[from] = 0.0;

    igraph_vector_int_t topological_order;
    igraph_vector_int_init(&topological_order, 0); // Inicializa para evitar problemas no destroy
    if (igraph_topological_sorting(&bush->grafo, &topological_order, IGRAPH_OUT) != IGRAPH_SUCCESS) {
        printf("Erro ao ordenar topologicamente o grafo.\n");
        exit(0); // Falha ao ordenar topologicamente, encerra o programa
    }
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);
    igraph_vector_print(&bush->flow);
    for (int i = 0; i < igraph_vector_int_size(&topological_order); ++i) {
        igraph_integer_t u = VECTOR(topological_order)[i];
        // Se 'u' não é alcançável a partir de 'from' (para menor caminho),
        // não pode contribuir para caminhos subsequentes.
        if (VECTOR(dist_shortest_local)[u] == DBL_MAX) continue;

        igraph_incident(&bush->grafo, &incident_edges, u, IGRAPH_OUT);
        for (long j = 0; j < igraph_vector_int_size(&incident_edges); ++j) {
            igraph_integer_t eid = VECTOR(incident_edges)[j];
            igraph_integer_t v_node = IGRAPH_TO(&bush->grafo, eid);
            igraph_real_t weight;

            weight = VECTOR(bush->flow)[eid]; // Acessa o peso da aresta
            // Relaxamento para o menor caminho
            if (VECTOR(dist_shortest_local)[u] + weight < VECTOR(dist_shortest_local)[v_node]) {
                VECTOR(dist_shortest_local)[v_node] = VECTOR(dist_shortest_local)[u] + weight;
                // Atualiza o predecessor para o menor caminho
                VECTOR(result->min_parents)[v_node] = u;
                if(v_node == to) result->min_cost = VECTOR(dist_shortest_local)[v_node];
                
            }

            // Relaxamento para o maior caminho
            // Só estender de 'u' se 'u' já tem um caminho mais longo válido (não -DBL_MAX)
            if (VECTOR(dist_longest_local)[u] != -DBL_MAX) {
                 if (VECTOR(dist_longest_local)[u] + weight > VECTOR(dist_longest_local)[v_node]) {
                    VECTOR(dist_longest_local)[v_node] = VECTOR(dist_longest_local)[u] + weight;
                    // Atualiza o predecessor para o maior caminho
                    VECTOR(result->max_parents)[v_node] = u;
                    if(v_node == to) result->max_cost = VECTOR(dist_longest_local)[v_node];
                }
            }
        }
        igraph_vector_int_clear(&incident_edges); // Reutiliza o vetor de arestas
    }
    igraph_vector_int_print(&result->min_parents);
    igraph_vector_int_print(&result->max_parents);
    igraph_vector_int_print(&topological_order);
    printf("Custo mínimo: %f\n", result->min_cost);
    printf("Custo máximo: %f\n", result->max_cost);
    // Define os custos finais para o nó 'to' na estrutura de resultado

    // Libera memória dos vetores locais
    igraph_vector_destroy(&dist_shortest_local);
    igraph_vector_destroy(&dist_longest_local);
    igraph_vector_int_destroy(&topological_order);
    igraph_vector_int_destroy(&incident_edges);

}

void init_bush(struct BUSH *bush, igraph_t *Grafo, igraph_vector_t *pesos, int fonte) {
    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t inbound,parents;
    
    igraph_vector_int_list_init(&vecs, 0);
    igraph_vector_int_list_init(&evecs, 0);
    igraph_vector_int_init(&inbound, 0);
    igraph_vector_int_init(&parents, 0);
    igraph_vector_init(&bush->flow,0);
    

    igraph_get_shortest_paths_dijkstra(
        Grafo, &vecs,&evecs,fonte,igraph_vss_all(),pesos, IGRAPH_OUT,&parents,&inbound
    );

    printf("%ld %ld\n",igraph_vector_int_list_size(&evecs),igraph_ecount(Grafo));
    igraph_vector_int_destroy(&inbound);
    igraph_vector_int_destroy(&parents);

    igraph_es_t edges;
    igraph_vector_int_t edges_vec;
    igraph_vector_int_init(&edges_vec, 0);
    int id;
    for (long int i = 0; i < igraph_vector_int_list_size(&evecs); i++) {
        igraph_vector_int_t* vec = igraph_vector_int_list_get_ptr(&evecs, i);
        for (long int j = 0; j < igraph_vector_int_size(vec); j++) {

            id = VECTOR(*vec)[j];
            igraph_vector_int_push_back(&edges_vec, id);
            
        }
    }
    igraph_vector_int_print(&edges_vec);

    igraph_es_vector(&edges, &edges_vec);

    igraph_subgraph_from_edges(
        Grafo, &bush->grafo, edges, false
    );
    printf("Grafo Bush:\n");
    printf("Number of edges in Bush: %ld\n", igraph_ecount(&bush->grafo));
    
    igraph_vector_int_list_destroy(&vecs);
    igraph_vector_int_list_destroy(&evecs);
    igraph_vector_int_destroy(&edges_vec);
    igraph_es_destroy(&edges);

    exit(0);

}

void Bushes(igraph_t *Grafo,struct OD_MATRIX* OD,int** edge_list,struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t *solucao) {
    int i,c = 0,alvo,fonte,j,converged_bush;
    igraph_vector_t fluxo,pesos;
    igraph_vector_init(&pesos, BPR_PARAMETERS->L);
    igraph_vector_fill(&pesos, 0);

    struct BUSH* bushes = (struct BUSH*) malloc(OD->size * sizeof(struct BUSH));

    for ( i = 0; i < OD->size; i++){
        printf("Iniciando Bush %d\n",i+1);
        fonte = OD->Elementos[i].fonte;
        bushes[i].n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
        init_bush(&bushes[i], Grafo, &pesos, i);
        igraph_vector_print(&bushes[i].flow);
        
    }
    bool IS_EPSLION_GLOBAL_CONVERGED = false;
    while (!IS_EPSLION_GLOBAL_CONVERGED){
        for ( i = 0; i < OD->size; i++){

            bool* bush_is_converged = (bool*) calloc(bushes[i].n_alvos , sizeof(bool));
            converged_bush = 0;
            fonte = OD->Elementos[i].fonte;
            while(converged_bush < bushes[i].n_alvos){
                for ( j = 0; j < bushes[i].n_alvos; j++) {

                    if(bush_is_converged[j]) continue;
                    struct min_max_bush path;
                    alvo = VECTOR(OD->Elementos[i].alvos)[j];
    
                    find_dag_shortest_longest_costs_and_parents(
                        &bushes[i], fonte, alvo, &path
                    );

                    bush_is_converged[j] = (path.max_cost - path.min_cost < EPSILON);
                    if(!bush_is_converged[j]){
                        printf("Caminho com custo diferente encontrado entre %d e %d\n", fonte, alvo);
                        printf("Custo mínimo: %f, Custo máximo: %f\n", path.min_cost, path.max_cost);
                        igraph_vector_int_print(&path.min_parents);
                        igraph_vector_int_print(&path.max_parents);
                    }
                    else converged_bush++;
                }
            }
            free(bush_is_converged);
            printf("Bush %d convergiu!!!!!!!!!!!!!!!!!\n", i+1);
            
            

            igraph_vector_destroy(&bushes[i].flow);
            igraph_destroy(&bushes[i].grafo);
            //exit(0); // Se não encontrar caminho, encerra o programa
            //printf("Bush %d convergiu com %d alvos\n", i, n_alvos);
    
        }
        
        struct BUSH* new_bushes = (struct BUSH*) malloc(OD->size * sizeof(struct BUSH));
        
        // Cria uma nova bush considerando a rede inteira
        /*for ( i = 0; i < OD->size; i++){
            printf("Iniciando Bush %d\n",i+1);
            fonte = OD->Elementos[i].fonte;
            new_bushes[i].n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
            init_bush(&new_bushes[i], Grafo, OD, edge_list, &solucao, i);
            bool is_converged = true;
            for ( j = 0; j < bushes[i].n_alvos; j++){
                struct min_max_bush path,path_complete;
                alvo = VECTOR(OD->Elementos[i].alvos)[j];

                find_dag_shortest_longest_costs_and_parents(
                    &bushes[i], fonte, alvo, &path
                );
                find_dag_shortest_longest_costs_and_parents(
                    &new_bushes[i], fonte, alvo, &path_complete
                );
                is_converged = (path.max_cost - path_complete.min_cost < EPSILON);
                if(!is_converged){
                    printf("Caminho com custo diferente encontrado entre %d e %d\n", fonte, alvo);
                    printf("Custo mínimo: %f, Custo máximo: %f\n", path.min_cost, path.max_cost);
                    igraph_vector_int_print(&path.min_parents);
                    igraph_vector_int_print(&path.max_parents);
                    igraph_vector_int_print(&path_complete.min_parents);
                    igraph_vector_int_print(&path_complete.max_parents);
                    break;
                }
            }
            if(!is_converged){

            }
        }*/

        free(new_bushes); // Libera a memória dos bushes antigos
        exit(0); // Encerra o programa após processar todas as bushes
    }
    


    
}

