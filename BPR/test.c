#include <igraph/igraph.h>
#include <stdio.h>
#include <stdlib.h> // Para rand() e srand()
#include <time.h>   // Para time()

int main() {
    igraph_t g;
    igraph_vector_t pesos_vector; // Vetor para armazenar os pesos (números reais)
    igraph_eit_t eit;             // Iterador de Arestas (Edge Iterator)

    // Inicializa a tabela de atributos, necessário para usar "weight"
    igraph_set_attribute_table(&igraph_cattribute_table);

    // Inicializa o gerador de números aleatórios
    srand(time(NULL));

    int num_vertices = 10;
    double prob_aresta = 0.4; // 40% de chance de cada par de vértices ter uma aresta

    printf("Criando um grafo de Erdős-Rényi com %d vértices e probabilidade %.2f.\n", num_vertices, prob_aresta);

    // 1. Criar o grafo de Erdős-Rényi
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, num_vertices, prob_aresta, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    long int num_arestas = igraph_ecount(&g);
    printf("O grafo foi gerado com %ld arestas.\n", num_arestas);

    // 2. Preparar e atribuir os pesos às arestas
    if (num_arestas > 0) {
        printf("Atribuindo pesos aleatórios entre 1.0 e 10.0 a cada aresta...\n\n");

        // Inicializa o vetor de pesos com o tamanho exato do número de arestas
        igraph_vector_init(&pesos_vector, num_arestas);

        // Preenche o vetor com pesos aleatórios
        for (long int i = 0; i < num_arestas; i++) {
            // Gera um double aleatório no intervalo [1.0, 10.0]
            double peso_aleatorio = 1.0 + ((double)rand() / RAND_MAX) * 9.0;
            VECTOR(pesos_vector)[i] = peso_aleatorio;
        }

        // Define o atributo "weight" para TODAS as arestas do grafo de uma só vez
        // SET_EAS = Set Edge Attribute (Numeric)
        SETEANV(&g, "weight", &pesos_vector);

        // O vetor de pesos foi copiado para o grafo, podemos limpá-lo
        igraph_vector_destroy(&pesos_vector);
    }

    // 3. Imprimir as arestas e seus pesos para verificação
    printf("Listando as arestas do grafo e seus pesos:\n");

    // Cria um iterador para percorrer todas as arestas
    igraph_eit_create(&g, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);

    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t edge_id = IGRAPH_EIT_GET(eit); // Pega o ID da aresta
        igraph_integer_t from_id, to_id;

        // Pega os IDs dos vértices que a aresta conecta
        igraph_edge(&g, edge_id, &from_id, &to_id);

        // Pega o valor do atributo "weight" para esta aresta específica (pelo seu ID)
        // EAN = Edge Attribute (Numeric)
        double peso = EAN(&g, "weight", edge_id);

        // Imprime a aresta e seu peso formatado com duas casas decimais
        printf("  Aresta: %ld -- %ld,  Peso: %.2f\n", (long)from_id, (long)to_id, peso);

        IGRAPH_EIT_NEXT(eit); // Avança para a próxima aresta
    }

    // 4. Limpeza da memória
    igraph_eit_destroy(&eit);
    igraph_destroy(&g);

    printf("\nGrafo destruído com sucesso.\n");

    return 0;
}