#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "igraph.h"


int get_file_size(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    int lines = 0;
    char ch;
    while (!feof(file)) {
        ch = fgetc(file);
        if (ch == '\n') {
            lines++;
        }
    }

    fclose(file);
    return lines;
}

int main() {
    int Nodes,Edges,i,j,u,v;
    double random,distance;
    int number_bairros = get_file_size("./OD/input/densidade.txt");
    Edges = get_file_size("./OD/input/edges_fortaleza.txt");
    Nodes = get_file_size("./OD/input/nodes_fortaleza.txt");
    int* densidade = (int*) malloc(number_bairros * sizeof(int));
    double* sizes = (double*) calloc(Nodes, sizeof(double));
    double* r = (double*) calloc(number_bairros, sizeof(double));
    int* bairros = (int*) calloc(Nodes, sizeof(int));

    igraph_vector_t length;
    igraph_vector_int_t edges;
    
    igraph_vector_init(&length, 0);
    igraph_vector_int_init(&edges, 0);
    
    srand(42);

    FILE* file = fopen("./OD/input/densidade.txt", "r");
    for (i = 0; i < number_bairros; i++) {
        if (fscanf(file, "%d\n",&v) != 1) {
            perror("Error reading file");
            fclose(file);
            exit(EXIT_FAILURE);
        }
        densidade[i] = v;
    }
    fclose(file);
    
    file = fopen("./OD/input/edges_fortaleza.txt", "r");
    for (i = 0; i < Edges; i++) {
        if (fscanf(file, "%d %d %lf\n", &u, &v, &distance) != 3) {
            perror("Error reading file");
            fclose(file);
            exit(EXIT_FAILURE);
        }
        igraph_vector_int_push_back(&edges, u);
        igraph_vector_int_push_back(&edges, v);
        igraph_vector_push_back(&length, distance);
    }
    fclose(file);

    
    FILE* nos = fopen("./OD/input/nodes_fortaleza.txt", "r");
    for (i = 0; i < Nodes; i++) {
        if (fscanf(nos, "%d %d\n",&u,&v) != 2) {
            perror("Error reading nos 2\n");
            fclose(nos);
            exit(EXIT_FAILURE);
        }
        random = ((double) rand() / RAND_MAX);
        r[v] += random;
        sizes[u] = random;
        bairros[u] = v;
    }
    for (i = 0; i < Nodes; i++) {
        sizes[i] = sizes[i]*densidade[bairros[i]]/r[bairros[i]];
    }
    free(r);
    free(densidade);
    free(bairros);
    fclose(nos);
    igraph_t Grafo;
    printf("Creating graph with %d nodes and %d edges\n", Nodes, Edges);
    igraph_empty(&Grafo, Nodes, IGRAPH_UNDIRECTED);
    
    //igraph_empty(&Grafo, 22, IGRAPH_UNDIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);

    igraph_vector_t shortest_paths;
    igraph_vector_init(&shortest_paths, 0);
    double diameter = 29797.872000;
    igraph_vector_t histogram;
    igraph_vector_init(&histogram, 0);
    igraph_vector_init_range(&histogram, 0, diameter+1);
    igraph_vector_null(&histogram);
    igraph_integer_t from, to;
    igraph_vector_int_t path, path_edge;
    igraph_vector_int_init(&path, 0);
    igraph_vector_int_init(&path_edge, 0);
    igraph_diameter_dijkstra(&Grafo, &length, &diameter, &from, &to,&path,&path_edge, IGRAPH_DIRECTED, IGRAPH_ALL);
    printf("diameter: %f %ld %ld\n", diameter,from,to);
    exit(0);
    /*for (i = 0; i < Nodes; i++) {
        igraph_matrix_t result;
        igraph_matrix_init(&result, 0, 0);
        igraph_vs_t vs;
        igraph_vs_1(&vs, i);
        igraph_distances_dijkstra(&Grafo, &result, vs, igraph_vss_all(),&length, IGRAPH_OUT);
        igraph_vs_destroy(&vs);
        for (j = 0; j < Nodes; j++){
            if(i == j) continue;
            distance = MATRIX(result, 0, j);
            if(MATRIX(result, 0, j) <= 0.) continue;
            if( MATRIX(result, 0, j) > diameter) continue;
            
            VECTOR(histogram)[(int)round(distance)]++;
            //printf("%e\n", (double) sizes[j]*sizes[i]/(distance*distance));
        }
        igraph_matrix_destroy(&result);
        printf("%d/%d\n", i+1, Nodes);

    }
    FILE* histogram_file = fopen("./OD/output/histogram.txt", "w");
    if (histogram_file == NULL) {
        perror("Unable to open file for writing");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < igraph_vector_size(&histogram); i++) {
        fprintf(histogram_file, "%d %f\n", i, VECTOR(histogram)[i]);
    }

    fclose(histogram_file);*/
    igraph_vector_destroy(&histogram);

    igraph_vector_destroy(&shortest_paths);


    igraph_destroy(&Grafo);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&length);
}
