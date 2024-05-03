#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "igraph.h"
#include "bib/simulate.h"

int main(int argc,char *argv[ ]){
    int arquivo = atoi(argv[1]);
    simulate(arquivo);
}