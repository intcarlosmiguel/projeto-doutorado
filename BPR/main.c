#include "bib/simulate.h"

int main(int argc,char *argv[ ]){
    int arquivo = atoi(argv[1]);
    int seed = atoi(argv[2]);
    simulate(arquivo,seed);
}