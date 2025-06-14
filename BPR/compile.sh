gcc -I./igraph/build/include -I./igraph/include main.c -o main -Ibib -L./igraph/build/src -ligraph -fopenmp -O3 -lstdc++ -lm
./main "./file/example/edges.txt" "./file/example/dial_matod.txt"
#gcc test.c -o test -Ibib -ligraph -fopenmp -O3 -lstdc++ -lm
#./test