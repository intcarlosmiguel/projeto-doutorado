gcc -I/home/carlos/igraph/build/include -I/home/carlos/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++ 
./main "./file/example/edges.txt" "./file/example/dial_matod.txt"