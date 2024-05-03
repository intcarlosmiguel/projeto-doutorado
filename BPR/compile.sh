gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++ -Wall -Werror

for i in {1..121}; do
    ./main $i
done