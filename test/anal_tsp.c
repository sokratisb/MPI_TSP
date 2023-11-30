#include <stdio.h>
#include <stdlib.h>
#include "../include/tsp.h"

int N;
int  *best_route;
float shortest_distance;

void    swap(int *x, int *y){ 
    int temp; 
    temp = *x; 
    *x = *y; 
    *y = temp; 
} 

void    anal_tsp(int *a, int l, int r){ 
    int i;
    float curr_dist = 0.0;
    if (l == r){
        a[N] = a[0];
        for( i=0; i<N+1; i++){
            printf("%d ", a[i]);
        }
        printf("\n");

        for( i=0; i<N; i++ ){
            if( i==(N-1) )
                curr_dist += read_distance(list_of_cities[a[i]], list_of_cities[a[0]]);
            else
                curr_dist += read_distance(list_of_cities[a[i]], list_of_cities[a[i+1]]);
        }
        if( curr_dist<shortest_distance ){
            shortest_distance = curr_dist;
            for( i=0; i<N; i++ )
                best_route[i] = a[i];
            best_route[i] = a[0];
        }
    }
    else{
        for (i = l; i <= r; i++){
            swap((a + l), (a + i));
            anal_tsp(a, l + 1, r);
            swap((a + l), (a + i));
        }
    }
}

void    init_best_route(int N){
    int i;
    best_route = malloc((N+1) * sizeof(int));
    for( i=0; i<N; i++ )
        best_route[i] = i;
    best_route[i] = best_route[0];
    shortest_distance = 0.0;
    for( i=0; i<N; i++ ){
        if( i==(N-1) )
            shortest_distance += read_distance(list_of_cities[i], list_of_cities[0]);
        else
            shortest_distance += read_distance(list_of_cities[i], list_of_cities[i+1]);
    }
}


void    init(int N, char *file){
    init_cities(N, file);
    init_distances(N);
    init_best_route(N);
}

void    print_best_route(unsigned int exec_time){
    int i;
    printf("\n\n*************************************************\n");
    printf("*\tExecution time :\t%u (sec)\n", exec_time/1000);
    printf("*\tBest route:\t\t\t\t*\n");
    for( i=0; i<N+1; i++ ){
        printf("*\t\t%d ", best_route[i]);
        if( i!=N )
            printf("\t>\t\t\t*\n");
    }
    printf("\t\t\t\t*");
    printf("\n*\tDistance = %f\t\t\t*\n", shortest_distance);
    printf("*************************************************\n\n");
}


int     main(int argc, char **argv){
    /* argv[1] = number of cities      *
     * argv[2] = cities file           */
    int i;
    unsigned int start, finish, elapsed;
    int *aux_array = malloc((N+1) * sizeof(int));

    N = atoi(argv[1]);
    init(N, argv[2]);

    for( i=0; i<N; i++ )
        aux_array[i] = list_of_cities[i]->ID;
    aux_array[i] = aux_array[0];

    start = get_clock_count();
    anal_tsp(aux_array, 0, N-1);
    finish  = get_clock_count();
    elapsed = finish - start;
    print_best_route(elapsed);
    return 0;
}