#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "../include/tsp.h"
#include "../include/graph.h"


int     main(int argc, char **argv){
    /* argv[1] = number of cities      *
     * argv[2] = cities file           */
    int N;
    unsigned int start, finish, elapsed;
    time_t t; srand((unsigned) time(&t));
    
    if( argc!=3 ){
        fprintf(stderr, "[Error]: Not enough/too many arguments:\n");
        fprintf(stderr, "         1: number of cities\n");
        fprintf(stderr, "         2: cities file\n");
        return 0;
    }
    N = atoi(argv[1]);
    init_seq(N, argv[2]);
    init_distances(N);
    g_init();
    start = get_clock_count();

    while( !done_flag ){
        calculate_fitnesses(N);
        normalize_fitnesses();
        // print_total_fitness();
        update_best_offspring(N);
        // print_pops_and_fits(N);
        create_new_generation(N);
        printf("generation count = %lu\n", gen_count);
        // print_best_offspring(N);
        // printf("\n********************************************************\n\n");
        g_show_salesman_paths(N);
    }
    finish  = get_clock_count();
    elapsed = finish - start;
    print_stats(N, elapsed);
    
    return 0;
}