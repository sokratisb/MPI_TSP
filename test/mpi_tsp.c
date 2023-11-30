#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include "../include/tsp.h"
#include "mpi.h"


#define MASTER 0

enum type{ INIT_TAG, RESULT_TAG, GEN_CNT_TAG };
MPI_Status status;

void    send_cities_to_slaves(int nprocs, int N){
    int dest, i;
    int size = N*sizeof(CITY);
    int *buffer = malloc(size);
    for( i=0; i<3*N; i+=3 ){
        buffer[i]   = list_of_cities[i/3]->ID;
        buffer[i+1] = list_of_cities[i/3]->x;
        buffer[i+2] = list_of_cities[i/3]->y;
    }

    for( i = 1; i < nprocs; i++ ){
        dest = i;
        MPI_Send(buffer, size, MPI_BYTE, dest, INIT_TAG, MPI_COMM_WORLD);
    }
}

void    receive_cities_from_master(int N){
    int i;
    int source = 0;
    int size = N*sizeof(CITY);
    int *buffer = malloc(size);

    list_of_cities = malloc(N*sizeof(CITY *));
    for( i=0; i<N; i++ )
        list_of_cities[i] = malloc(sizeof(CITY));

    MPI_Recv(buffer, size, MPI_BYTE, source, INIT_TAG, MPI_COMM_WORLD, &status);
    for( i=0; i<3*N; i+=3 ){
        list_of_cities[i/3]->ID = buffer[i];
        list_of_cities[i/3]->x  = buffer[i+1];
        list_of_cities[i/3]->y  = buffer[i+2];
    }
}

void    master_init(int N, int nprocs, char *file){
    init_cities(N, file);
    init_distances(N);
    send_cities_to_slaves(nprocs, N);
    init_gens(N);
    shuffle_populations(N);
}

void    receive_best_off_from_slaves(int N){
    int i;
    int size = (N+1)*3*sizeof(int) + sizeof(float) + sizeof(int) ;
    float *buffer = malloc(size);
    float slave_fitness;
    int gen_increment;
    int recv_flag;

    MPI_Iprobe(MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &recv_flag, &status);
    while(recv_flag){
        MPI_Recv(buffer, size, MPI_BYTE, MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
        slave_fitness = buffer[3*(N+1)];
        gen_increment = (int)buffer[3*(N+1)+1];
        printf("[MASTER]: Received new (local) best from {SLV_%d}:\n\t  ", status.MPI_SOURCE);
        for( i=0; i<3*(N+1); i+=3){
            printf("%d", (int)buffer[i]);
            if( i!=(3*(N+1)-1) )
                printf(">");
        }
        printf("\n[MASTER]: Arrived fit = %f -- best = %f -- gen incr = %d\n", slave_fitness, best_fitness, gen_increment);
        if( slave_fitness > best_fitness ){
            best_fitness = slave_fitness;
            for( i=0; i<(3*(N+1)); i+=3 ){
                best_offspring[i/3]->ID = (int)buffer[i];
                best_offspring[i/3]->x  = (int)buffer[i+1];
                best_offspring[i/3]->y  = (int)buffer[i+2];
                gen_count = 0;
                done_flag = 0;
            }
        }
        else{
           gen_count += gen_increment;
           total_generations += gen_increment;
        }
        MPI_Iprobe(MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &recv_flag, &status);
    }
}

void    check_for_gen_cnt_messages(void){
    int recv_flag;
    int dummy;

    MPI_Iprobe(MPI_ANY_SOURCE, GEN_CNT_TAG, MPI_COMM_WORLD, &recv_flag, &status);
    while(recv_flag){
        MPI_Recv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, GEN_CNT_TAG, MPI_COMM_WORLD, &status);
        gen_count += 100;
        total_generations += 100;
        printf("[MASTER]: Generation count +100 => gen_cnt = %lu\n", gen_count);
        MPI_Iprobe(MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &recv_flag, &status);
    }
}

void    conditionally_print_new_best(int N){
    if( gen_count==0 )
        print_best_offspring(N);
}

void    master_loop(int N){
    while( !done_flag ){
        calculate_fitnesses(N);
        normalize_fitnesses();
        update_best_offspring(N);
        conditionally_print_new_best(N);
        create_new_generation(N);
        receive_best_off_from_slaves(N);
        check_for_gen_cnt_messages();
        total_generations++;
        average_generations++;
    }
}

void    send_new_best_to_master(int N, int mynum){
    int i;
    int dest = 0;
    int size = (N+1)*3*sizeof(int) + sizeof(float) + sizeof(int);
    float *buffer = malloc(size);
    unsigned long gen_increment = (total_generations - last_update_gen_count)%101;

    for( i=0; i<(3*(N+1)); i+=3 ){
        buffer[i+0] = best_offspring[i/3]->ID;
        buffer[i+1] = best_offspring[i/3]->x;
        buffer[i+2] = best_offspring[i/3]->y;
    }
    buffer[3*(N+1)]   = best_fitness;
    buffer[3*(N+1)+1] = gen_increment;
    printf("[SLV_%d]:  Sending to master: fit = %f, incr = %d\n", mynum, buffer[3*(N+1)], (int)buffer[3*(N+1)+1]);
    MPI_Send(buffer, size, MPI_BYTE, dest, RESULT_TAG, MPI_COMM_WORLD);
}

void    conditionally_send_gen_count_msg(void){
    unsigned long diff = total_generations - last_update_gen_count;
    int dummy = 0; int dest = 0;

    if( diff % 101==100 ){
        MPI_Send(&dummy, 1, MPI_INT, dest, GEN_CNT_TAG, MPI_COMM_WORLD);
    }

}

void    slave_init(int N){
    receive_cities_from_master(N);
    init_distances(N);
    init_gens(N);
    shuffle_populations(N);
}

void    slave_loop(int N, int mynum){
    int i;
    int update_best_off_res;
    while(1){
        calculate_fitnesses(N);
        normalize_fitnesses();
        update_best_off_res = update_best_offspring(N);
        if( update_best_off_res ){
            send_new_best_to_master(N, mynum);
            last_update_gen_count = total_generations;
        }
        else
            conditionally_send_gen_count_msg();
        create_new_generation(N);
        total_generations++;
    }
}

int     main(int argc, char **argv){
    int mynum, nprocs, N;
    unsigned int start, finish, elapsed;
    time_t t; srand((unsigned) time(&t) + (unsigned)getpid());

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &mynum );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    if( argc!=3 ){
        fprintf(stderr, "[Error]: Not enough/too many arguments:\n");
        fprintf(stderr, "         1: number of cities\n");
        fprintf(stderr, "         2: cities file\n");
        return 0;
    }
    N = atoi(argv[1]);

    if( mynum==MASTER){
        master_init(N, nprocs, argv[2]);
        start = get_clock_count();
        master_loop(N);        
    }
    else{
        slave_init(N);
        slave_loop(N, mynum);
    }
    finish  = get_clock_count();
    elapsed = finish - start;
    print_stats(N, elapsed);
    int errorcode;
    MPI_Abort(MPI_COMM_WORLD, errorcode);
}