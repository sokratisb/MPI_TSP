#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "../include/tsp.h"

void    init_cities(int N, char *file){
    int i;
    FILE *fp = fopen(file, "r");
    list_of_cities = malloc(N*sizeof(CITY *));
    for( i=0; i<N; i++ ){
        list_of_cities[i] = malloc(sizeof(CITY));
        fscanf(fp, "%d %d", &(list_of_cities[i]->x), &(list_of_cities[i]->y));
        list_of_cities[i]->ID = i;
    }
    fclose(fp);
}

void    init_gens(int N){
    int i, j;

    best_fitness   = -10.0; // Something small to compare the first calculated fitness with.
    fitness        = malloc(MAX_POPUL*sizeof(float));
    norm_fitness   = malloc(MAX_POPUL*sizeof(float));
    curr_gen       = malloc(MAX_POPUL*sizeof(CITY **));
    next_gen       = malloc(MAX_POPUL*sizeof(CITY **));
    best_offspring = malloc((N+1)*sizeof(CITY *));
    for( i=0; i<MAX_POPUL; i++ ){
        curr_gen[i] = malloc((N+1)*sizeof(CITY *));
        next_gen[i] = malloc((N+1)*sizeof(CITY *));
        for( j=0; j<N+1; j++ ){
            curr_gen[i][j]    = malloc(sizeof(CITY));
            next_gen[i][j]    = malloc(sizeof(CITY));
            best_offspring[j] = malloc(sizeof(CITY));
        }
    }
    for( i=0; i<MAX_POPUL; i++ ){
        for( j=0; j<N; j++ ){
            curr_gen[i][j]->x  = list_of_cities[j]->x;
            curr_gen[i][j]->y  = list_of_cities[j]->y;
            curr_gen[i][j]->ID = j;
        }
        curr_gen[i][j]->x  = list_of_cities[0]->x;
        curr_gen[i][j]->y  = list_of_cities[0]->y;
        curr_gen[i][j]->ID = j;
    }
    gen_count = 0;
    done_flag = 0;
    best_distance = 0.0;
    total_generations = 0;
    average_generations = 0;
    last_update_gen_count = 0;
}

void    init_distances(int N){
    int i,j;
    CITY *city0, *city1;
    
    distances = malloc(N*sizeof(float *));
    for( i=0; i<N; i++ )
        distances[i] = malloc(N*sizeof(float));

    for( i=0; i<N; i++ ){
        city0 = list_of_cities[i];
        for( j=0; j<N; j++ ){
            city1 = list_of_cities[j];
            if( j<i )
                distances[i][j] = distances[j][i];
            else
                distances[i][j] = calculate_distance(city0, city1);
            // printf("[%d] -> [%d]  = %f\n", i, j, distances[i][j]);
        }
    }
}

void    init_seq(int N, char *file){
    init_cities(N, file);
    init_gens(N);
    shuffle_populations(N);
}

void    shuffle_populations(int N){
    int i, p, index1, index2;
    int tmp;

    for( p=0; p<MAX_POPUL; p++ ){
        for( i=0; i<MAX_SHUFFLE; i++ ){
            index1 = rand()%N;
            index2 = rand()%N;
            if( index1==index2 ){
                index2 = (index1+1)%N;
            }
            tmp = curr_gen[p][index1]->x;
            curr_gen[p][index1]->x = curr_gen[p][index2]->x;
            curr_gen[p][index2]->x = tmp;
            tmp = curr_gen[p][index1]->y;
            curr_gen[p][index1]->y = curr_gen[p][index2]->y;
            curr_gen[p][index2]->y = tmp;
            tmp = curr_gen[p][index1]->ID;
            curr_gen[p][index1]->ID = curr_gen[p][index2]->ID;
            curr_gen[p][index2]->ID = tmp;
        }
        curr_gen[p][N]->x  = curr_gen[p][0]->x;
        curr_gen[p][N]->y  = curr_gen[p][0]->y;
        curr_gen[p][N]->ID = curr_gen[p][0]->ID;
    }
    
}

float   calculate_distance(CITY *chrom1, CITY *chrom2){
    float x0 = (float)(chrom1->x);
    float y0 = (float)(chrom1->y);
    float x1 = (float)(chrom2->x);
    float y1 = (float)(chrom2->y);

    float dist = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1) );
    return dist;
}

float   read_distance(CITY *chrom1, CITY *chrom2){
    int index0 = chrom1->ID;
    int index1 = chrom2->ID;
    return distances[index0][index1];
}

void   calculate_fitnesses(int N){
    int i,j;
    float total_dist;
    CITY *chrom1, *chrom2;
    
    for( i=0; i<MAX_POPUL; i++ ){
        total_dist = 0.0;
        for( j=0; j<N; j++ ){
            chrom1 = curr_gen[i][j];
            chrom2 = curr_gen[i][j+1];
            total_dist += read_distance(chrom1, chrom2);
        }
        if( total_dist==0.0)
            fitness[i] = 1.0/(total_dist + 1.0);
        else
            fitness[i] = 1.0/(total_dist);
    }
}

void    normalize_fitnesses(void){
    int i;
    float sum = 0.0;
    for( i=0; i<MAX_POPUL; i++ ){
        sum += fitness[i];
    }
    for( i=0; i<MAX_POPUL; i++ ){
        norm_fitness[i] = fitness[i] / sum;
    }
}

void    update_best_distance(int N){
    int i;
    best_distance = 0.0;
    for( i=0; i<N; i++)
        best_distance += read_distance(best_offspring[i], best_offspring[i+1]);
}

int    update_best_offspring(int N){
    int i, j;
    int new_best_flag = 0; 
    for( i=0; i<MAX_POPUL; i++ ){
        if( fitness[i] > best_fitness ){
            best_fitness = fitness[i] + 0.000001;
            for( j=0; j<N+1; j++ ){
                best_offspring[j]->ID = curr_gen[i][j]->ID;
                best_offspring[j]->x  = curr_gen[i][j]->x;
                best_offspring[j]->y  = curr_gen[i][j]->y;
            }
            update_best_distance(N);
            new_best_flag = 1;
        }
    }
    if( new_best_flag ){
        gen_count = 0;
        return 1;
    }
    else{
        gen_count++;
        if( gen_count>=MAX_GEN_CNT )
            done_flag = 1;
        return 0;
    }
}

void    pick_fit_parents(int N){
    int i,j,k;
    int index = 0;
    float sum;
    float choice;

/* IMPORTANT! WHEN CHOICE CLOSE TO 1.00000 index GOES OUT OF *
 * RANGE BECAUSE, SINCE WE ARE USING float, THE SUM OF ALL   *
 * NORMALIZED FITNESS DOESN'T ADD UP TO EXACTLY 1.00000      *
 * TO FIX IT EITHER SUBTRACT FROM choice OR USE double.      */

    for( i=0; i<MAX_POPUL; i++ ){
        sum    = norm_fitness[0];
        choice = (float)rand()/(float)(RAND_MAX) - 0.00001;
        // printf("choice = %f\n", choice);
        while( choice > sum){
            index++;
            sum += norm_fitness[index];
        }
        for( j=0; j<N+1; j++ ){
            next_gen[i][j]->ID = curr_gen[index][j]->ID;
            next_gen[i][j]->x  = curr_gen[index][j]->x;
            next_gen[i][j]->y  = curr_gen[index][j]->y;
        }
        index = 0;
    }
}

void    copy_next_gen_to_curr(int N){
    int i,j;
    for( i=0; i<MAX_POPUL; i++ ){
        for( j=0; j<N+1; j++ ){
            curr_gen[i][j]->ID = next_gen[i][j]->ID;
            curr_gen[i][j]->x  = next_gen[i][j]->x;
            curr_gen[i][j]->y  = next_gen[i][j]->y;
        }
    }
}

void    choose_mutation_indeces(int N, int *index1, int *index2){
    *index1 = rand()%N;
    if( (*index1)>=N/2 ){
        *index2 = (N/2) + (rand()%(N/2));
        if( (*index1)==(*index2) ){
            *index2 = ((*index1)+1)%(N/2) + (N/2);
        }
    }
    else{
        *index2 =(rand()%(N/2));
        if( (*index1)==(*index2) ){
            *index2 = ((*index1)+1)%(N/2);
        }
    }
}

void    swap_chroms(int N, int offspring, int index1, int index2){
    int tmp;
    tmp = next_gen[offspring][index1]->ID;
    next_gen[offspring][index1]->ID = next_gen[offspring][index2]->ID;
    next_gen[offspring][index2]->ID = tmp;
    tmp = next_gen[offspring][index1]->x;
    next_gen[offspring][index1]->x = next_gen[offspring][index2]->x;
    next_gen[offspring][index2]->x = tmp;
    tmp = next_gen[offspring][index1]->y;
    next_gen[offspring][index1]->y = next_gen[offspring][index2]->y;
    next_gen[offspring][index2]->y = tmp;

    next_gen[offspring][N]->x  = next_gen[offspring][0]->x;
    next_gen[offspring][N]->y  = next_gen[offspring][0]->y;
    next_gen[offspring][N]->ID = next_gen[offspring][0]->ID;
}

void    mutate_half_population(int N){
    int i,j, tmp;
    int index1, index2;
    
    for( i=0; i<MAX_POPUL/2; i++ ){
        for( j=0; j<N; j++ ){
            if( rand()%100 >MUTATION_RATE ){
                choose_mutation_indeces(N, &index1, &index2);
                // printf("BEFORE:\tnext_gen[%d] = {", i);
                // for( j=0; j<N; j++ )
                //     printf("%d,", next_gen[i][j]->ID);
                // printf("} index1=%d, index2 = %d\n", index1, index2);
                swap_chroms(N, i, index1, index2);
                // printf("AFTER:\tnext_gen[%d] = {", i);
                // for( j=0; j<N; j++ )
                //     printf("%d,", next_gen[i][j]->ID);
                // printf("}\n");
            }
        }
    }
}

void    mutate_best_offspring(int N){
    int i, j;
    int index1, index2;
    for( i=(MAX_POPUL-100); i<MAX_POPUL; i++ ){
        for( j=0; j<N; j++ ){
            if( rand()%100 >MUTATION_RATE ){
                choose_mutation_indeces(N, &index1, &index2);
                next_gen[i][index1]->ID = best_offspring[index2]->ID;
                next_gen[i][index2]->ID = best_offspring[index1]->ID;
                next_gen[i][index1]->x  = best_offspring[index2]->x;
                next_gen[i][index2]->x  = best_offspring[index1]->x;
                next_gen[i][index1]->y  = best_offspring[index2]->y;
                next_gen[i][index2]->y  = best_offspring[index1]->y;
                for( j=0; j<N; j++ ){
                    if( (j!=index1) && (j!=index2) ){
                        next_gen[i][j]->ID = best_offspring[j]->ID;
                        next_gen[i][j]->x  = best_offspring[j]->x;
                        next_gen[i][j]->y  = best_offspring[j]->y;
                    }
                }
                next_gen[i][N]->x  = next_gen[i][0]->x;
                next_gen[i][N]->y  = next_gen[i][0]->y;
                next_gen[i][N]->ID = next_gen[i][0]->ID;
            }
        }
        
    }
}

void    generate_crossover_range(int *index1, int *index2, int N){
    *index1 = rand()%(N-1);
    *index2 = *index1 + rand()%(N - *index1);
    if( (*index1==0) && (*index2>N-3) )
        *index1 = 2;
    else if( *index1==*index2 )
        *index2 += 1;
}

void    create_crossover_offspring(int p1_index, int p2_index, int ng_index, int coi1, int coi2, int N){
    int entry_index, i,j;
    /* Copy elements in crossover range */
    for( i=coi1; i<=coi2; i++){
        next_gen[ng_index][i]->ID = curr_gen[p2_index][i]->ID;
        next_gen[ng_index][i]->x  = curr_gen[p2_index][i]->x;
        next_gen[ng_index][i]->y  = curr_gen[p2_index][i]->y;
    }
    /* Fill the rest of the array with elements for parent_1 */
    if( coi1==0 )
            entry_index = coi2+1;
        else
            entry_index = 0;

    for( i=0; i<N; i++){
        for( j=coi1; j<=coi2; j++ ){
            if( curr_gen[p1_index][i]->ID==curr_gen[p2_index][j]->ID )
                break;
        }
        if( j==coi2+1 ){
            next_gen[ng_index][entry_index]->ID = curr_gen[p1_index][i]->ID;
            next_gen[ng_index][entry_index]->x  = curr_gen[p1_index][i]->x;
            next_gen[ng_index][entry_index]->y  = curr_gen[p1_index][i]->y;
            entry_index++;
            if( entry_index==coi1 )
                entry_index = coi2+1;
        }
    }
    next_gen[ng_index][N]->ID = next_gen[ng_index][0]->ID;
    next_gen[ng_index][N]->x  = next_gen[ng_index][0]->x;
    next_gen[ng_index][N]->y  = next_gen[ng_index][0]->y;

}

void    crossover_other_half(int N){
    int i,j,k;
    int index1, index2;
    int parent1, parent2;
    int input_index;
    int city1, city2;

    
    for( i=MAX_POPUL/2; i<(MAX_POPUL-200); i++ ){
        parent2 = rand()%MAX_POPUL;
        if( parent2==i )
            parent2 = i - MAX_POPUL/2;
        generate_crossover_range( &index1, &index2, N );
        create_crossover_offspring( i, parent2, i, index1, index2, N);
    }
}

void    crossover_best_offspr(int N){
    int index1, index2, entry_index, i,j,k;
    int parent_2;
    
    for( i=MAX_POPUL-200; i<MAX_POPUL-100; i++ ){
        generate_crossover_range( &index1, &index2, N );
        parent_2 = rand()%(MAX_POPUL-200);
        for( j=index1; j<=index2; j++){
            next_gen[i][j]->ID = best_offspring[j]->ID;
            next_gen[i][j]->x  = best_offspring[j]->x;
            next_gen[i][j]->y  = best_offspring[j]->y;
        }
        /* Fill the rest of the array with elements for parent_1 */
        if( index1==0 )
                entry_index = index2+1;
            else
                entry_index = 0;

        for( k=0; k<N; k++){
            for( j=index1; j<=index2; j++ ){
                if( best_offspring[k]->ID==curr_gen[parent_2][j]->ID )
                    break;
            }
            if( j==index2+1 ){
                next_gen[i][entry_index]->ID = curr_gen[parent_2][k]->ID;
                next_gen[i][entry_index]->x  = curr_gen[parent_2][k]->x;
                next_gen[i][entry_index]->y  = curr_gen[parent_2][k]->y;
                entry_index++;
                if( entry_index==index1 )
                    entry_index = index2+1;
            }
        }
        next_gen[i][N]->ID = next_gen[i][0]->ID;
        next_gen[i][N]->x  = next_gen[i][0]->x;
        next_gen[i][N]->y  = next_gen[i][0]->y;
    }
}

void    create_new_generation(int N){
    pick_fit_parents(N);
    mutate_half_population(N);
    mutate_best_offspring(N);
    crossover_other_half(N);
    copy_next_gen_to_curr(N);
}




/********************************************** DEBUG functs ************************************************/


void    print_cities(int N){
    int i;
    for( i=0; i<N; i++ ){
        printf("%d. (%d,%d)\n", list_of_cities[i]->ID, list_of_cities[i]->x, list_of_cities[i]->y);
    }
}

void    print_pops_and_fits(int N){
    int i, j;
    for( i=0; i<MAX_POPUL; i++ ){
        for( j=0; j<N; j++ ){
            printf("%d.(%d,%d)", curr_gen[i][j]->ID, curr_gen[i][j]->x, curr_gen[i][j]->y);
            if( j!=(N-1) )
                printf("->");
        }
        printf("\tFit = %f\n", fitness[i]);
    }
}

void    print_best_offspring(int N){
    int i;
    printf("[MASTER]: Found new best:\n\t\t");
    for( i=0; i<(N+1); i++ ){
        printf("%d", best_offspring[i]->ID);
        if( i!=N )
            printf(">");
    }
    printf("\n\t  Fitness = %f\n", best_fitness);
}

void    print_total_fitness(void){
    int i;
    float sum = 0.0;
    for( i=0; i<MAX_POPUL; i++ ){
        sum += norm_fitness[i];
    }
    printf("Total fitness = %f\n", sum);
}

void    print_stats(int N, int exec_time){
    int i;
    // FILE *fp = fopen("log.txt", "a");
    printf("\n\n*************************************************\n");
    printf("*\tExecution time    :\t%d (sec)\t*\n", exec_time/1000);
    printf("*\tTotal generations :\t%lu\t\t*\n", total_generations);
    printf("*\tMaster gen count  :\t%lu\t\t*\n", average_generations);
    printf("*\tOptimal path      :\t\t\t*\n*\t\t\t\t\t\t*\n");
    for( i=0; i<N+1; i++ ){
        printf("*\t\t%d\t", best_offspring[i]->ID);
        if( i!=N )
            printf(">\t\t\t*\n");
    }
    printf("\t\t\t*\n*\t\t\t\t\t\t*\n");
    printf("*\tTotal distance    :\t%f\t*\n", best_distance);
    printf("*\tTotal fitness     :\t%f\t*\n", best_fitness);
    printf("*************************************************\n\n");
    // fprintf(fp, "%lu (sec)\t%lu (gens)\t\tFit: %f, Dist: %f\n", exec_time/1000, total_generations, best_fitness, best_distance);
    // fclose(fp);
}

unsigned int get_clock_count(void){
    struct timeval the_time;
    gettimeofday(&the_time, NULL);
    return (unsigned int)( ((the_time.tv_sec - 879191283) * 1000) + (the_time.tv_usec / 1000) );
}