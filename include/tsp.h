#define MAX_POPUL     10000
#define MAX_SHUFFLE   100
#define MAX_GEN_CNT   2000
#define MUTATION_RATE 90

typedef struct city{
    int x;
    int y;
    int ID;
}CITY;

CITY  **list_of_cities;
CITY ***curr_gen;
CITY ***next_gen;
float  *fitness;
float  *norm_fitness;
float   best_fitness;
float   best_distance;
float **distances;
CITY  **best_offspring;
unsigned long  gen_count;
unsigned short done_flag;
unsigned long  total_generations;
unsigned long  average_generations;
unsigned long  last_update_gen_count;

void    init_cities(int N, char *file);
void    init_gens(int N);
void    init_distances(int N);
void    init_seq(int N, char *file);
void    shuffle_populations(int N);
float   calculate_distance(CITY *chrom1, CITY *chrom2);
float   read_distance(CITY *chrom1, CITY *chrom2);
void    calculate_fitnesses(int N);
void    normalize_fitnesses(void);
void    update_best_distance(int N);
int     update_best_offspring(int N);
void    pick_fit_parents(int N);
void    choose_mutation_indeces(int N, int *index1, int *index2);
void    swap_chroms(int N, int offspring, int index1, int index2);
void    mutate_half_population(int N);
void    mutate_best_offspring(int N);
void    generate_crossover_range(int *index1, int *index2, int N);
void    create_crossover_offspring(int p1_index, int p2_index, int ng_index, int coi1, int coi2, int N);
void    crossover_other_half(int N);
void    crossover_best_offspr(int N);
void    create_new_generation(int N);

void    print_cities(int N);
void    print_pops_and_fits(int N);
void    print_best_offspring(int N);
void    print_total_fitness(void);
void    print_stats(int N, int exec_time);

unsigned int get_clock_count(void);