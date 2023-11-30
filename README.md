# MPI_TSP
Approximate solution to the Traveling Salesman Problem (TSP) using the Message Passing Interface (MPI) for parallelization

Included in the test folder is the analytical solution to the TSP (anal_tsp.c), a non-graphical (no allegro 5 visuals)
version of the parallel Genetic Algorithm solution (mpi_tsp.c), the graphical version of the GA (mpi_graph.c) and a sequential
solution to the TSP using the Genetic algorithm library (sequential_tsp.c). How to run:

make all
mpirun -n P ./mpi_tsp N cities_file.txt

P = num of procs
N = number of cities

eg: 
mpirun -n 4 ./mpi_tsp 20 cities.txt

- Master and Slaves have their own table of permutations.
- All perform distance calculations for every permutation
in the population table (they all access the city-to-city
distances through a table they have created individually)
- Mutation and Crossover on the existing permutations cre-
ates the new population.
- If a slave finds a local best route, it sends a message
to the master containing:
i) the route itself
ii) the route's fitness
iii) a generation count increment value (how many generations
passed since the last +100 gen count message).
- It is highly probable to send a follow-up new best messa-
ge if one has just been found, since in the next generation
mutation and crossover takes place on the local best as
well (trying to find small tweaks to the best route to make
it even shorter).
- If a slave doesn't find a local best route for 100 gene-
rations, a message is sent to the master. Upon receiving,
the master increments both the total generation counter
and the termination condition counter by 100.
- The termination condition counter resets if a new total
best has been found by either the master or the slaves.
- If the counter reaches a certain threshold, the program
terminates

Terminal prints:
- Master found new total best
- Slave found new local best
- Master received new local best
- Master received +100 message
