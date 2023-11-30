all: mpi_tsp.c tsp.c
	export PATH=/home/sokratisb/mpi/bin/:$PATH
	mpicc mpi_graph.c tsp.c -o mpi_tsp -lm `pkg-config --cflags --libs allegro-5 allegro_primitives-5`

clean:
	rm mpi_tsp
