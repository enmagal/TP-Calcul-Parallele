#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include "mpi.h"

int main(int argc, char **argv)
{
    int rang=-1, nbprocs=0;
    char processor_name[100];
    int cpu_id, namelen;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rang );
    MPI_Comm_size( MPI_COMM_WORLD, &nbprocs );
    printf( " Hello from process %d of %d\n ", rang, nbprocs);

    MPI_Get_processor_name(processor_name, &namelen);
    printf("Le nom du processeur : %s\n", processor_name);
    cpu_id = sched_getcpu();
    printf(" CPU Id : %d\n", cpu_id);

    MPI_Finalize();

    return EXIT_SUCCESS;
}