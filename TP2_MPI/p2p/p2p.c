#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char **argv)
{
    int rang, nbprocs, dest=0, source, etiquette = 50;
    MPI_Status statut;
    char message[10000];
    char message_suiv[100];
    char message_prec[10000];

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rang );
    MPI_Comm_size( MPI_COMM_WORLD, &nbprocs );
    printf("Processeur : %d \n", rang);
    if ( rang != 0 ) {
        MPI_Recv(message_prec, 10000, MPI_CHAR, MPI_ANY_SOURCE, etiquette, MPI_COMM_WORLD, &statut );

        sprintf(message_suiv, "Bonjour de la part de P%d!\n" , rang);
        sprintf(message, "%s;%s", message_prec, message_suiv);

        MPI_Send(message, strlen(message)+1, MPI_CHAR, rang-1, etiquette, MPI_COMM_WORLD );
        printf("J'ai reçu : %s \n", message_prec);
    }
    else{
        sprintf(message_suiv, "Bonjour de la part de P%d!\n" , rang);

        MPI_Send(message_suiv, strlen(message)+1, MPI_CHAR, nbprocs-1, etiquette, MPI_COMM_WORLD );
        MPI_Recv(message_prec, 10000, MPI_CHAR, MPI_ANY_SOURCE, etiquette, MPI_COMM_WORLD, &statut );
        
        printf("J'ai reçu : %s \n", message_prec);
    }

    MPI_Finalize();
    
    return EXIT_SUCCESS ;
}