#include "pingPong-a.h"

void ping_pong_a(int argc, char **argv)
{
    int myRank, numProcs;
    int pingCount = 42;
    int pongCount = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status  status;
    if (myRank == 0)
    {
        printf("Sending Ping (# %i)\n", pingCount);
        MPI_Send(&pingCount, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        printf("Awaiting pong \n");
        MPI_Recv(&pongCount, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        printf("Received Pong (# %i)\n", pongCount);
    }
    else
    {
        MPI_Recv(&pingCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        pongCount ++;
        pongCount *= -1;
        printf("Sending pong \n");
        MPI_Send(&pingCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        printf("Sending Pong (# %i)\n", pongCount);
        printf("Received Ping (# %i)\n", pingCount);
    }

    MPI_Finalize();
}
